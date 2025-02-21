library(match2C)
library(mvtnorm)
library(tidyverse)
library(dplyr)
library(tableone)
library(tidyr)
library(xtable)
library(exact2x2)
library(sensitivitymult)
library(rbounds)

############################### Helper function #################################################
# Checking and correcting data inconsistencies
debug_contr_error <- function (dat, subset_vec = NULL) {
  if (!is.null(subset_vec)) {
    ## step 0
    if (mode(subset_vec) == "logical") {
      if (length(subset_vec) != nrow(dat)) {
        stop("'logical' `subset_vec` provided but length does not match `nrow(dat)`")
      }
      subset_log_vec <- subset_vec
    } else if (mode(subset_vec) == "numeric") {
      ## check range
      ran <- range(subset_vec)
      if (ran[1] < 1 || ran[2] > nrow(dat)) {
        stop("'numeric' `subset_vec` provided but values are out of bound")
      } else {
        subset_log_vec <- logical(nrow(dat))
        subset_log_vec[as.integer(subset_vec)] <- TRUE
      } 
    } else {
      stop("`subset_vec` must be either 'logical' or 'numeric'")
    }
    dat <- base::subset(dat, subset = subset_log_vec)
  } else {
    ## step 1
    dat <- stats::na.omit(dat)
  }
  if (nrow(dat) == 0L) warning("no complete cases")
  ## step 2
  var_mode <- sapply(dat, mode)
  if (any(var_mode %in% c("complex", "raw"))) stop("complex or raw not allowed!")
  var_class <- sapply(dat, class)
  if (any(var_mode[var_class == "AsIs"] %in% c("logical", "character"))) {
    stop("matrix variables with 'AsIs' class must be 'numeric'")
  }
  ind1 <- which(var_mode %in% c("logical", "character"))
  dat[ind1] <- lapply(dat[ind1], as.factor)
  ## step 3
  fctr <- which(sapply(dat, is.factor))
  if (length(fctr) == 0L) warning("no factor variables to summary")
  ind2 <- if (length(ind1) > 0L) fctr[-ind1] else fctr
  dat[ind2] <- lapply(dat[ind2], base::droplevels.factor)
  ## step 4
  lev <- lapply(dat[fctr], base::levels.default)
  nl <- lengths(lev)
  ## return
  list(nlevels = nl, levels = lev)
}

# McNemar Test for outcome analysis
mcnemar_test <- function(dt_treated, dt_control, str = 'mt30stat'){
  pp=0
  pn=0
  np=0
  nn=0
  treated_outcome = dt_treated[, str]
  control_outcome = dt_control[, str]
  
  # Delete missing outcomes
  ind_missing = which(is.na(treated_outcome) | is.na(control_outcome))
  if (length(ind_missing) > 0) {
    treated_outcome = treated_outcome[-ind_missing]
    control_outcome = control_outcome[-ind_missing]
  }
  n = length(treated_outcome)
  for(i in seq(1,n,1))
  {
    if(treated_outcome[i]==1 & control_outcome[i]==1){pp = pp+1}
    else if(treated_outcome[i]==1 & control_outcome[i]==0){pn = pn+1}
    else if(treated_outcome[i]==0 & control_outcome[i]==1){np = np+1}
    else if(treated_outcome[i]==0 & control_outcome[i]==0){nn = nn+1}
  }
  contingency_mtx <- matrix(c(pp,pn,np,nn),ncol = 2,byrow = T, 
                            dimnames = list(c("t_1","t_0"),c("c_1","c_0")))
  print(contingency_mtx)
  res = exact2x2::mcnemar.exact(contingency_mtx)
  return(c(100*mean(treated_outcome),
           100*mean(control_outcome),
           res$estimate, res$conf.int, res$p.value))
}

################ Data processing ####################
# Read the data
dat = read.csv("comparison3.csv")
dat[is.na(dat)] = 'NA' # Replace NA with 'NA' string

# Create treatment indicator 'Z' for high-volume hospital and binary mortality outcome
dat = dat %>%
  mutate(Z = transfer_high,
         mt30stat = ifelse(mt30stat == 'Dead', 1, 0)) 

# Filter for relevant surgery years and hospital classes
dat = dat %>%
  filter(surgyear == "2019"|surgyear == "2018") %>%
  filter(((pennclass == 'Class A') & (Z == 0)) | (pennclass == 'Class B')) 

# Select time-unvarying covariates
dat_X = dat %>%
  dplyr::select(age, hdef, hdef_missing, bmi, gender, hypertn, race, cancer, status, heartfailacute) %>%
  mutate(hdef_missing = ifelse(hdef_missing == 0, "No", "Yes"))

# Check which categorical variables have multiple values. 
# If a variable has only one level, then exclude it
nlevel_vars = debug_contr_error(dat_X)
cat_vars_to_include = names(nlevel_vars$nlevels[which(nlevel_vars$nlevels >= 2)])

# Select continuous variables and categorical variables with multiple levels
vars_cont = c('age', 'bmi', 'hdef') 
vars_to_include = c(vars_cont, cat_vars_to_include) 
dat_X = dat_X[, vars_to_include]

# Create model matrix
modelmatrix = model.matrix(as.formula(paste("~", paste(vars_to_include, collapse="+"))),
                           data = dat_X)
modelmatrix = cbind(modelmatrix, Z = dat$Z)

dat$pennclass = ifelse(dat$pennclass == "Class A", 1, 0)

########################## Standard Matching ################################

# Estimate the propensity score using logistic regression for treatment assignment
pscore = glm(Z ~., family = binomial(link = "logit"), 
             data = as.data.frame(modelmatrix[,-1]))$fitted.values
dat$propensity = pscore

# Prepare the covariate matrix for matching
X_match_left = as.data.frame(modelmatrix[,-1]) # Remove intercept
X_match_left = X_match_left[,!duplicated(as.list(X_match_left))] # Remove duplicated columns
X_match_left = as.matrix(X_match_left) # Convert to matrix
X_match_left = X_match_left[, !colnames(X_match_left) %in% c("Z")] # Exclude treatment variable 

# Closely matching using robust Mahalanobis distance on the left side
dist_list_left <- create_list_from_scratch(dat$Z, 
                                           X_match_left, 
                                           p = pscore,
                                           caliper_low = 0.5,
                                           k = 1000,
                                           method = 'robust maha')

# Scale distances for further matching processing
dist_list_left$d <- 1e2 * dist_list_left$d

# Perform matching using propensity scores on the right side
dist_list_right = create_list_from_scratch(Z = dat$Z, 
                                           X = pscore,
                                           p = pscore,
                                           caliper_low = 1,
                                           k = 1000,
                                           method = 'L1') #L1 distance RHS

# Perform near fine balance matching for three categorical variables
vars_fb3 = c('distearloc', 'disretloc', 'distalextloc')
X_fb3 = as.matrix(dat[, vars_fb3])
dist_list_right_fb3 = create_list_from_scratch(Z = dat$Z, 
                                               X = X_fb3,
                                               p = pscore,
                                               caliper_low = 1,
                                               k = 1000,
                                               method = 'L1')

################# Directional Penalty on a Single Time-Varying Covariate: Penn Class #################

# Custom distance function
pennclass_dir_dist <- function(X_control, X_treated_i){
  d_1 = sweep(X_control, 2, as.matrix(X_treated_i)) 
  return(d_1)
}

# Distance between treated and control units based on Penn Class
dist_list_right_dir_pennclass = create_list_from_scratch(Z = dat$Z, 
                                                         X =  as.matrix(dat[, "pennclass"]),
                                                         p = pscore,
                                                         caliper_low = 1,
                                                         k = 1000,
                                                         method = 'other',
                                                         dist_func = pennclass_dir_dist)

# Function for matching with direction penalty 
transfer_match <- function(epsilon_pennclass_LB = 0.05, epsilon_pennclass_UB = 0.1,
                           lambda_ps = 10, a = 5e2, b = 3000, c = 5000){
  
  # Calculate the directional penalties for lower and upper bounds
  d_LB = epsilon_pennclass_LB - dist_list_right_dir_pennclass$d
  d_UB = dist_list_right_dir_pennclass$d - epsilon_pennclass_UB
  
  # Update the matching distances with penalties
  dist_list_right$d = 10*(lambda_ps*dist_list_right$d + a*dist_list_right_fb3$d) + b*d_LB + c*d_UB
  
  # Adjust negative distances (if any) by shifting them to ensure all distances are positive
  if (any(dist_list_right$d < 0)){
    offset_d = -min(dist_list_right$d)
    dist_list_right$d = dist_list_right$d + offset_d
  }
  
  start_time <- Sys.time()
  
  # Pair matching using the updated distance matrices 
  matching_output = match_2C_list(Z = dat$Z, 
                                  dataset = dat,
                                  dist_list_1 = dist_list_left, # close matching 
                                  dist_list_2 = dist_list_right, # fine balancing
                                  lambda = 1,
                                  controls = 1)
  end_time <- Sys.time()
  end_time - start_time # Return the computation time
  
  # Filter the matched pairs and sort by matched set
  dt_matched = matching_output$data_with_matched_set_ind %>% filter(!is.na(matched_set))
  dt_treated = dt_matched %>% filter(Z == 1)
  dt_treated = dt_treated[order(dt_treated$matched_set),]
  dt_control = dt_matched %>% filter(Z == 0)
  dt_control = dt_control[order(dt_control$matched_set),]
  
  # Calculate the proportion of Penn Class A in the control group
  control_ClassA_proportion = sum(dt_control$pennclass)/length(dt_control$pennclass)
  
  # Balance table for the matched data
  vars_fb = c("distearloc2", "disretloc2", "distalextloc2") 
  print(xtable(print(tableone::CreateTableOne(vars = c("pennclass",vars_to_include, vars_fb, "propensity","mt30stat"), 
                                        factorVars = c(cat_vars_to_include, vars_fb3),
                                        data = rbind(dt_treated, dt_control), 
                                        includeNA = T, strata = 'Z', test = T))))
  
  # Perform outcome analysis using McNemar's test for matched pairs
  resB = mcnemar_test(dt_treatedB, dt_controlB, str = 'mt30stat')
  
  # Output the McNemar test results
  resB = as.data.frame(resB)
  rownames(resB) = c("100*mean(treated_outcome)",
                     "100*mean(control_outcome)",
                     "estimate", "CI.lower", "CI.higher", "p.value")
  print(xtable::xtable(resB, digits = 5))
  
  # Output the distribution of Penn Class in both treated and control groups
  print(table(dt_treated$pennclass))
  print(table(dt_control$pennclass))
  print(paste0("control_ClassA_proportion = ", control_ClassA_proportion))
}

##### Sensitivity analysis for a single time-varying confounder: Penn Class

transfer_match(epsilon_pennclass_LB = 0.19, epsilon_pennclass_UB = 0.21,
               lambda_ps = 10, a = 5e2, b = 4800, c = 5026)

