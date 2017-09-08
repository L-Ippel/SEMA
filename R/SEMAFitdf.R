#' Fit multilevel models in a data stream III
#' 
#' @description Fit multilevel models online on a data set
#'
#' @details This function fits the multilevel models online, or row-by-row
#'   on a data set. Similar to \code{\link{sema_fit_set}} and 
#'   \code{\link{sema_fit_one}} the algorithm updates the model parameters a 
#'   data point at a time. However, instead of these two functions, this 
#'   function fits the multilevel model on a data set and it uses 
#'   \code{formula}.
#'   
#' @param formula A symbolic representation of the model, formula is used 
#'   similar to lme4's \code{lmer}: 
#'   \code{response ~ fixed effects + (random effects | grouping variable)} 
#' @param data_frame A data frame consisting of the variables mentioned in the
#'   formula.
#' @param intercept This indicates whether there is a column in data frame with
#'   1's.
#' @param print_every Do you want the results printed to the consule? The
#'   default is NA, meaning no printing, if a number is privided the function
#'   prints a summary of the model every 'print_every' data points.
#' @param store_every Do you want to store results during the data stream? The
#'   default is NA, i.e., no results are stored, if a number is privided the
#'   function stores the fixed effects, random effects variance and residual
#'   variance in seperate data frames every 'store_every' data points.
#' @param start_resid_var This is optional if the user wants to provide a start
#'   value of the residual variance, default start value is 1.
#' @param start_random_var This is optional if the user wants to provide a
#'   start values of the variance of the random effects covariates, default
#'   start value is 1. NOTE, if start values are provided make sure that the
#'   length of the vector of start values matches the number of random effects.
#' @param start_fixed_coef This is optional if the user wants to provide start
#'   values of the fixed effects, default is set to NULL such that sema_fit_one
#'   can create the vector of start values matching the number of fixed
#'   effects. NOTE, if start values are provided make sure that the length of
#'   the vector of start values matches the number of fixed effects.
#' @param prior_n If starting values are provided, prior_n determines the
#'   weight of the starting value of the residual variance, default is 0.
#' @param prior_j If starting values are provided, prior_j determins the weight
#'   of the starting value of the variance of the random effects and the fixed
#'   effects, default is 0.
#' @keywords online multilevel models
#' @export
#' @examples 
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#' test_data <- build_dataset(n = 1500, 
#'                            j = 200, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## fit a multilevel model: 
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5  | id), 
#'                    data_frame = test_data, intercept = TRUE)
#' @return A list with updated global parameters (model),
#'   a list with lists of all units parameters and contributions (unit),
#'   if store_every is a number 3 data frames \code{fixed_coef_df},
#'   \code{random_var_df}, \code{resid_var_df}.

sema_fit_df <- function(formula,
                        data_frame = data.frame(),
                        intercept = FALSE,
                        print_every = NA,
                        store_every = NA,
                        start_resid_var = 1,
                        start_random_var = 1,
                        start_fixed_coef = NULL,
                        prior_n = 0,
                        prior_j = 0){
  if(!is.data.frame(data_frame)){
    stop("data_frame is not a data frame")
  }
  var_vec     <- all.vars(formula)
  n_var_form  <- length(var_vec)
  data_y_vec  <- data_frame[, match(var_vec[1], names(data_frame))]
  data_id_vec <- data_frame[, match(var_vec[n_var_form], names(data_frame))]

  char_form   <- as.character(formula)[3]
  if(stringr::str_count(char_form, stringr::fixed("(")) > 1){
    stop("incorrect formula")
  }
  temp_formula <- stringr::str_split(char_form, stringr::fixed("("), 
                                     simplify = T)
  temp_fixed   <- stringr::str_split(temp_formula[1], stringr::fixed(" + "),
                                     simplify = T)
  temp  <- stringr::str_split(temp_formula[2], stringr::fixed(" + "),
                              simplify = T)
  
  if(ncol(temp)==1){
    temp_random <- stringr::str_split(temp[, ncol(temp)], 
                                                      stringr::fixed(" | "),
                                                      simplify = T)[1]
  }
  if(ncol(temp)!=1){
    temp_random  <- c(temp[1, 1:(ncol(temp) - 1)],
                    stringr::str_split(temp[, ncol(temp)], 
                                       stringr::fixed(" | "),
                                       simplify = T)[1])
  }
  if(temp_random[1] == 1 | temp_fixed[1] == 1){
    data_frame$intercept <- 1
    temp_random[1] <- temp_fixed[1] <- "intercept"
  }
  fixed_covar <- match(temp_fixed[, 1:(ncol(temp_fixed) - 1)],
                       names(data_frame))
  random_covar <- match(temp_random, names(data_frame))
  if(sum(is.na(c(fixed_covar, random_covar))) > 0){
    stop("variables in formula are not present in data frame")
  }
  if(intercept == FALSE){
    data_fixed_df  <- cbind(1, data_frame[, fixed_covar])
    data_random_df <- cbind(1, data_frame[, random_covar])
  }
  if(intercept == TRUE){
    data_fixed_df  <- data_frame[, fixed_covar]
    data_random_df <- as.data.frame(data_frame[, random_covar])
    if(data_fixed_df[1, 1] != 1){
      stop("first column should be intercept, check whether
          intercept == TRUE is correctly specified or change
           order of the variables")
    }
  }
  id_records		     <- list(NA)
  class(id_records)  <- c("list", "sema")
  id_vector		       <- c()
  res                <- NULL
  print              <- FALSE
  if(!is.na(store_every)){
    length_store <- length(seq(store_every, nrow(data_frame), store_every))
    fixed_coef   <- as.data.frame(n = seq(store_every, nrow(data_frame),
                                          store_every),
                                  matrix(NA, ncol = ncol(data_fixed_df) + 1,
                                         nrow =  length_store))
    random_var   <- as.data.frame(n = seq(store_every, nrow(data_frame),
                                          store_every),
                                  matrix(NA, ncol = ncol(data_random_df) + 1,
                                         nrow =  length_store))
    resid_var <-  as.data.frame(n = seq(store_every, nrow(data_frame),
                                        store_every),
                                matrix(NA, ncol = 2, nrow = length_store))
  }
  for(i in 1:nrow(data_frame)){
    id <- data_id_vec[i]
    if(!is.element(id, id_vector)){
      id_vector		 <- c(id_vector, id)
      temp_id			 <- which(id_vector == id)
      id_suff_stat <- NULL
    }
    else{
      temp_id		   <- which(id_vector == id)
      id_suff_stat <- id_records[[temp_id]]
    }
    
    res		<- sema_fit_one(data_fixed = as.numeric(data_fixed_df[i, ]),
                         data_random = as.numeric(data_random_df[i, ]),
                         data_y      = data_y_vec[i],
                         theta       = res$model,
                         theta_j     = id_suff_stat,
                         id          = id,
                         print       = print)
    if( !is.na(print_every) &
      (res$model$n+1) %% print_every == 0) {print <- TRUE }
    
    id_records[[temp_id]]	<- res$unit
    if(!is.na(store_every) & res$model$n %% store_every == 0){
      fixed_coef[res$model$n / store_every, ] <- store_fixed_coef(
        object = res$model)
      random_var[res$model$n / store_every, ] <- store_random_var(
        object = res$model)
      resid_var[res$model$n / store_every, ] <-  store_resid_var(
        object = res$model)
    }
  }
  if(!is.na(store_every)){ 
    names(fixed_coef) <- c("n", temp_fixed[, 1:(ncol(temp_fixed) - 1)])
    names(random_var) <- c("n", temp_random)
    names(resid_var)  <- c("n", "residual_variance")
    
    final <- list("formula"       = formula,
                  "model"         = res$model,
                  "unit"          = id_records,
                  "fixed_coef_df" = fixed_coef,
                  "random_var_df" = random_var,
                  "resid_var_df"  = resid_var)
  } 
  else{
    final <- list("formula"       = formula,
                  "model"         = res$model,
                  "unit"          = id_records)
  }
  class(final) <- c("list","sema")
  return(final)
}
