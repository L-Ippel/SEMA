#' SEMA predict function 
#' 
#' @description This function returns the predicted scores based on the current 
#'   parameter estimates and the data currently entering of this unit.
#'
#' @param theta The current state of model parameters.
#' @param theta_j The current state of the individual parameters.
#' @param data_fixed The current observation fixed effects variables.
#' @param data_random The current observation random effects variables.
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
#' ## predict for the last row of the dataset:                     
#' predict_unit_outcome <- predict(theta = get_theta(model = m1), 
#'                                theta_j = get_theta_j(model = m1, 
#'                                  id = test_data$id[2000]),
#'                                data_fixed = as.numeric(
#'                                  test_data[2000, 3:7]),
#'                                data_random = as.numeric(
#'                                  test_data[2000, c(3, 5, 6)]))
#' @export
#' @keywords prediction


predict	<- function(theta = get_theta(),
                         theta_j = get_theta_j(),
                         data_fixed,
                         data_random){
  if(is.null(theta_j)){
    return(theta$y)
  }
  else{
    return(data_fixed %*% theta$fixed_coef_hat +
             data_random %*% theta_j$random_coef)
  }
}

#' get model statistics from Sema object
#' 
#' @description This function returns a list with all model statistics. The
#'   function is for instance used in the predict function. 
#'   
#' @param model sema object
#' @export 
get_theta <- function(model){
  return(model$model)
}


#' get unit statistics from Sema object
#' 
#' @description This function returns a list with all unit statistics. The
#'   function is for instance used in the predict function. 
#'   
#' @param model sema object
#' @param id unit's identifier
#' 
#' @export 
get_theta_j <- function(model, id){
  return(model$unit[[id]])
}

