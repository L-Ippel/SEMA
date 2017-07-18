#' SEMA predict function 
#' 
#' @description This function returns the predicted scores based on the current 
#'   parameter estimates and the data currently entering of this unit.
#'
#' @param theta The current state of model parameters.
#' @param theta_j The current state of the individual parameters.
#' @param data_fixed The current observation fixed effects variables.
#' @param data_random The current observation random effects variables.
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

get_theta <- function(model){
  return(model$model)
}

get_theta_j <- function(model, id){
  return(model$unit[[id]])
}

