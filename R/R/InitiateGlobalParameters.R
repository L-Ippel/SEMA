#' Sema model list 
#' 
#' This function creates a list with all required objects at the global level
#'   needed by the sema algorithm 
#'   
#' @details The function creates a list with objects used by the sema 
#'   algorithm, most are not used for interpretation. The user can provide 
#'   start values to enhance the algorithm's performance. In that case, the
#'   function also requires the weight of those start values, 
#'   when no starting values are given default values are used.
#'
#' @param n_fixed Number of fixed effects.
#' @param n_random Number of random effects.
#' @param start_resid_var Start values residual variance, default start
#'   value equals 1.
#' @param start_random_var Start values variance of random effects,
#'   default start values equal 1. Make sure that the length of 
#'   \code{start_random_var} matches the number of random effects. 
#' @param start_fixed_coef The default is set to NULL, when no start values 
#'   are provided, \code{start_fixed_coef} is set to \code{rep(1, n_fixed)} 
#'   within the \code{create_theta_main} function. Note that the length of the
#'   start values must match the number of fixed effects.  
#' @param prior_n A scalar indicating the weight of the start value of the 
#'   residual variance.
#' @param prior_j A scalar indicating the weight of the start values of the  
#'  fixed effects coefficients and random effects variance.
#' @keywords model method start sema fitting
#' @export
#' @examples 
#' ## Create a list of objects required to fit the multilevel model using
#' ## sema:
#' ## NOTE: default start values fixed effect coefficients, residual variance 
#' ## and random effects variance is equal to 1. When this function is used 
#' ## outside the sema_fit functions, prior_n and prior_j, i.e., weight given 
#' ## to the start values should be entered. Within the sema_fit functions 
#' ## defaults are given. 
#' 
#' model_statistics <- create_theta_main(n_fixed = 5,
#'                                       n_random = 3,
#'                                       start_resid_var = 1,
#'                                       start_random_var = 1,
#'                                       start_fixed_coef = NULL,
#'                                       prior_n = 0,
#'                                       prior_j = 0)
#' @return A list which contains: the fixed effects coefficients, 
#'    \code{fixed_coef_hat}; the Complete Data Sufficient Statistic of the 
#'    fixed effects, \code{t1}; the variance of the random effects 
#'    \code{random_var_hat}; the Complete Data Sufficent Statistics of the 
#'    variance of random effects, \code{t2}; the residual variance, 
#'    \code{resid_var_hat}; the Complete Data Sufficiant Statistic of
#'    the residual variance, \code{t3}, the number of observations, \code{n};
#'    the number of individuals, \code{j}; the average dependent variable, 
#'    which can be used for prediction purposes but it is not needed for the
#'    algorithm, \code{y}; the square matrix of squared fixed effects, which is
#'    only updated as long as the matrix is not yet invertible, once inverted,
#'    only the inverted matrix is updated, \code{x_sq}; the product of fixed 
#'    effects and dependent variable, \code{xy_vector}; the inverse of the 
#'    \code{x_sq} matrix, \code{x_inv}.
#'
create_theta_main <- function(n_fixed,
                              n_random,
                              start_resid_var = 1,
                              start_random_var = 1,
                              start_fixed_coef = NULL,
                              prior_n,
                              prior_j){
  theta        <- list()
  class(theta) <- c("list", "sema")
  theta$fixed_coef_hat  <- start_fixed_coef
  theta$t1	            <- start_fixed_coef * prior_j

  theta$random_var_hat	<- matrix(diag(c(start_random_var), n_random),
                                 nrow = n_random,
                                 ncol = n_random)
  theta$t2		          <- theta$random_var_hat * prior_j

  theta$resid_var_hat	<- start_resid_var
  theta$t3            <- start_resid_var * prior_n

  theta$n	<- 0							
  theta$j	<- 0						
  theta$y	<- 0
  theta$x_sq 	    <- matrix(0, nrow = n_fixed, ncol = n_fixed)
  theta$xy_vector	<- rep(0, n_fixed)
  theta$x_inv	    <- NULL
  return(theta)
}

#' try_theta is a function which creates the model parameters and objects
#'   related to them. If the input is not a list, this function creates a list
#'   which is also a "sema" class with all needed objects stored within.
#' @param theta Either an empty object or an object with class list and sema.
#' @param n_fixed Number of fixed effects.
#' @param n_random Number of random effects.
#' @param start_resid_var Starting values residual variance, default starting
#'   value equals 1.
#' @param start_random_var Starting values variance of random effects,
#'   default starting values equal 1.
#' @param start_fixed_coef Starting values fixed effect coefficients,
#'   default starting values equal 1.
#' @param prior_n Weight starting value residual variance.
#' @param prior_j Weight starting values fixed effects coefficients and random
#'   effects variance.
#' @return If theta is a list, it returns theta as it was supplied, else
#'   a list with needed objects is returned.

try_theta <- function(theta,
                      n_fixed,
                      n_random,
                      start_resid_var,
                      start_random_var,
                      start_fixed_coef,
                      prior_n,
                      prior_j){
  if(!is.list(theta)){
    return(create_theta_main(n_fixed          = n_fixed,
                             n_random         = n_random,
                             start_resid_var  = start_resid_var,
                             start_random_var = start_random_var,
                             start_fixed_coef = start_fixed_coef,
                             prior_n          = prior_n,
                             prior_j          = prior_j))
  }
  else{
    return(theta)
  }
}
