#' Fit multilevel models in a data stream I
#' 
#' @description Fit a multilevel model online, row-by-row, without storing
#'   data points
#'   
#' @details This is the main function to fit the multilevel models in a data 
#'   stream. The function takes in one observation which consists of the id
#'   number of the unit the data of the fixed effects covariates, the values
#'   of the random effects covariates, the response or outcome, and the current
#'   state of the model parameters. Currently the algorithm can fit models including
#'   fixed effects at level 1 and 2 and random intercepts and slopes for
#'   continuous outcomes. The user manages storage and retrieval of unit's
#'   parameters. This function is also used in \code{\link{sema_fit_set}} and
#'   \code{\link{sema_fit_df}}. 
#' @seealso \code{\link{sema_fit_set}}, \code{\link{sema_fit_df}},
#' \code{\link{summary_sema}}, \code{\link{ranef}},
#'   \code{\link{store_resid_var}}, \code{\link{store_random_var}}, 
#'   \code{\link{store_fixed_coef}}
#'
#' @param data_fixed A vector with the data of the fixed effects covariates.
#' @param data_random A vector with the data of the random effects covariates.
#' @param data_y A scalar with the response of this unit.
#' @param id A scalar which identifies the unit of this data point.
#' @param theta_j A list with this unit's parameters and contributions to the
#'   sufficient statistics.
#' @param theta A list with model parameters and sufficient statistics.
#' @param print The default is FALSE, if TRUE the function
#'   prints a summary of the model.
#' @param start_resid_var A scalar, optional if the user wants to provide a
#'   start value of the residual variance, default start value is 1.
#' @param start_random_var A vector, optional if the user wants to provide a
#'   start values of the variance of the random effects covariates, default
#'   start value is 1. NOTE, if start values are provided make sure that the
#'   length of the vector of start values matches the number of random effects.
#' @param start_fixed_coef A vector, optional if the user wants to provide
#'   start values of the fixed effects, default is set to NULL such that
#'   \code{sema_fit_one} creates the vector of start values matching the number
#'   of fixed effects. NOTE, if start values are provided make sure that the
#'   length of the vector of start values matches the number of fixed effects.
#' @param prior_n A scalar, if starting values are provided, prior_n determines
#'   the weight of the starting value of the residual variance, default is 0.
#' @param prior_j A scalar, if starting values are provided, prior_j determines
#'   the weight of the starting value of the variance of the random effects and
#'   the fixed effects, default is 0.
#' @keywords online multilevel models method fitting stream
#' @export
#' @return A list with a list with updated unit level parameters for one unit
#'   and a list with updated global parameters.

sema_fit_one <- function(data_fixed,
                         data_random,
                         data_y,
                         id,
                         theta_j,
                         theta,
                         print = FALSE,
                         start_resid_var  = 1,
                         start_random_var = 1,
                         start_fixed_coef = NULL,
                         prior_n = 0,
                         prior_j = 0){
  if(!is.list(theta)){
    if(is.null(start_fixed_coef)){
      start_fixed_coef <- rep(1, length(data_fixed))
    }
  }

  theta <- try_theta(theta            = theta,
                     n_fixed          = length(data_fixed),
                     n_random         = length(data_random),
                     start_resid_var  = start_resid_var,
                     start_random_var = start_random_var,
                     start_fixed_coef = start_fixed_coef,
                     prior_n          = prior_n,
                     prior_j          = prior_j)

  theta_j <- try_theta_j(theta_j  = theta_j,
                         n_fixed  = length(data_fixed),
                         n_random = length(data_random),
                         ids      = id)

  theta_j	<- update_id(parameters_id = theta_j,
                       data_fixed    = data_fixed,
                       data_random   = data_random,
                       data_y        = data_y)

  theta$n         <- theta$n + 1
  if(theta_j$n_j == 1){
    theta$j <- theta$j + 1
  }
  theta$xy_vector <- theta$xy_vector + data_fixed * data_y
  theta$y	        <- update_average(old = theta$y, obs = data_y, n = theta$n)

  theta_j$c_inv	<- compute_c_inv(z_sq       = theta_j$z_sq,
                                 resid_var  = theta$resid_var_hat,
                                 random_var = theta$random_var_hat)

  theta_j$random_coef	<- compute_random_coef(c_inv      = theta_j$c_inv,
                                             fixed_coef = theta$fixed_coef_hat,
                                             zy         = theta_j$zy,
                                             zx         = theta_j$zx_mat)

  temp 		      <- online_e_step(theta_main = theta, theta_id = theta_j)
  theta$t1	    <- temp$t1
  theta_j$t1_j	<- temp$t1_j
  theta$t2	    <- temp$t2
  theta_j$t2_j	<- temp$t2_j
  theta$t3  	  <- temp$t3
  theta_j$t3_j	<- temp$t3_j

  if(is.null(theta$x_inv)){
    theta$x_sq	  <- update_product(old   = theta$x_sq,
                                   part_a = data_fixed,
                                   part_b = t(data_fixed))
    theta$x_inv	  <- try_solve(x_sq = theta$x_sq)
  }
  else{
    theta$x_inv	          <- update_x_inv(data_fixed  = data_fixed,
                                          x_inv       = theta$x_inv)
    theta$fixed_coef_hat	<- compute_fixed_coef(x_inv = theta$x_inv,
                                                xy    = theta$xy_vector,
                                                t1    = theta$t1)
    theta$random_var_hat	<- compute_random_var(t2    = theta$t2,
                                                units = theta$j)
    theta$resid_var_hat	  <- compute_resid_var(t3     = theta$t3,
                                               n      = theta$n)
  }
  if(print){
    print(summary_sema(x = theta))
  }
  final <- list(unit = theta_j, model = theta)
  return(final)
}
