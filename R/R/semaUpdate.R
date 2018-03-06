#' Update all contributions to Complete Data Sufficient Statistics
#'  
#' @description Allow SEMA to run additional 'sweeps' over all units, during 
#' the data stream. 
#'
#' @details This function can be used in combination with the
#'  \code{\link{sema_fit_one}}. In cases where either units do not return 
#'  in the data stream frequently, or when the data generating process 
#'  gradually changes, additional sweeps or updates where the required objects
#'  are updated for all units, instead of only for the most recently entered 
#'  unit. 
#'
#' @param theta_jList A list with the objects for every known unit in the data.
#' @param theta A list with all model objects. 
#' @keywords online multilevel models method fitting stream
#' @export

sema_update <- function(theta_jList,
                        theta,
                        J = theta$j,
                        n = theta$n){
  theta$t1   <- rep(0, length(theta$t1))
  theta$t2   <- matrix(0, nrow = nrow(theta$t2), ncol = ncol(theta$t2))
  theta$t3   <- 0
  random.inv <- solve(theta$random_var_hat)
  
  for (i in 1:J) {
    theta_jList[[i]]$c_inv       <- solve(theta_jList[[i]]$z_sq +
                                            as.numeric(theta$resid_var_hat) * 
                                            random.inv)
    theta_jList[[i]]$random_coef <- theta_jList[[i]]$c_inv %*%
      (as.matrix(theta_jList[[i]]$zy) -
          t(theta_jList[[i]]$zx_mat) %*%
          as.matrix(theta$fixed_coef_hat))
    
    theta_jList[[i]]$t1_j        <- theta_jList[[i]]$zx_mat %*%
      as.matrix(theta_jList[[i]]$random_coef)
    
    theta_jList[[i]]$t2_j        <- theta_jList[[i]]$random_coef %*%
      t(theta_jList[[i]]$random_coef) +
      as.numeric(theta$resid_var_hat) * theta_jList[[i]]$c_inv
    
    theta_jList[[i]]$t3_j        <- theta_jList[[i]]$y_sq +
      t(theta$fixed_coef_hat) %*% theta_jList[[i]]$x_sq %*% 
      theta$fixed_coef_hat  +
      t(theta_jList[[i]]$random_coef) %*% theta_jList[[i]]$z_sq %*% 
      theta_jList[[i]]$random_coef -
      2 * theta_jList[[i]]$xy %*% theta$fixed_coef_hat  -
      2 * t(theta_jList[[i]]$zy) %*%  theta_jList[[i]]$random_coef +
      2 * t(theta$fixed_coef_hat) %*% theta_jList[[i]]$zx_mat %*%  
      theta_jList[[i]]$random_coef +
      theta$resid_var * sum(diag(theta_jList[[i]]$c_inv %*% 
                                   theta_jList[[i]]$z_sq))
    
    theta$t1 <- theta$t1 + theta_jList[[i]]$t1_j
    theta$t2 <- theta$t2 + theta_jList[[i]]$t2_j
    theta$t3 <- theta$t3 + theta_jList[[i]]$t3_j
    
  }
  theta$fixed_coef_hat <- theta$x_inv %*% (theta$xy_vector - theta$t1)
  theta$random_var_hat <- theta$t2 / J
  theta$resid_var_hat  <- theta$t3 / n
  class(theta)         <- c("list", "sema")
  results              <- list(model = theta, unit = theta_jList) 
  class(results)       <- c("list", "sema")
  return(results)
}