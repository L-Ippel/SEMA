#' Test whether matrix is invertible
#' 
#' @description \code{try_solve} is a function which tries to invert the 
#'   squared fixed effects covariates matrix.
#' @param x_sq A square matrix of fixed effects variables squared.
#' @return If succesful: the inverted matrix. If not succesful: NULL
try_solve <- function(x_sq){
  tryCatch(solve(x_sq),
           error = function(e){
             return(NULL)
             })
}

#' update_x_inv is the online update function of the inverse matrix.
#' @param data_fixed A fixed covariates vector from the new data point.
#' @param x_inv The current inverted matrix.
#' @return The updated inverted matrix.

update_x_inv <- function(data_fixed,
                         x_inv){
  return(
    x_inv - ((x_inv %*% data_fixed %*% t(data_fixed) %*% x_inv) /
               as.numeric((1 + t(data_fixed) %*% x_inv %*% data_fixed))))
}

#' compute_fixed_coef computes the coefficients of the fixed effects,
#'   see from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
#'   EQUATION 14.10.
#' @param x_inv The inverted matrix of the fixed effects variables.
#' @param xy The product of fixed effects variables and dependent variable.
#' @param t1 A vector with Complete Data Sufficient Statistics of the
#'   fixed effects.
#' @return A vector with fixed effects coefficients.

compute_fixed_coef <- function(x_inv,
                               xy,
                               t1){
  return(x_inv %*% (xy - t1)) #estimation of fixed coefficients 1xp
}

#' compute_random_var computes the variance of the random effects
#'   see from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
#'   EQUATION 14.10.
#' @param t2 A matrix with Complete Data Sufficient Statisticsof the
#'   random effects.
#' @param units The total number of units (scalar).
#' @return A matrix with variance of the random effects.

compute_random_var <- function(t2,
                               units){
  return(t2 / units)
}

#' compute_resid_var computes the residual variance.
#' @param t3 The Complete Data Sufficient Statistic of the residual variance,
#'   a scalar.
#' @param n The number of observations.
#' @return The residual variance (scalar).

compute_resid_var <- function(t3,
                              n){
  return(as.numeric(t3 / n))
}
