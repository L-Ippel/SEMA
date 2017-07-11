
#' compute_c_inv computes the uncertainty of the random effects coefficients.
#'   NOTE, this function contains a matrix inversion of the variance of the
#'   random effects, when this matrix is large, the SEMA algorithm slows down
#'   from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
#'   EQUATION 14.15b

#' @param z_sq Square matrix of the squared random effects variables.
#' @param resid_var Residual variance, a scalar.
#' @param random_var Variance of the random effects, a matrix.
#' @return Inverted matrix with updated uncertainty of the random effects
#'   coefficients.

compute_c_inv <- function(z_sq,
                          resid_var,
                          random_var){
  return(solve(z_sq + resid_var * solve(random_var)))
}

#' compute_random_coef computes the random effects: intercepts and slopes
#'   from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
#'   EQUATION 14.15a
#'
#' @param c_inv The inverse of the uncertainty matrix,
#' @seealso \code{\link{compute_c_inv}}.
#' @param fixed_coef A vector with coefficients of the fixed effects.
#' @param zy The product of random effects variables and dependent variable,
#'   a vector.
#' @param zx The product of random effects variables and fixed effect 
#'   variables, a matrix.
#' @return A vector with predicted random effects coefficients.

compute_random_coef <- function(c_inv,
                                fixed_coef,
                                zy,
                                zx){
  return(c_inv %*% (t(zy) - t(zx) %*% fixed_coef))
}

#' compute_t1_j computes the individual contribution to the
#'   complete data sufficient statistics of the fixed effects.
#' @param zx The product of the random effects and fixed effects variables,
#'   a matrix.
#' @param random_coef A vector of random effects coefficients,
#' @seealso  \code{\link{compute_random_coef}}.
#' @return A vector of individual contributions to the CDSS of the fixed
#'   effects.

compute_t1_j <- function(zx,
                         random_coef){
  return(zx %*% random_coef)
}

#' compute_t2_j computes the individual contribution to the complete data
#'   sufficient statistics of the random effects.
#' @param random_coef A vector of random effects coefficients,
#'   @seealso  \code{\link{compute_random_coef}}.
#' @param resid_var The residual variance, a scalar.
#' @param c_inv The uncertainty matrix of the random effects coefficients,
#'   @seealso  \code{\link{compute_random_coef}}.
#' @return A matrix of individual contributions to the CDSS of variance the
#'   random effects.

compute_t2_j <- function(random_coef,
                         resid_var,
                         c_inv){
  return(random_coef %*% t(random_coef) + resid_var * c_inv)
}

#' compute_t3_j computes the individual contribution to the complete data
#'   sufficient statistics of the residual variance.
#' @param  parameters_j A list with unit parameters:
#'   y_sq, a scalar of #' the sum of the squared responses; x_sq,
#'   a matrix with outer product of the fixed effects covariates;
#'   z_sq, a matrix with outer product of the random effects covariates;
#'   xy, a vector with the product of fixed effects covariates and the response
#'   variable;
#'   zy, a vector with the product of random effects covariates and the response
#'   variable;
#'   xz, product of fixed effects covariates and random effects covariates,
#'   a matrix;
#'   c_inv, the inverse of the uncertainty matrix of the
#'   random effects coefficients,
#'   @seealso \code{\link{compute_c_inv}}.
#' @param random_coef random effects coefficients,
#'   @seealso  \code{\link{compute_random_coef}}.
#' @param fixed_coef A vector of fixed effects coefficients.
#' @param resid_var The residual variance, a scalar.
#' @return A scalar with unit's contribution to the Complete Data
#'   Sufficient Statistics of the residual variance.

compute_t3_j <- function(parameters_j,
                         random_coef,
                         fixed_coef,
                         resid_var){
  fixed_coef_sq_matrix     <- fixed_coef %*% t(fixed_coef)
  fixed_random_coef_matrix <- fixed_coef %*% random_coef
  random_coef_sq_matrix    <- t(random_coef) %*% random_coef
  return(parameters_j$y_sq +
           sum(fixed_coef_sq_matrix * parameters_j$x_sq) +
           sum(random_coef_sq_matrix * parameters_j$z_sq) -
           2 * sum(t(fixed_coef) * parameters_j$xy) -
           2 * sum(random_coef * parameters_j$zy) +
           2 * sum(parameters_j$zx * fixed_random_coef_matrix) +
           resid_var * sum(diag(parameters_j$c_inv %*% parameters_j$z_sq)))
}

#' online_e_step includes the latter three functions, for short hand notation
#'   this function retrieves all the needed information from the global storage
#'   and individual storage.
#'   See from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
#'   EQUATION 14.23
#' @param theta_main A list current state of model parameters.
#' @param theta_id A list with information of the current entering individual.
#' @return A list with contributions and new CDSS,
#'   these do still need to be split and returned to the right objects.

online_e_step <- function(theta_main,
                          theta_id){
  theta_main$t1	<- theta_main$t1 - theta_id$t1_j
  theta_id$t1_j	<- compute_t1_j(zx = theta_id$zx_mat,
                                random_coef = theta_id$random_coef)
  theta_main$t1	<- theta_main$t1 + theta_id$t1_j
  theta_main$t2	<- theta_main$t2 - theta_id$t2_j
  theta_id$t2_j	<- compute_t2_j(random_coef = theta_id$random_coef,
                                resid_var = theta_main$resid_var_hat,
                                c_inv  = theta_id$c_inv)
  theta_main$t2	<- theta_main$t2 + theta_id$t2_j

  theta_main$t3	<- theta_main$t3 - theta_id$t3_j
  theta_id$t3_j	<- compute_t3_j(parameters_j = theta_id,
                                random_coef  = t(theta_id$random_coef),
                                fixed_coef   = theta_main$fixed_coef_hat,
                                resid_var    = theta_main$resid_var_hat)
  theta_main$t3	<- theta_main$t3 + theta_id$t3_j

  return(list(t1   = theta_main$t1,
              t1_j = theta_id$t1_j,
              t2   = theta_main$t2,
              t2_j = theta_id$t2_j,
              t3   = theta_main$t3,
              t3_j = theta_id$t3_j))
}
