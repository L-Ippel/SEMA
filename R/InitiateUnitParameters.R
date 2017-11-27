#' Sema unit list 
#' 
#' This function creates a list with objects required at the unit level 
#' @param n_fixed The number of fixed effects.
#' @param n_random The number of random effects.
#' @param ids The id label.
#' @export
#' @examples 
#' ## We create a list with objects used by the sema_fit functions
#' ## when this function is used outside the sema_fit functions, ids, 
#' ## which is the identifying label must be given 
#'  
#' unit_statistics <- build_theta_j(n_fixed = 5,
#'                                   n_random = 3,
#'                                   ids = 1)
#' @return A list containing the following objects per individual:
#'    id label \code{id}; the number of observations of this unit, \code{n_j};
#'    the square matrix for the data of the fixed effects covariates, 
#'    \code{x_sq}; the sum of the squared response, \code{y_sq}; the product of
#'    the fixed effects covariates with the reponse, \code{xy}, a row vector; 
#'    the square matrix for the data of the random effects covariates
#'    \code{z_sq}; the product of fixed effects covariates and random effects 
#'    covariates \code{zx_mat}, dimensions: number of fixed effects x number of 
#'    random effects; the product of random effects covariates and dependent 
#'    variable \code{zy}, a row vector.

build_theta_j <- function(n_fixed,
                          n_random,
                          ids){
  theta_j             <- list()
  class(theta_j) <- c("list", "sema")
  theta_j$id          <- ids
  theta_j$n_j         <- 0
  theta_j$random_coef <- rep(0, n_random)
  theta_j$z_sq        <- matrix(0, nrow  = n_random, ncol = n_random)
  theta_j$zx_mat      <- matrix(0, nrow  = n_fixed, ncol  = n_random)
  theta_j$c_inv       <- matrix(NA, nrow = n_random, ncol = n_random)
  theta_j$y_sq        <- 0
  theta_j$x_sq        <- matrix(0, n_fixed, n_fixed)
  theta_j$xy          <- rep(0, n_fixed)
  theta_j$zy          <- matrix(0, nrow = n_random)
  theta_j$t1_j        <- rep(0, n_fixed)
  theta_j$t2_j        <- matrix(0, nrow = n_random, ncol = n_random)
  theta_j$t3_j        <- 0
  return(theta_j)
}

#' try_theta_j is a function which creates the unit parameters and objects
#'   related to them. If the input is not a list, this function creates a list
#'   which is also a "sema" class with all needed objects stored within.
#' @param theta_j Either an empty object or an object with class list and sema.
#' @param n_fixed Number of fixed effects.
#' @param n_random Number of random effects.
#' @param ids The id label
#' @return If theta_j is a list, it returns theta_j as it was supplied, else
#'   a list with needed objects is returned.
try_theta_j <- function(theta_j,
                        n_fixed,
                        n_random,
                        ids){
  if(!is.list(theta_j)){
    return(build_theta_j(n_fixed   = n_fixed,
                         n_random  = n_random,
                         ids       = ids))
  }
  else{
    return(theta_j)
  }
}
