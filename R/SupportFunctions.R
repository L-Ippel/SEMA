#' In this file all the functions needed by the SEMA algorithm are listed.
#'   NOTE that we use observation to refer to the level 1 data and
#'   individual or unit to refer to the units at level 2.
#'   So, the observations are clustered within individuals,
#'   but other groupings are dealt with the equally.
#'   NOTE 2, these functions are all read in by the SEMA algorithm meaning
#'   that you do not need to worry about transponing certain vectors or
#'   matrices, it has been taken care of in the main function.

#' update_average is the online update of a sample mean
#' @param old The old (current) mean.
#' @param obs The new observation.
#' @param n The total number of observations.
#' @return The updated sample mean.

update_average <- function(old,
                           obs,
                           n){
  return(old + (obs - old) / n)
}

#' update_product is the online update of a matrix product.
#'   NOTE: if you want to update the X'X matrix, part_a = x,
#'   part_b = t(x), where x = rx1 vector make sure the dimensions of the
#'   vectors of part_a and part_b match to compute the matrix product
#' @param old The old (current) matrix.
#' @param part_a A new observation, left or pre multiply part.
#' @param part_b A new observation, right or post multiply part.
#' @return The updated matrix product.

update_product <- function(old,
                           part_a,
                           part_b){
  return(old + part_a %*% part_b)
}

#' update_id is the online update of all unit level objects
#'   which do not depend on model parameters.
#' @param parameters_id A list of current state of the parameters,
#'   @seealso \code{\link{build_theta_j}}.
#' @param data_fixed A vector with data of the fixed effects covariates.
#' @param data_random A vector with data of random effects covariates.
#' @param data_y A scalar with response or dependent variable.
#' @return A list with updated objects of one unit.

update_id <- function(parameters_id,
                      data_fixed,
                      data_random,
                      data_y){
  parameters_id$n_j    <- parameters_id$n_j + 1
  parameters_id$x_sq	 <- update_product(old   = parameters_id$x_sq,
                                        part_a = data_fixed,
                                        part_b = t(data_fixed))
  parameters_id$y_sq	 <- update_product(old   = parameters_id$y_sq,
                                        part_a = data_y,
                                        part_b = data_y)
  parameters_id$xy		 <- update_product(old   = parameters_id$xy,
                                        part_a = data_y,
                                        part_b = data_fixed)
  parameters_id$z_sq	 <- update_product(old   = parameters_id$z_sq,
                                       part_a  = data_random,
                                       part_b  = t(data_random))
  parameters_id$zx_mat <- update_product(old   = parameters_id$zx_mat,
                                        part_a = data_fixed,
                                        part_b = t(data_random))
  parameters_id$zy		 <- update_product(old   = parameters_id$zy,
                                       part_a  = data_y,
                                       part_b  = data_random)
  return(parameters_id)
}
