#' sema output 
#' 
#' @description \code{is_sema} checks whether the input is sema output
#' 
#' @param x Sema output which should be a list. 
#' @return LOGICAL
#' @export
is_sema <- function(x){
  inherits(x, "sema")
}

#' Interpreting sema output
#' 
#' @description Returns a list with the current model parameter estimates
#' 
#' @details The output of the sema_fit functions are usually large, and 
#'   difficult to read, lists. In order to interpret the output of the sema_fit
#'   functions, \code{summary_sema} returns a general overview of the model 
#'   parameters, the fixed effects coefficients, the random effects variances
#'   and covariances, and the residual variance. The \code{summary_sema} 
#'   function is also incorpprated within the \code{sema_fit_one} function, 
#'   such that with the argument \code{print_every} the user can see a summary
#'   of the updated model parameters every X data points. 
#' @param x A sema model output.
#' @keywords summary model sema
#' @seealso \code{\link{sema_store_fixed_coef}}, 
#'    \code{\link{sema_store_random_var}}, \code{\link{sema_store_resid_var}}, 
#'    \code{\link{sema_ranef}}
#' @return A list with sample size, number of units, the coefficients of the
#'   fixed effects, the variance of the random effects and the residual
#'   variance.
#' @export
summary_sema <- function(x){
  if(is.null(x$model) & !is.null(x$n)){
    x$model <- x
  }
  if(!is_sema(x) & !is_sema(x$model)){
    stop("this object is not sema output")
  }
  sema_overview   <- list()

  sema_overview$sample_size        <- x$model$n
  sema_overview$number_of_units    <- x$model$j

  sema_overview$fixed_coef         <- as.numeric(x$model$fixed_coef_hat)
  sema_overview$var_random_effects <- diag(x$model$random_var_hat)
  n_random                         <- length( sema_overview$var_random_effects)

  x$model$random_var_hat[upper.tri(x$model$random_var_hat, diag = T)] <- NA

  sema_overview$cov_random_effects <- as.table(x$model$random_var_hat)
  dimnames(sema_overview$cov_random_effects) <- list(1:n_random, 1:n_random)

  sema_overview$residual_var       <- x$model$resid_var_hat

  cat('==============================================================', '\n')
  return(sema_overview)
}

#' Extract the random effects coefficients. 
#' 
#' @description This function loops over the list with unit parameters and 
#'   returns a dataframe with the random coefficients, centered around zero.
#' @details The \code{sema_ranef} function works similar to lme4's 
#'   \code{ranef}. The function uses \code{ldply} function from the 
#'   \code{plyr} package to extract all coefficients. 
#' @param x The sema output, a list with unit parameters.
#' @export
#' @keywords coefficients model 
#' @return A data frame with random effect coefficients per unit, sorted by
#'   entry.
sema_ranef <- function(x){
  if(!is_sema(x)){
    stop("this object is not sema output")
  }
  ranef_df <- plyr::ldply(.data = x, .fun = function(x){
    return(as.numeric(rbind(x$id, x$random_coef)))
  })
  names(ranef_df) <- c("id", 1:(ncol(ranef_df) - 1))
  return(ranef_df)
}

#' Store fixed effects coefficients
#' 
#' @description This function extracts the fixed effects coefficients and saves
#'   them in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{sema_ranef}}, \code{\link{sema_store_random_var}}, 
#'   \code{\link{sema_store_resid_var}}
#' @param sema_model The sema output, a list containing the model parameters.
#' @keywords save coefficients
#' @export
#' @return A data frame with the number of observations and the fixed effects
#'   coefficients
sema_store_fixed_coef <- function(sema_model){
  sema_coef_fixed <- matrix(c(sema_model$n,
                              as.numeric(sema_model$fixed_coef_hat)),
                            nrow = 1)
  return(sema_coef_fixed)
}

#' Store the variance of the random effects
#' 
#' @description This function extracts the variance of the random effects and 
#' saves them in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{sema_ranef}}, \code{\link{sema_store_fixed_coef}}, 
#'   \code{\link{sema_store_resid_var}}
#' @param sema_model The sema output, a list containing the model parameters.
#' @export
#' @keywords save coefficients
#' @return A data frame with the number of observations and the random effects
#'   variance.
sema_store_random_var <- function(sema_model){
  sema_coef_random <- as.data.frame(
                        matrix(c(sema_model$n,
                                 as.numeric(diag(sema_model$random_var_hat))),
                        nrow = 1))
  return(sema_coef_random)
}

#' store residual variance
#' 
#' @description This function extracts the residual variance and saves
#'   it in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{sema_ranef}}, \code{\link{sema_store_fixed_coef}},
#'  \code{\link{sema_store_random_var}} 
#' @param sema_model The sema output, a list containing the model parameters.
#' @export
#' @keywords save coefficients
#' @return A data frame with the number of observations and the residual
#'   variance.
sema_store_resid_var <- function(sema_model){
  sema_resid_var <- as.data.frame(
    matrix(c(sema_model$n,
             as.numeric(sema_model$resid_var_hat)),
           nrow = 1))
  return(sema_resid_var)
}
