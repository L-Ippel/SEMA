#' sema output 
#' 
#' @description \code{is_sema} checks whether the input is sema output
#' 
#' @param x Sema output which should be a list. 
#' @return LOGICAL
#' @export
is.sema <- function(x){
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
#' @seealso \code{\link{store_fixed_coef}}, 
#'    \code{\link{store_random_var}}, \code{\link{store_resid_var}}, 
#'    \code{\link{ranef}}
#' @return A list with sample size, number of units, the coefficients of the
#'   fixed effects, the variance of the random effects and the residual
#'   variance.
#' @export
#' @examples
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' summary_sema(m1)
summary_sema <- function(x){
  if(is.null(x$model) & !is.null(x$n)){
    x$model <- x
  }
  if(!is.sema(x) & !is.sema(x$model)){
    stop("this object is not sema output")
  }
  sema_overview   <- list()
  if(!is.null(x$formula)){sema_overview$formula <- x$formula}
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
#' @details The \code{ranef} function works similar to lme4's 
#'   \code{ranef}. The function uses \code{ldply} function from the 
#'   \code{plyr} package to extract all coefficients. 
#' @param x The sema output, a list with unit parameters.
#' @export
#' @examples
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' ranef(m1)
#' @keywords coefficients model 
#' @return A data frame with random effect coefficients per unit, sorted by
#'   entry.
ranef <- function(x){
  if(!is.sema(x)){
    stop("this object is not sema output")
  }
  if(!is.null(x$unit)) {x <- x$unit}
  ranef_df <- plyr::ldply(.data = x, .fun = function(x){
    return(as.numeric(rbind(x$id, x$random_coef)))
  })
  names(ranef_df) <- c("id", 1:(ncol(ranef_df) - 1))
  return(ranef_df)
}

#' Extract the fixed effects coefficients. 
#' 
#' @description This function returns the estimated fixed coefficients.
#' @param x The sema output, a list with unit parameters.
#' @export
#' @examples
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' fixef(m1)
#' @keywords coefficients model 
#' @return A data frame with fixed effect coefficients.
fixef <- function(x){
  if(!is.sema(x)){
    stop("this object is not sema output")
  }
  return(data.frame("fixed coefficients" = x$model$fixed_coef_hat))
}
#' Store fixed effects coefficients
#' 
#' @description This function extracts the fixed effects coefficients and saves
#'   them in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{ranef}}, \code{\link{store_random_var}}, 
#'   \code{\link{store_resid_var}}
#' @param object The sema output, a list containing the model parameters.
#' @keywords save coefficients
#' @examples
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' ## to subtract the fixed effects from the m1 object:
#' store_fixed_coef(m1)

#' @export
#' @return A data frame with the number of observations and the fixed effects
#'   coefficients
store_fixed_coef <- function(object){
  sema_coef_fixed <- data.frame(matrix(c(object$model$n,
                              as.numeric(object$model$fixed_coef_hat)),
                            nrow = 1))
  names(sema_coef_fixed) <- c("n", 1:(ncol(sema_coef_fixed)-1))
  return(sema_coef_fixed)
}

#' Store the variance of the random effects
#' 
#' @description This function extracts the variance of the random effects and 
#' saves them in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{ranef}}, \code{\link{fixef}}, 
#'   \code{\link{store_fixed_coef}}, 
#'   \code{\link{store_resid_var}}
#' @param object The sema output, a list containing the model parameters.
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' ## to subtract the variance of the random effects from the m1 object:
#' store_random_var(m1)
#' @export
#' @keywords save coefficients
#' @return A data frame with the number of observations and the random effects
#'   variance.
store_random_var <- function(object){
  sema_coef_random <- data.frame(
                        matrix(c(object$model$n,
                                 as.numeric(diag(object$model$random_var_hat))),
                        nrow = 1))
  names(sema_coef_random) <- c("n", 1:(ncol(sema_coef_random)-1))
  return(sema_coef_random)
}

#' store residual variance
#' 
#' @description This function extracts the residual variance and saves
#'   it in a data frame with the number of data points seen so far. 
#' 
#' @seealso \code{\link{ranef}}, \code{\link{store_fixed_coef}},
#'  \code{\link{store_random_var}} 
#' @param object The sema output, a list containing the model parameters.
#' @export
#' @examples 
#' ## First we create a dataset, consisting of 2500 observations from 20 
#' ## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The 
#' ## variance of the random effects equals 1, 4, and 9. Lastly the 
#' ## residual variance equals 4:
#'   
#' test_data <- build_dataset(n = 2500, 
#'                            j = 20, 
#'                            fixed_coef = 1:5, 
#'                            random_coef_sd = 1:3, 
#'                            resid_sd = 2)
#'                            
#' ## Next, we fit a simple model to these data                           
#' m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5 | id), 
#'                     data_frame = test_data, 
#'                     intercept = TRUE)
#' ## to subtract the residual variance from the m1 object:
#' store_resid_var(m1)
#' @keywords save coefficients
#' @return A data frame with the number of observations and the residual
#'   variance.
store_resid_var <- function(object){
  sema_resid_var <- data.frame(
    matrix(c(object$model$n,
             as.numeric(object$model$resid_var_hat)),
           nrow = 1))
  names(sema_resid_var) <- c("n", "resid-var")
  
  return(sema_resid_var)
}

