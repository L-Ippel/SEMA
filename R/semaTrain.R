

genCovarMatrix <-function(cor.random.var, 
                          varX){
  covar.level2 <- matrix(NA, nrow = length(varX), ncol = length(varX))
  for(i in 1:(length(varX))){
    for(j in 1:length(varX)){
      if(i == j){
        covar.level2[i, j] <- varX[i]
      }
      else{
        covar.level2[i, j] <- covar.level2[j, i] <- 
          cor.random.var * sqrt(varX[i]) * sqrt(varX[j])
      }
    }
  }
  return(covar.level2)
}

#' Fit Multilevel model with EM  
#' 
#' @description To provide SEMA with good starting values, it might help to 
#'   use the first part of the data as a training set. 
#'
#' @details This function fits the multilevel models offline to a small part
#'   of the data set. 
#'   
#' @param data A data frame
#' @param id Integer, indicating in which column the grouping variable 
#'   can be found.
#' @param y Integer, indicating in which column the dependent variable 
#'   can be found.
#' @param start.res This is optional if the user wants to provide a start
#'   value of the residual variance, default start value is 1.
#' @param start.random This is optional if the user wants to provide a
#'   start values of the variance of the random effects covariates, default
#'   start value is 1. NOTE, if start values are provided make sure that the
#'   length of the vector of start values matches the number of random effects.
#' @param start.fixed This is optional if the user wants to provide start
#'   values of the fixed effects, if no start values are provided, 0 is used.
#'   NOTE, if start values are provided make sure that the length of
#'   the vector of start values matches the number of fixed effects.
#' @param start.cor Start value for correlation between the random effects. 
#' @param max.iter Integer, maximum number of iterations for EM algorithm. 
#' @param crit.value A threshold or stopping criterion, if the maximum change
#'   in parameter values from one iteration to the next is less than this 
#'   value, the algorithm stops.    
#' @keywords Expectation Maximization algorithm 
#' @export
emAlgorithm <- function(data,
                        id,
                        y, 
                        start.fixed  = c(90, rnorm(13)),
                        start.random = runif(5, 0, 4), 
                        data.random  = c(3, 9:12),
                        start.cor    = .15, 
                        start.res    = 1,
                        max.iter     = 20, 
                        crit.value   = .0001){
  if(!is.data.frame(data)){
    stop("data is not a data frame")
  } 
  names(data)[id] <- "id"
  names(data)[y] <- "y"
  data <- as.data.frame(cbind("id" = data$id, "y" = data$y, data[, -c(id, y)]))
  
  J <- length(unique(data$id))
  if(is.null(start.fixed)){
    start.fixed <- rep(0, (ncol(data)-2))
  }
  if(is.null((start.random))){
    start.random <- rep(1, length(data.random))
  }
  if(is.null(start.res)){
    start.res <- 1
  }
  if(!is.matrix(start.random)){
    start.random <- genCovarMatrix(cor.random.var = start.cor, 
                                   varX           = start.random)
  }
  n        <- nrow(data)
  theta_j  <- list()
  
  Xmat.inv <- solve(t(as.matrix(data[, -c(1,2)])) %*% as.matrix(data[, -c(1, 2)]))
  XY.vec   <- t(as.matrix(data[, -c(1, 2)]))%*% as.matrix(data[, 2])
  id.statistics <- unique(data$id)
  
  for(t in 1:J){
    temp             <- list()
    temp$id          <- id.statistics[t]
    temp$n_j         <- nrow(data[data$id == id.statistics[t], ])
    temp$z_sq        <- matrix(0, nrow = length(data.random), 
                               ncol = length(data.random)) + t(as.matrix(
      data[data$id == id.statistics[t], data.random]))%*%
      as.matrix(data[data$id == id.statistics[t], data.random])
    
    temp$zx_mat      <- matrix(0, nrow = nrow(Xmat.inv), ncol = length(data.random)) +
      as.matrix(t(data[data$id == id.statistics[t], -c(1, 2)])) %*%
      as.matrix(data[data$id == id.statistics[t], data.random])
    temp$y_sq        <- as.matrix(t(data$y[data$id == id.statistics[t]])) %*% 
      as.matrix(data$y[data$id == id.statistics[t]])
    temp$x_sq        <- as.matrix(t(
      data[data$id == id.statistics[t], -c(1, 2)])) %*%
      as.matrix(data[data$id == id.statistics[t], -c(1, 2)])
    temp$xy          <- as.matrix(t(data$y[data$id == id.statistics[t]])) %*% 
      as.matrix(data[data$id == id.statistics[t], -c(1,2)]) 
    
    temp$zy          <- t(as.matrix(data[data$id == id.statistics[t], data.random])) %*%
      as.matrix(data[data$id == id.statistics[t], 2])
    
    
    theta_j[[t]]     <- temp
  }
  for(i in 1:max.iter){
    T1         <- rep(0, nrow(Xmat.inv))
    T2         <- matrix(0, nrow = length(data.random), ncol = length(data.random))
    T3         <- 0
    random.inv <- solve(start.random)
    for(s in 1:J){
      theta_j[[s]]$c_inv <- solve( theta_j[[s]]$z_sq + start.res * random.inv)
      
      theta_j[[s]]$random_coef <- theta_j[[s]]$c_inv %*% 
        (theta_j[[s]]$zy - t(theta_j[[s]]$zx_mat)  %*% as.matrix(start.fixed))
      theta_j[[s]]$t1_j <- theta_j[[s]]$zx_mat %*% theta_j[[s]]$random_coef
      theta_j[[s]]$t2_j <- theta_j[[s]]$random_coef %*% t(theta_j[[s]]$random_coef) + 
        start.res * theta_j[[s]]$c_inv
      theta_j[[s]]$t3_j <- t((data$y[data$id == id.statistics[s]] - 
                                as.matrix(data[data$id == id.statistics[s], -c(1, 2)]) %*% 
                                start.fixed - 
                                as.matrix(data[data$id == id.statistics[s], data.random]) %*% 
                                theta_j[[s]]$random_coef)) %*% 
        (data$y[data$id == id.statistics[s]] - as.matrix(
          data[data$id == id.statistics[s], -c(1, 2)] ) %*% 
           start.fixed -as.matrix(
             data[data$id == id.statistics[s], data.random]) %*% 
           theta_j[[s]]$random_coef) + start.res * sum(diag(theta_j[[s]]$c_inv  %*% 
                                                              theta_j[[s]]$z_sq))
      T1 <- T1 + theta_j[[s]]$t1_j
      T2 <- T2 + theta_j[[s]]$t2_j
      T3 <- T3 + theta_j[[s]]$t3_j
    }
    fixed  <- Xmat.inv %*% (XY.vec - T1)
    random <- T2 / J
    res    <- as.numeric(T3 / n)
    if(abs(max(fixed - start.fixed))   < crit.value & 
       abs(max(random - start.random)) < crit.value &
       abs(res - start.res)            < crit.value) break
    else{
      start.fixed  <- fixed
      start.random <- random
      start.res    <- res
    }
  }
  theta <- list() 
  theta$fixed_coef_hat  <- fixed
  theta$t1	            <- T1
  
  theta$random_var_hat	<- random
  theta$t2		          <- T2
  
  theta$resid_var_hat	  <- res
  theta$t3              <- T3
  
  theta$n	              <- n							
  theta$j	              <- J						
  
  theta$xy_vector	      <- XY.vec
  theta$x_inv	          <- Xmat.inv
  
  
  return(list(model   = theta, 
              unit    = theta_j, 
              id_list = id.statistics))
}
