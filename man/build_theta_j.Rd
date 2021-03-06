% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitiateUnitParameters.R
\name{build_theta_j}
\alias{build_theta_j}
\title{Sema unit list}
\usage{
build_theta_j(n_fixed, n_random, ids)
}
\arguments{
\item{n_fixed}{The number of fixed effects.}

\item{n_random}{The number of random effects.}

\item{ids}{The id label.}
}
\value{
A list containing the following objects per individual:
   id label \code{id}; the number of observations of this unit, \code{n_j};
   the square matrix for the data of the fixed effects covariates, 
   \code{x_sq}; the sum of the squared response, \code{y_sq}; the product of
   the fixed effects covariates with the reponse, \code{xy}, a row vector; 
   the square matrix for the data of the random effects covariates
   \code{z_sq}; the product of fixed effects covariates and random effects 
   covariates \code{zx_mat}, dimensions: number of fixed effects x number of 
   random effects; the product of random effects covariates and dependent 
   variable \code{zy}, a row vector.
}
\description{
This function creates a list with objects required at the unit level
}
\examples{
## We create a list with objects used by the sema_fit functions
## when this function is used outside the sema_fit functions, ids, 
## which is the identifying label must be given 
 
unit_statistics <- build_theta_j(n_fixed = 5,
                                  n_random = 3,
                                  ids = 1)
}
