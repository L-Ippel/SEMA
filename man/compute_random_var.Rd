% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mstep.R
\name{compute_random_var}
\alias{compute_random_var}
\title{compute_random_var computes the variance of the random effects
  see from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
  EQUATION 14.10.}
\usage{
compute_random_var(t2, units)
}
\arguments{
\item{t2}{A matrix with Complete Data Sufficient Statisticsof the
random effects.}

\item{units}{The total number of units (scalar).}
}
\value{
A matrix with variance of the random effects.
}
\description{
compute_random_var computes the variance of the random effects
  see from raudenbush and bryk (2002) Hierarchial linear models, 2nd edition,
  EQUATION 14.10.
}
