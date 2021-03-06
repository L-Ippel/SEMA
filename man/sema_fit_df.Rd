% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semaFitdf2.R
\name{sema_fit_df}
\alias{sema_fit_df}
\title{Fit multilevel models in a data stream III}
\usage{
sema_fit_df(formula, data_frame = data.frame(), intercept = FALSE,
  print_every = NA, store_every = NA, start_resid_var = 1,
  start_random_var = 1, start_fixed_coef = 1:5, start_cor = 0.15,
  update = NULL, train = NULL, threshold = 1e-04, max_iter = 800,
  prior_n = 0, prior_j = 0)
}
\arguments{
\item{formula}{A symbolic representation of the model, formula is used
similar to lme4's \code{lmer}:
\code{response ~ fixed effects + (random effects | grouping variable)}}

\item{data_frame}{A data frame consisting of the variables mentioned in the
formula.}

\item{intercept}{This indicates whether there is a column in data frame with
1's.}

\item{print_every}{Do you want the results printed to the consule? The
default is NA, meaning no printing, if a number is privided the function
prints a summary of the model every 'print_every' data points.}

\item{store_every}{Do you want to store results during the data stream? The
default is NA, i.e., no results are stored, if a number is privided the
function stores the fixed effects, random effects variance and residual
variance in seperate data frames every 'store_every' data points.}

\item{start_resid_var}{This is optional if the user wants to provide a start
value of the residual variance, default start value is 1.}

\item{start_random_var}{This is optional if the user wants to provide a
start values of the variance of the random effects covariates, default
start value is 1. NOTE, if start values are provided make sure that the
length of the vector of start values matches the number of random effects.}

\item{start_fixed_coef}{This is optional if the user wants to provide start
values of the fixed effects, default is set to NULL such that sema_fit_one
can create the vector of start values matching the number of fixed
effects. NOTE, if start values are provided make sure that the length of
the vector of start values matches the number of fixed effects.}

\item{start_cor}{This is a starting value for the correlations between the
random effects.}

\item{update}{The default is NULL, when an integer is provided
\code{\link{sema_update}} is called to do a full update to recompute all
contributions to the complete data suffient statistics.}

\item{train}{The default value is \code{NULL}, meaning that there SEMA is
fit to the data without a training set. When a different value is
provided, this indicate the first number of rows which are used for the
training set. See \code{\link{emAlgorithm}} for a full description of
training SEMA.}

\item{threshold}{In case of a training set, this thresholds determines
when the EM algorithm should terminate. When the parameter estimates
change less than this threshold, EM algorithm terminates.}

\item{max_iter}{In case of a training set, you can fix the number of
iterations of the EM algorithm.}

\item{prior_n}{If starting values are provided, prior_n determines the
weight of the starting value of the residual variance, default is 0.}

\item{prior_j}{If starting values are provided, prior_j determins the weight
of the starting value of the variance of the random effects and the fixed
effects, default is 0.}
}
\value{
A list with updated global parameters (model),
  a list with lists of all units parameters and contributions (unit),
  if store_every is a number 3 data frames \code{fixed_coef_df},
  \code{random_var_df}, \code{resid_var_df}.
}
\description{
Fit multilevel models online on a data set
}
\details{
This function fits the multilevel models online, or row-by-row
  on a data set. Similar to \code{\link{sema_fit_set}} and
  \code{\link{sema_fit_one}} the algorithm updates the model parameters a
  data point at a time. However, instead of these two functions, this
  function fits the multilevel model on a data set and it uses
  \code{formula}.
}
\examples{
## First we create a dataset, consisting of 2500 observations from 20
## units. The fixed effects have the coefficients 1, 2, 3, 4, and 5. The
## variance of the random effects equals 1, 4, and 9. Lastly the
## residual variance equals 4:
test_data <- build_dataset(n = 1500,
                           j = 200,
                           fixed_coef = 1:5,
                           random_coef_sd = 1:3,
                           resid_sd = 2)

## fit a multilevel model:
m1 <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5  | id),
                   data_frame = test_data, intercept = TRUE)
}
\keyword{models}
\keyword{multilevel}
\keyword{online}
