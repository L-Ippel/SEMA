library(testthat)
library(SEMA)

context("data generating")

test_data <- build_dataset(n = 1000, 
                           j = 15, 
                           fixed_coef = 1:5 ,
                           random_coef_sd = c(1,2),
                           resid_sd = 5,
                           n_level_2_var = 2
                           )
test_that("build_dataset creates long data set", {
  expect_is(test_data, "data.frame")
  expect_equal(sum(is.na(test_data)), 0)
})