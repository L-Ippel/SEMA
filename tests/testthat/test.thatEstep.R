library(testthat)
library(SEMA)
context("testing E step functions")


z_sq       <- matrix(1:9, nrow = 3)
resid_var  <- 5
random_var <- diag(3)
test_coef  <- matrix(1:4, nrow = 4)
zy         <- matrix(sqrt(diag(z_sq)) / 3, ncol = 3)
zx         <- matrix(1:12, nrow = 4)

test_c_inv <-  compute_c_inv(z_sq, resid_var, random_var) 
test_that("Estep c_inv", {
  expect_is(test_c_inv, "matrix")
  expect_equal(sum(is.na(test_c_inv)), 0)
  expect_length(test_c_inv[1, ], 3)
  expect_length(test_c_inv[, 1], 3)
  
})

test_random_coef <- compute_random_coef(c_inv = test_c_inv, 
                    fixed_coef = test_coef,
                    zy,
                    zx)
test_that("Estep random coefficient", {
  expect_is(test_random_coef, "matrix")
  expect_equal(sum(is.na(test_random_coef)), 0)
  expect_length(test_random_coef[, 1], 3)
})

test_t1 <- compute_t1_j(zx = zx, random_coef = test_random_coef)
test_that("Estep CDSSt1", {
  expect_is(test_t1, "matrix")
  expect_equal(sum(is.na(test_t1)), 0)
  expect_length(test_t1[, 1], 4)
})

test_t2 <- compute_t2_j(random_coef = test_random_coef,
                        resid_var = 5,
                        c_inv  = test_c_inv)
test_that("Estep CDSSt2", {
  expect_is(test_t2, "matrix")
  expect_equal(sum(is.na(test_t2)), 0)
  expect_length(test_t2[, 1], 3)
  expect_length(test_t2[1, ], 3)
  expect_gt((test_t2[1,1]), 0)
  expect_gt((test_t2[2,2]), 0)
  expect_gt((test_t2[3,3]), 0)
})

parameters_j <- list(y_sq = 49,
                     c_inv = test_c_inv,
                     x_sq = 1:4 %*% t(1:4),
                     z_sq = z_sq,
                     xy = matrix((1:4) * 3, nrow = 1),
                     zy = zy, 
                     zx = zx)

test_t3 <- compute_t3_j(parameters_j,
                        random_coef = t(test_random_coef),
                        fixed_coef = test_coef,
                        resid_var = 5)
test_that("Estep CDSSt3", {
  expect_is(test_t3, "numeric")
  expect_equal(sum(is.na(test_t3)), 0)
})
