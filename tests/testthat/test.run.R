library(testthat)
library(SEMA)

context("running full SEMA algorithm")
set.seed(653)
test_data <- build_dataset(n = 1500, 
                           j = 200, 
                           fixed_coef = 1:5, 
                           random_coef_sd = 1:3, 
                           resid_sd = 2)


check <- sema_fit_df(formula = y ~ 1 + V3 + V4 + V5 + V6 + (1 + V4 + V5  | id), 
                     data_frame = test_data, intercept = TRUE)

test_that("sema_fit_df fits a simple model ",
{
  expect_is(check, "list")
  expect_length(check, 3)
  expect_equal(sum(is.na(check$model)), 0)
  expect_equal(sum(is.na(check$unit)), 0)
})

check2 <- NULL 
data_fixed_var <- c(3:7)
data_random_var <- c(3)
for(i in 1:nrow(test_data)){
  check2 <- sema_fit_set(data_fixed = test_data[i, data_fixed_var],
                         data_random = test_data[i, data_random_var],
                         data_y = test_data$y[i],
                         id = test_data$id[i],
                         theta_list = check2, 
                         print = FALSE)
}

test_that("sema_fit_set fits a simple model ",
          {
            expect_is(check2, "list")
            expect_length(check2, 3)
            expect_equal(sum(is.na(check2$model)), 0)
            expect_equal(sum(is.na(check2$unit)), 0)
          })


id_records		     <- list(NA)
id_vector		       <- c()
check3                <- NULL
print              <- FALSE

for(i in 1:nrow(test_data)){
  id <- test_data$id[i]
  if(!is.element(id, id_vector)){
    id_vector		 <- c(id_vector, id)
    temp_id			 <- which(id_vector == id)
    id_suff_stat <- NULL
  }
  else{
    temp_id		   <- which(id_vector == id)
    id_suff_stat <- id_records[[temp_id]]
  }
  
  check3		<- sema_fit_one(data_fixed = as.numeric(test_data[i, data_fixed_var]),
                       data_random = as.numeric(test_data[i, data_random_var]),
                       data_y      = test_data$y[i],
                       theta       = check3$model,
                       theta_j     = id_suff_stat,
                       id          = test_data$id[i],
                       print       = print)

  id_records[[temp_id]]	<- check3$unit
}

test_that("sema_fit_one fits a simple model ",
          {
            expect_is(check3, "list")
            expect_length(check3, 2)
            expect_equal(sum(is.na(check3$model)), 0)
            expect_equal(sum(is.na(check3$unit)), 0)
          })


