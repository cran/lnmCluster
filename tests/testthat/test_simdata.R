context("simdata BIC and loglik")
library(lnmCluster)

W_count <- simdata
range_G <- c(1:2)
range_Q <- c(2:3)


#model selection without permutation
test_that("model selection selects the best model", {
  x <- lnmbiclust(W_count, range_G, range_Q, "GUU", "BIC")
  expect_equal(x$best_model$BIC,-5703.645,tolerance=0.1)
  expect_equal(max(x$best_model$overall_loglik),-2790.956, tolerance = 0.1)
  expect_equal(dim(x$all_fitted_model)[1],4,tolerance=0)
})


test_that("model selection selects the best model in lnmfa", {
  x <- lnmfa(W_count, range_G, range_Q, "GUU", "BIC")

  expect_equal(x$best_model$BIC,-5613.997,tolerance=0.1)
  expect_equal(max(x$best_model$overall_loglik),-2696.332, tolerance = 0.1)
  expect_equal(dim(x$all_fitted_model)[1],4,tolerance=0)
})



#fit one model
test_that("test to fit one model", {
  x <- lnmbiclust(W_count, 2, 2, "GUU", "BIC")

  expect_equal(x$BIC,-5750.454,tolerance=0.1)
  expect_equal(max(x$overall_loglik),-2773.783, tolerance = 0.1)
})


test_that("test to fit one model in lnmfa", {
  x <- lnmfa(W_count, 2, 2, "GUU", "BIC")

  expect_equal(x$BIC,-5613.997,tolerance=0.1)
  expect_equal(max(x$overall_loglik),-2696.332, tolerance = 0.1)
})



#with permutation
test_that("test to fit permutation UUU model", {
  x <- lnmbiclust(W_count, 2, c(2:3), "UUU", "BIC",permutation=TRUE)#should select 33

  expect_equal(x$best_model$BIC,-5757.848,tolerance=0.1)
  expect_equal(max(x$best_model$overall_loglik),-2755.347, tolerance = 0.1)
  expect_equal(dim(x$all_fitted_model)[1],4,tolerance=0)
})
