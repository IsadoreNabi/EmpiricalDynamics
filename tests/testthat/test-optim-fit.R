# =============================================================================
# Tests for the optim_fit class
#
# `fit_with_optim()` is the last fallback of `fit_specified_equation()`. The
# object it returns used to have no methods at all, so the moment that fallback
# was taken every prediction died with "no applicable method for 'predict'".
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

make_data <- function(n = 60) {
  set.seed(404)
  d <- data.frame(x = seq(-3, 3, length.out = n))
  d$y <- 3.7 * d$x^2 + 5 + stats::rnorm(n, sd = 0.5)
  d
}

test_that("the general optimiser produces an object with methods", {
  d <- make_data()

  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  expect_s3_class(eq$fit, "optim_fit")

  # The exact inverse of the defect: the class now has methods.
  methods_present <- vapply(
    c("predict", "coef", "fitted", "residuals", "print"),
    function(generic) !is.null(utils::getS3method(generic, "optim_fit",
                                                  optional = TRUE)),
    logical(1))
  expect_true(all(methods_present))
})

test_that("predict without newdata returns the fitted values", {
  d <- make_data()
  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  expect_equal(stats::predict(eq$fit), stats::fitted(eq$fit))
  expect_length(stats::fitted(eq$fit), nrow(d))
})

test_that("predict on new rows evaluates the equation there", {
  d <- make_data()
  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  held_out <- d[1:10, , drop = FALSE]
  got <- stats::predict(eq$fit, newdata = held_out)

  expect_length(got, 10L)
  expect_true(all(is.finite(got)))
  expect_equal(got, stats::fitted(eq$fit)[1:10])
})

test_that("coefficients are named and residuals are what they claim to be", {
  d <- make_data()
  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  coefs <- stats::coef(eq$fit)
  expect_named(coefs, c("a", "b"))
  # The truth is a = 3.7, b = 5; the optimiser should be near it.
  expect_equal(unname(coefs), c(3.7, 5), tolerance = 0.2)

  expect_equal(stats::residuals(eq$fit), d$y - stats::fitted(eq$fit))
})

test_that("a constant equation predicts a value per row rather than one value", {
  d <- make_data()

  eq <- fit_specified_equation("a", data = d, response = "y",
                               start = list(a = 1), method = "optim")

  expect_length(stats::fitted(eq$fit), nrow(d))
  expect_length(stats::predict(eq$fit, newdata = d[1:5, , drop = FALSE]), 5L)
})

test_that("printing an optim_fit works", {
  d <- make_data()
  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  expect_output(print(eq$fit), "optim_fit")
  expect_output(print(eq$fit), "Coefficients")
})

test_that("predicting through the symbolic_equation wrapper reaches the method", {
  d <- make_data()
  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1), method = "optim")

  got <- stats::predict(eq, newdata = d[1:5, , drop = FALSE])
  expect_length(got, 5L)
  expect_true(all(is.finite(got)))
})

test_that("Levenberg-Marquardt is not silently abandoned when there are no weights", {
  skip_if_not_installed("minpack.lm")
  d <- make_data()

  # `nlsLM()` decides whether it was given weights with `missing()`, so handing
  # it `weights = NULL` used to break every single fit and drop the package to
  # its fallback. The primary method must now go through without a word.
  expect_silent(
    eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                                 start = list(a = 1, b = 1), method = "LM")
  )
  expect_equal(unname(stats::coef(eq)), c(3.7, 5), tolerance = 0.2)
})

test_that("weights are honoured by the fit", {
  skip_if_not_installed("minpack.lm")
  d <- make_data()
  w <- rep(1, nrow(d))
  w[d$x < 0] <- 100

  unweighted <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                                       start = list(a = 1, b = 1), method = "LM")
  weighted <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                                     start = list(a = 1, b = 1), method = "LM",
                                     weights = w)

  expect_false(isTRUE(all.equal(stats::coef(unweighted), stats::coef(weighted))))
})
