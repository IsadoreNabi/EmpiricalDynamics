# =============================================================================
# Regression tests: the fit explore_dynamics() publishes under fit_<pred>
#
# The referent is algebra applied to coefficients recovered from data built
# here, not the package's own output. A cubic in x whose coefficients are
# known is fitted at negligible noise; the published coefficients have to
# reproduce them, in ascending powers, and the published range has to be the
# range of the predictor that was handed in.
#
# The reason the fit is published at all is that a name does not answer what a
# caller needs to ask of a curve. "quadratic" is compatible both with a
# relation that rises over the whole range observed and with one that turns
# inside it, and those are opposite answers about direction. The last test
# here is that distinction, made from the published coefficients alone.
# =============================================================================

test_that("the published fit reproduces the coefficients it was built from", {
  set.seed(20260819)
  x <- seq(-2, 3, length.out = 400)
  y <- 1.5 - 2 * x + 0.75 * x^2 + 0.25 * x^3 + stats::rnorm(400, sd = 1e-4)
  d <- data.frame(x = x, y = y)

  r <- suppressMessages(explore_dynamics(d, target = "y", predictors = "x"))
  f <- r$statistics$fit_x

  expect_false(is.null(f))
  expect_identical(f$form, r$statistics$nonlin_x)
  expect_identical(f$form, "cubic")
  expect_identical(f$degree, 3L)
  expect_length(f$coefficients, 4L)
  expect_equal(f$coefficients, c(1.5, -2, 0.75, 0.25), tolerance = 1e-3)
  expect_equal(f$predictor_range, range(x))
  expect_named(f$aic, c("linear", "quadratic", "cubic"))
  expect_identical(unname(which.min(f$aic)), 3L)
  expect_identical(f$degree, unname(which.min(f$aic)))
})

test_that("the coefficients are in ascending powers, so the curve can be rebuilt", {
  set.seed(20260819)
  x <- stats::runif(400, -1, 4)
  y <- 3 + 0.5 * x - 1.25 * x^2 + stats::rnorm(400, sd = 1e-4)
  d <- data.frame(x = x, y = y)

  f <- suppressMessages(
    explore_dynamics(d, target = "y", predictors = "x"))$statistics$fit_x
  expect_identical(f$degree, 2L)

  ## The contract of the return value, applied: the fitted value at any point
  ## is the coefficients against ascending powers. If the order were reversed
  ## or shifted, this would not agree with the data anywhere.
  at <- c(-0.5, 0, 1, 2.5)
  rebuilt <- vapply(at, function(z) sum(f$coefficients * z^(0:f$degree)),
                    numeric(1))
  expect_equal(rebuilt, 3 + 0.5 * at - 1.25 * at^2, tolerance = 1e-3)
})

test_that("a truly linear relation rebuilds, whichever degree the AIC picks", {
  ## The word is NOT asserted here, and the reason is measured rather than
  ## assumed. On data that is linear by construction the AIC comparison picks
  ## "linear" about seventy-eight per cent of the time -- 78.8, 79.8 and 77.4
  ## at n = 120 for noise of 0.05, 0.5 and 2, and 77.6, 77.6 and 73.8 at
  ## n = 300 -- because two useless powers improve the likelihood by a
  ## chi-square on two degrees of freedom, which beats the penalty of four with
  ## a probability that does not shrink as n grows. The remaining fifth is
  ## called quadratic or cubic with tiny higher coefficients.
  ##
  ## So the invariant a caller can rely on is not the word: it is that the
  ## published fit reproduces the relation and keeps one slope sign over the
  ## range observed. That is what is asserted, and it is exactly what the word
  ## alone could not have told anyone.
  set.seed(20260819)
  x <- stats::runif(300, -3, 3)
  d <- data.frame(x = x, y = -4 + 2 * x + stats::rnorm(300, sd = 0.5))

  f <- suppressMessages(
    explore_dynamics(d, target = "y", predictors = "x"))$statistics$fit_x
  expect_identical(f$degree, unname(which.min(f$aic)))
  expect_length(f$coefficients, f$degree + 1L)

  at <- c(-2, 0, 2)
  rebuilt <- vapply(at, function(z) sum(f$coefficients * z^(0:f$degree)),
                    numeric(1))
  expect_equal(rebuilt, -4 + 2 * at, tolerance = 0.3)

  dcoef <- f$coefficients[-1L] * seq_len(f$degree)
  slopes <- vapply(seq(f$predictor_range[1L], f$predictor_range[2L],
                       length.out = 50L),
                   function(z) sum(dcoef * z^(seq_len(f$degree) - 1L)),
                   numeric(1))
  expect_true(all(slopes > 0))
})

test_that("every predictor gets its own fit, and none is confused with another", {
  set.seed(20260819)
  n <- 300
  a <- stats::runif(n, -2, 2)
  b <- stats::runif(n, -2, 2)
  d <- data.frame(a = a, b = b, y = 2 * a - 3 * b + stats::rnorm(n, sd = 1e-4))

  st <- suppressMessages(
    explore_dynamics(d, target = "y", predictors = c("a", "b")))$statistics
  expect_equal(st$fit_a$coefficients[2L], 2, tolerance = 0.2)
  expect_equal(st$fit_b$coefficients[2L], -3, tolerance = 0.2)
  expect_equal(st$fit_a$predictor_range, range(a))
  expect_equal(st$fit_b$predictor_range, range(b))
})

test_that("the same shape word covers a curve that turns and one that does not", {
  ## This is why the coefficients had to be published. Both relations below are
  ## called "quadratic". The first turns inside the range observed and has no
  ## direction over it; the second is the same parabola seen on one side of its
  ## vertex and rises throughout. The word cannot tell them apart; the sign of
  ## the derivative at the two ends of the published range can.
  set.seed(20260819)
  turning <- data.frame(x = stats::runif(400, -3, 3))
  turning$y <- turning$x^2 + stats::rnorm(400, sd = 1e-3)
  rising <- data.frame(x = stats::runif(400, 1, 4))
  rising$y <- rising$x^2 + stats::rnorm(400, sd = 1e-3)

  fit_of <- function(d) suppressMessages(
    explore_dynamics(d, target = "y", predictors = "x"))$statistics$fit_x
  slope_ends <- function(f) {
    dcoef <- f$coefficients[-1L] * seq_len(f$degree)
    vapply(f$predictor_range,
           function(z) sum(dcoef * z^(seq_len(f$degree) - 1L)), numeric(1))
  }

  ft <- fit_of(turning); fr <- fit_of(rising)
  expect_identical(ft$form, "quadratic")
  expect_identical(fr$form, "quadratic")
  expect_identical(ft$form, fr$form)

  expect_lt(prod(slope_ends(ft)), 0)  # the derivative changes sign: it turns
  expect_gt(prod(slope_ends(fr)), 0)  # it keeps one sign: it rises throughout
})
