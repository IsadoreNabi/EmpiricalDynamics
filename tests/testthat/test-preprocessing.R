# =============================================================================
# Tests for preprocessing.R
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

# Generate test data
set.seed(123)
n <- 100
test_time <- seq(0, 10, length.out = n)
test_signal <- sin(test_time) + 0.1 * stats::rnorm(n)

# -----------------------------------------------------------------------------
# Tests for derivative computation
# -----------------------------------------------------------------------------

test_that("compute_derivative works with finite differences", {
  # Use 't' instead of 'time' to match function signature
  result <- compute_derivative(test_signal, t = test_time, method = "fd")
  
  expect_type(result, "list")
  expect_true("derivative" %in% names(result))
  expect_equal(length(result$derivative), n)
  # Interior points should be non-NA
  expect_false(any(is.na(result$derivative[5:(n-5)])))
})

test_that("compute_derivative works with Savitzky-Golay", {
  testthat::skip_if_not_installed("signal")
  
  result <- compute_derivative(test_signal, t = test_time, method = "savgol")
  
  expect_type(result, "list")
  expect_equal(length(result$derivative), n)
})

test_that("compute_derivative works with spline", {
  result <- compute_derivative(test_signal, t = test_time, method = "spline")
  
  expect_type(result, "list")
  expect_equal(length(result$derivative), n)
  expect_false(any(is.na(result$derivative)))
})

test_that("compute_derivative TVR requires CVXR", {
  testthat::skip_if_not_installed("CVXR")
  
  result <- compute_derivative(test_signal, t = test_time, method = "tvr", lambda = 1)
  
  expect_type(result, "list")
  expect_equal(length(result$derivative), n)
  expect_s3_class(result, "tvr_derivative")
})

test_that("compute_derivative validates inputs", {
  expect_error(compute_derivative(test_signal, method = "invalid"))
  expect_error(compute_derivative("not numeric"))
})

test_that("compute_derivative spectral warns for non-periodic", {
  # Non-periodic data
  trend_signal <- test_time^2
  expect_warning(
    compute_derivative(trend_signal, t = test_time, method = "spectral"),
    regexp = "periodic|Gibbs"
  )
})

# -----------------------------------------------------------------------------
# Tests for variable specification
# -----------------------------------------------------------------------------

test_that("specify_variables creates correct structure", {
  df <- data.frame(
    time = 1:10,
    Z = stats::rnorm(10),
    X = stats::rnorm(10),
    Y = stats::rnorm(10)
  )
  
  result <- specify_variables(
    data = df,
    endogenous = "Z",
    coupled = "X",
    exogenous = "Y",
    time_col = "time"
  )
  
  # specify_variables returns the dataframe with an attribute "var_spec"
  expect_s3_class(result, "data.frame")
  expect_s3_class(result, "specified_data")
  
  spec <- attr(result, "var_spec")
  expect_type(spec, "list")
  
  expect_equal(spec$endogenous, "Z")
  expect_equal(spec$coupled, "X")
  expect_equal(spec$exogenous, "Y")
  expect_equal(spec$time_col, "time")
})

test_that("specify_variables validates column existence", {
  df <- data.frame(Z = 1:10)
  
  expect_error(
    specify_variables(df, endogenous = "nonexistent"),
    regexp = "not found"
  )
})

# -----------------------------------------------------------------------------
# Tests for derivative diagnostics
# -----------------------------------------------------------------------------

test_that("suggest_differentiation_method returns valid suggestion", {
  suggestion <- suggest_differentiation_method(test_signal, test_time)
  
  expect_type(suggestion, "list")
  expect_true("recommended" %in% names(suggestion))
  expect_true(suggestion$recommended %in% c("tvr", "savgol", "spline", "fd", "spectral"))
})

test_that("compare_differentiation_methods runs without error", {
  testthat::skip_if_not_installed("ggplot2")
  
  # Should not throw error
  result <- compare_differentiation_methods(
    test_signal, test_time,
    methods = c("fd", "spline"),
    plot = FALSE # Avoid plotting during tests
  )
  
  expect_type(result, "list")
})

# -----------------------------------------------------------------------------
# Tests for transformations
# -----------------------------------------------------------------------------

test_that("create_transformations generates expected columns", {
  df <- data.frame(
    Z = c(1, 2, 3, 4, 5),
    X = c(2, 4, 6, 8, 10)
  )
  
  # When 'variables' is provided without 'transformations', it uses defaults (log, sqrt, sq)
  result <- create_transformations(
    df,
    variables = c("Z", "X")
  )
  
  expect_true("log_Z" %in% names(result))
  expect_true("sqrt_X" %in% names(result))
  expect_true("sq_Z" %in% names(result))
  
  expect_equal(result$log_Z, log(abs(df$Z) + 1e-10))
  expect_equal(result$sqrt_X, sqrt(abs(df$X)))
})

# -----------------------------------------------------------------------------
# Tests for sampling diagnostics
# -----------------------------------------------------------------------------

test_that("diagnose_sampling_frequency provides diagnostics", {
  result <- diagnose_sampling_frequency(test_signal, test_time)
  
  expect_type(result, "list")
  expect_true("dt_mean" %in% names(result))
  expect_true("recommendation" %in% names(result) || "message" %in% names(result))
})