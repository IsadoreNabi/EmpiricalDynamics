# =============================================================================
# Tests for symbolic_search.R
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

# Generate test data with known dynamics: dZ/dt = a + b*Z
set.seed(456)
n <- 50
test_data <- data.frame(
  Z = seq(0, 5, length.out = n),
  X = stats::rnorm(n)
)
# True: dZ = 1 + 0.5 * Z
test_data$dZ <- 1 + 0.5 * test_data$Z + 0.1 * stats::rnorm(n)

# -----------------------------------------------------------------------------
# Tests for fit_specified_equation
# -----------------------------------------------------------------------------

test_that("fit_specified_equation works with linear model", {
  # minpack.lm is required for the default method "LM"
  testthat::skip_if_not_installed("minpack.lm")
  
  result <- fit_specified_equation(
    "a + b * Z",
    data = test_data,
    response = "dZ"
  )
  
  expect_s3_class(result, "symbolic_equation")
  expect_true("a" %in% names(coef(result)))
  expect_true("b" %in% names(coef(result)))
  
  # Check coefficients are reasonable
  expect_true(abs(coef(result)["a"] - 1) < 0.5)
  expect_true(abs(coef(result)["b"] - 0.5) < 0.3)
})

test_that("fit_specified_equation works with provided starting values", {
  testthat::skip_if_not_installed("minpack.lm")
  
  result <- fit_specified_equation(
    "a + b * Z",
    data = test_data,
    response = "dZ",
    start = list(a = 0.5, b = 0.3)
  )
  
  expect_s3_class(result, "symbolic_equation")
})

test_that("fit_specified_equation handles nonlinear expressions", {
  testthat::skip_if_not_installed("minpack.lm")
  
  # Generate data with quadratic term
  test_data2 <- test_data
  test_data2$dZ <- 0.5 + 0.3 * test_data2$Z - 0.05 * test_data2$Z^2 + 0.05 * stats::rnorm(n)
  
  result <- fit_specified_equation(
    "a + b * Z + c * Z^2",
    data = test_data2,
    response = "dZ"
  )
  
  expect_s3_class(result, "symbolic_equation")
  expect_true("c" %in% names(coef(result)))
})

test_that("fit_specified_equation validates inputs", {
  expect_error(
    fit_specified_equation("a + b * Z", data = test_data, response = "nonexistent")
  )
})

test_that("predict works for fitted equations", {
  testthat::skip_if_not_installed("minpack.lm")
  
  result <- fit_specified_equation(
    "a + b * Z",
    data = test_data,
    response = "dZ"
  )
  
  pred <- predict(result, newdata = test_data)
  
  expect_length(pred, nrow(test_data))
  expect_false(any(is.na(pred)))
})

# -----------------------------------------------------------------------------
# Tests for symbolic search (R genetic)
# -----------------------------------------------------------------------------

test_that("symbolic_search_r_genetic runs and returns results", {
  testthat::skip_on_cran()  # Skip on CRAN due to computation time
  
  # Note: The function signature is symbolic_search(target, predictors, ...)
  # Parameters like population_size are internal/default in the R backend or 
  # handled differently, so we stick to the public API.
  
  result <- symbolic_search(
    target = test_data$dZ,
    predictors = test_data["Z"], # Pass as dataframe
    backend = "r_genetic",
    constraints = list(max_complexity = 5),
    n_runs = 1
  )
  
  expect_s3_class(result, "symbolic_search_result")
  expect_true("pareto_front" %in% names(result))
  expect_true("all_equations" %in% names(result))
})

# -----------------------------------------------------------------------------
# Tests for Pareto selection
# -----------------------------------------------------------------------------

test_that("select_equation chooses appropriate model", {
  # Create mock search result
  mock_result <- list(
    pareto_front = data.frame(
      complexity = c(3, 5, 8, 12),
      mse = c(0.5, 0.3, 0.25, 0.24), # Changed from rmse to mse to match select_equation logic if needed
      rmse = c(sqrt(0.5), sqrt(0.3), sqrt(0.25), sqrt(0.24)),
      equation = c("a", "a+b*Z", "a+b*Z+c*Z^2", "a+b*Z+c*Z^2+d*Z^3"),
      stringsAsFactors = FALSE
    ),
    all_equations = list(
      list(string = "a", complexity = 3, mse = 0.5, rmse = sqrt(0.5)),
      list(string = "a+b*Z", complexity = 5, mse = 0.3, rmse = sqrt(0.3))
    ),
    data = test_data,
    predictors = "Z"
  )
  class(mock_result) <- "symbolic_search_result"
  
  # Test different selection criteria
  eq_knee <- select_equation(mock_result, criterion = "knee")
  # select_equation returns a symbolic_equation object (list), not a character
  expect_s3_class(eq_knee, "symbolic_equation")
  
  eq_min <- select_equation(mock_result, criterion = "min_complexity")
  # Check the string component
  expect_equal(eq_min$string, "a")
})

# -----------------------------------------------------------------------------
# Tests for initial value estimation
# -----------------------------------------------------------------------------

test_that("estimate_initial_values provides reasonable estimates", {
  result <- estimate_initial_values(
    "a + b * Z",
    data = test_data,
    response = "dZ"
  )
  
  expect_type(result, "list")
  expect_true("a" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true(all(sapply(result, is.numeric)))
})

# -----------------------------------------------------------------------------
# Tests for coef and summary methods
# -----------------------------------------------------------------------------

test_that("coef method extracts coefficients correctly", {
  testthat::skip_if_not_installed("minpack.lm")
  
  eq <- fit_specified_equation(
    "a + b * Z",
    data = test_data,
    response = "dZ"
  )
  
  cf <- coef(eq)
  expect_type(cf, "double")
  expect_named(cf)
})

test_that("summary method works for symbolic_equation", {
  testthat::skip_if_not_installed("minpack.lm")
  
  eq <- fit_specified_equation(
    "a + b * Z",
    data = test_data,
    response = "dZ"
  )
  
  # Should not throw error
  expect_output(summary(eq))
})