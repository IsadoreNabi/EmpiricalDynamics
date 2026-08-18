# =============================================================================
# Tests for cross-validating an equation that came out of a search
#
# A discovered equation arrives with its constants written in as literals, which
# is exactly the shape `get_pareto_set()` hands back and the shape that used to
# fail on every single candidate. The equations below are built by hand in that
# shape so that the test needs neither Julia nor a search to run.
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

discovered <- function(expression, complexity = 5L) {
  structure(
    list(string = expression, expression = expression,
         complexity = complexity),
    class = "symbolic_equation")
}

cv_data <- function(n = 90) {
  set.seed(7)
  d <- data.frame(x = seq(-3, 3, length.out = n))
  d$Z3 <- seq(0.5, 2, length.out = n)
  d$y <- 3.7 * d$x^2 + 5 + stats::rnorm(n, sd = 0.4)
  d
}

test_that("a discovered equation is scored on every fold", {
  d <- cv_data()

  cv <- cross_validate(discovered("(x * x) * 3.9305288261198923"),
                       data = d, response = "y", k = 5, method = "block",
                       verbose = FALSE)

  expect_s3_class(cv, "cv_result")
  expect_length(cv$rmse, 5L)
  # The whole point: not one fold falls over.
  expect_true(all(is.finite(cv$rmse)))
  expect_true(is.finite(cv$mean_rmse))
})

test_that("the constants are what gets re-estimated, and it is reported", {
  d <- cv_data()

  cv <- cross_validate(discovered("(x * x) * 3.9305288261198923"),
                       data = d, response = "y", k = 3, method = "block",
                       verbose = FALSE)

  expect_true(cv$refitted)
  expect_type(cv$parameterization, "character")
  expect_true(grepl("c1", cv$parameterization, fixed = TRUE))
  expect_identical(cv$refit_engine, "r")
})

test_that("an equation with no constant is scored without being refitted", {
  d <- cv_data()

  cv <- cross_validate(discovered("x"), data = d, response = "y",
                       k = 3, method = "block", verbose = FALSE)

  expect_false(cv$refitted)
  expect_null(cv$parameterization)
  expect_true(all(is.finite(cv$rmse)))
})

test_that("a bare constant equation is scored too", {
  d <- cv_data()

  expect_silent(
    cv <- cross_validate(discovered("-0.049732085687063235"),
                         data = d, response = "y", k = 3, method = "block",
                         verbose = FALSE)
  )
  expect_true(all(is.finite(cv$rmse)))
  expect_true(cv$refitted)
})

test_that("the forms a real front contains all survive", {
  d <- cv_data()

  fronts <- c(
    "(x * x) * 3.9305288261198923",
    "-0.049732085687063235",
    "abs(2.5743296383495817e-14)",
    "x",
    "x * 8.52617055983837",
    "-0.035493998791350026 + (8.526144599521835 * x)",
    "sin(x) + (8.512849039087632 * x)",
    "((x * -0.27984256967296445) + 10.403285992498093) * x",
    "(4.0763988901487735 * x) - Z3",
    "square(x) * 2.1",
    "inv(Z3) * 1.5"
  )

  for (expression in fronts) {
    cv <- cross_validate(discovered(expression), data = d, response = "y",
                         k = 3, method = "block", verbose = FALSE)
    expect_true(all(is.finite(cv$rmse)), info = expression)
  }
})

test_that("the fields the object already promised are still there", {
  d <- cv_data()

  cv <- cross_validate(discovered("(x * x) * 3.93"), data = d, response = "y",
                       k = 4, method = "block", verbose = FALSE)

  for (field in c("rmse", "mae", "r_squared", "mean_rmse", "sd_rmse",
                  "mean_mae", "mean_r_squared", "predictions", "fold_indices",
                  "method", "k")) {
    expect_true(field %in% names(cv), info = field)
  }
  expect_equal(cv$k, 4L)
  expect_identical(cv$method, "block")
  expect_length(cv$predictions, 4L)
  expect_true(all(c("actual", "predicted", "residual") %in%
                    names(cv$predictions[[1]])))
})

test_that("an equation the user specified still takes its own coefficients", {
  skip_if_not_installed("minpack.lm")
  d <- cv_data()

  eq <- fit_specified_equation("a * x^2 + b", data = d, response = "y",
                               start = list(a = 1, b = 1))
  cv <- cross_validate(eq, data = d, response = "y", k = 3, method = "block",
                       verbose = FALSE)

  expect_true(all(is.finite(cv$rmse)))
  expect_true(cv$refitted)
  expect_null(cv$parameterization)
})

test_that("weights are checked and used", {
  d <- cv_data()
  eq <- discovered("(x * x) * 3.93")

  expect_error(cross_validate(eq, data = d, response = "y", k = 3,
                              weights = c(1, 2), verbose = FALSE),
               "one value per row")
  expect_error(cross_validate(eq, data = d, response = "y", k = 3,
                              weights = rep(-1, nrow(d)), verbose = FALSE),
               "non-negative")

  w <- rep(1, nrow(d))
  w[d$x < 0] <- 50
  weighted <- cross_validate(eq, data = d, response = "y", k = 3,
                             method = "block", weights = w, verbose = FALSE)
  plain <- cross_validate(eq, data = d, response = "y", k = 3,
                          method = "block", verbose = FALSE)

  expect_true(all(is.finite(weighted$rmse)))
  expect_false(isTRUE(all.equal(weighted$rmse, plain$rmse)))
})

test_that("cross-validation still works for lm and nls equations", {
  d <- cv_data()

  model <- stats::lm(y ~ x, data = d)
  cv <- cross_validate(model, data = d, response = "y", k = 3,
                       method = "block", verbose = FALSE)

  expect_s3_class(cv, "cv_result")
  expect_true(all(is.finite(cv$rmse)))
})

test_that("the julia engine is refused politely when it is not there", {
  d <- cv_data()

  # Whatever the machine has, asking for the Julia engine on something that is
  # not a symbolic equation is a mistake this package names itself.
  expect_error(
    cross_validate(stats::lm(y ~ x, data = d), data = d, response = "y",
                   k = 3, refit_engine = "julia", verbose = FALSE),
    "symbolic equations only")
})

test_that("the julia engine agrees with the R engine when it is available", {
  skip_on_cran()
  skip_if_not_installed("JuliaCall")
  skip_if_not(EmpiricalDynamics:::ed_julia_refit_ready(),
              "the Julia refitting engine is not available")

  d <- cv_data()
  eq <- discovered("(x * x) * 3.9305288261198923")

  in_r <- cross_validate(eq, data = d, response = "y", k = 3, method = "block",
                         verbose = FALSE)
  in_julia <- cross_validate(eq, data = d, response = "y", k = 3,
                             method = "block", refit_engine = "julia",
                             verbose = FALSE)

  expect_true(all(is.finite(in_julia$rmse)))
  expect_identical(in_julia$refit_engine, "julia")
  # Two independent optimisers on the same objective: they must land together.
  expect_equal(in_julia$rmse, in_r$rmse, tolerance = 1e-4)
})

test_that("the julia engine reports a constant it could not move", {
  skip_on_cran()
  skip_if_not_installed("JuliaCall")
  skip_if_not(EmpiricalDynamics:::ed_julia_refit_ready(),
              "the Julia refitting engine is not available")

  d <- cv_data()

  # The search's restarts perturb the constants multiplicatively, so a constant
  # left at zero cannot move, and `abs(c)` at zero is a kink the optimiser does
  # not leave. What is under test is that this is said out loud rather than
  # passed off as a refit.
  stuck <- discovered("abs(-4.163336342344337e-16)")

  expect_warning(
    cv <- cross_validate(stuck, data = d, response = "y", k = 3,
                         method = "block", refit_engine = "julia",
                         verbose = FALSE),
    "did not move off its starting constants")

  expect_true(length(cv$stalled_folds) > 0L)
  expect_true(all(is.finite(cv$rmse)))

  # The R engine has no such trouble and says nothing.
  expect_silent(
    in_r <- cross_validate(stuck, data = d, response = "y", k = 3,
                           method = "block", verbose = FALSE))
  expect_length(in_r$stalled_folds, 0L)
})

test_that("bootstrap and sensitivity work on a discovered equation", {
  d <- cv_data()
  eq <- discovered("(x * x) * 3.9305288261198923")

  boot <- bootstrap_parameters(eq, data = d, response = "y", n_boot = 15)
  expect_true(is.data.frame(boot))
  expect_equal(nrow(boot), 1L)
  expect_true(all(is.finite(boot$estimate)))
  expect_true(boot$ci_lower <= boot$estimate + 1e-6)
  expect_true(boot$ci_upper >= boot$estimate - 1e-6)

  sens <- sensitivity_analysis(eq, data = d, response = "y")
  expect_true(is.data.frame(sens))
  expect_equal(nrow(sens), 1L)
  expect_true(all(is.finite(sens$sensitivity)))
})
