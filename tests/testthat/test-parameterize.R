# =============================================================================
# Tests for parameterize_equation()
#
# The forms exercised here are the ones a real Pareto front actually contains:
# a bare literal, a signed literal, scientific notation, a literal nested inside
# a unary operator, an equation with no literal at all, and variable names that
# carry digits.
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

test_that("a fitted constant becomes a free parameter starting at its value", {
  got <- parameterize_equation("(x * x) * 3.9305288261198923")

  expect_equal(got$n_parameters, 1L)
  expect_named(got$start, "c1")
  expect_equal(got$start$c1, 3.9305288261198923)
  expect_false(grepl("3.93", got$expression, fixed = TRUE))
  expect_true(grepl("c1", got$expression, fixed = TRUE))
})

test_that("every literal gets its own parameter", {
  got <- parameterize_equation("-0.0354 + (8.5261 * YmX)")

  expect_equal(got$n_parameters, 2L)
  expect_named(got$start, c("c1", "c2"))
  expect_equal(unname(unlist(got$start)), c(-0.0354, 8.5261))
})

test_that("a leading minus is folded into the starting value, not left as an operator", {
  got <- parameterize_equation("-0.049732085687063235")

  expect_equal(got$n_parameters, 1L)
  expect_equal(got$start$c1, -0.049732085687063235)
  # The sign belongs to the coefficient: the expression is just the parameter.
  expect_equal(trimws(got$expression), "c1")
})

test_that("scientific notation and literals inside unary calls are picked up", {
  got <- parameterize_equation("abs(2.5743296383495817e-14)")

  expect_equal(got$n_parameters, 1L)
  expect_equal(got$start$c1, 2.5743296383495817e-14)
  expect_true(grepl("abs(c1)", got$expression, fixed = TRUE))
})

test_that("an equation with no literal is returned untouched", {
  got <- parameterize_equation("YmX")

  expect_equal(got$n_parameters, 0L)
  expect_length(got$start, 0L)
  expect_equal(trimws(got$expression), "YmX")
})

test_that("digits inside a variable name are not constants", {
  got <- parameterize_equation("(4.0763988901487735 * sinX) - Z3")

  expect_equal(got$n_parameters, 1L)
  expect_true(grepl("Z3", got$expression, fixed = TRUE))
  expect_equal(got$start$c1, 4.0763988901487735)
})

test_that("an integer exponent is structure and is left alone by default", {
  protected <- parameterize_equation("x^2 * 3.5")
  expect_equal(protected$n_parameters, 1L)
  expect_true(grepl("x^2", protected$expression, fixed = TRUE))

  released <- parameterize_equation("x^2 * 3.5", protect_exponents = FALSE)
  expect_equal(released$n_parameters, 2L)
  expect_false(grepl("x^2", released$expression, fixed = TRUE))
})

test_that("generated names never collide with the columns of the data", {
  data <- data.frame(c1 = 1:5, x = 1:5)
  got <- parameterize_equation("x * 2.5 + c1", data = data)

  expect_equal(got$n_parameters, 1L)
  expect_false("c1" %in% names(got$start))
  expect_true(grepl("c1", got$expression, fixed = TRUE))  # the column survives
})

test_that("generated names never collide with names already in the expression", {
  got <- parameterize_equation("c1 * 2.5")

  expect_equal(got$n_parameters, 1L)
  expect_false("c1" %in% names(got$start))
})

test_that("the parameterised equation computes exactly what the original did", {
  set.seed(11)
  data <- data.frame(
    x    = stats::runif(20, 1, 3),
    Z3   = stats::runif(20),
    sinX = stats::runif(20),
    absX = stats::runif(20)
  )

  originals <- c(
    "(x * x) * 3.93",
    "-0.0497",
    "abs(2.57e-14)",
    "((4.07 * sinX) - Z3) * (2.76 - absX)",
    "x^2 * 3.5",
    "square(x) * 2.1 + inv(x)",
    "7.3075 * (sinX + (-0.2468 * Z3))"
  )

  for (original in originals) {
    got <- parameterize_equation(original, data = data)

    before <- rep_len(as.numeric(eval(
      parse(text = original),
      envir = EmpiricalDynamics:::create_eval_env(data))), nrow(data))
    after <- rep_len(as.numeric(eval(
      parse(text = got$expression),
      envir = EmpiricalDynamics:::create_eval_env(data, got$start))), nrow(data))

    expect_equal(after, before, info = original)
  }
})

test_that("a symbolic_equation can be handed over directly", {
  equation <- structure(
    list(string = "(x * x) * 3.93", expression = "(x * x) * 3.93"),
    class = "symbolic_equation")

  expect_equal(parameterize_equation(equation)$n_parameters, 1L)
})

test_that("the prefix is honoured", {
  got <- parameterize_equation("x * 2.5", prefix = "theta")
  expect_named(got$start, "theta1")
})

test_that("bad input is refused rather than guessed at", {
  expect_error(parameterize_equation(c("a", "b")), "single character string")
  expect_error(parameterize_equation("x * "), "Could not parse")
  expect_error(parameterize_equation("x * 2", prefix = ""), "non-empty")
})
