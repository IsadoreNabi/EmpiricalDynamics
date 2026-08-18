# =============================================================================
# Regression tests over real Pareto fronts
#
# The fixture holds 100 equations from two real runs of this package -- a
# simulated SDE and the Lorenz system -- together with the (thinned) data they
# were discovered on. Every one of these equations used to fail
# cross-validation; the criterion below was fixed before the sweep was run and
# is simply that none of them does: every fold of every equation returns a
# finite error.
#
# See data-raw/make_fixtures.R for where the fixture comes from.
# =============================================================================

library(testthat)
library(EmpiricalDynamics)

fronts_fixture <- function() {
  path <- testthat::test_path("fixtures", "ed_fronts.rds")
  skip_if_not(file.exists(path), "the front fixture is not available")
  readRDS(path)
}

as_discovered <- function(expression, complexity) {
  structure(list(string = expression, expression = expression,
                 complexity = complexity),
            class = "symbolic_equation")
}

# A deterministic spread across the fronts, for the run that happens everywhere.
sampled_rows <- function(equations, how_many = 3L) {
  n <- nrow(equations)
  unique(round(seq(1, n, length.out = min(how_many, n))))
}

score_case <- function(case, rows, engine = "r") {
  vapply(rows, function(i) {
    cv <- cross_validate(
      as_discovered(case$equations$expression[i], case$equations$complexity[i]),
      data = case$data, response = case$response,
      k = 3, method = "block", refit_engine = engine, verbose = FALSE)
    if (all(is.finite(cv$rmse))) cv$mean_rmse else NA_real_
  }, numeric(1))
}

test_that("the fixture is the shape the tests expect", {
  cases <- fronts_fixture()

  expect_gte(length(cases), 4L)
  expect_equal(sum(vapply(cases, function(x) nrow(x$equations), integer(1))), 100L)
  for (case in cases) {
    expect_true(all(c("name", "response", "data", "equations") %in% names(case)))
    expect_true(case$response %in% names(case$data))
    expect_true(all(c("expression", "complexity") %in% names(case$equations)))
  }
})

test_that("a spread of real discovered equations is scored on every fold", {
  cases <- fronts_fixture()

  for (case in cases) {
    rows <- sampled_rows(case$equations)
    scores <- suppressWarnings(score_case(case, rows))
    failed <- case$equations$expression[rows][is.na(scores)]
    expect_true(length(failed) == 0L,
                info = paste(case$name, paste(failed, collapse = " | ")))
  }
})

test_that("every equation of every real front is scored, with no exception", {
  skip_on_cran()
  cases <- fronts_fixture()

  failures <- character(0)
  for (case in cases) {
    rows <- seq_len(nrow(case$equations))
    scores <- suppressWarnings(score_case(case, rows))
    failures <- c(failures,
                  paste0(case$name, ": ", case$equations$expression[rows][is.na(scores)]))
  }
  failures <- failures[!endsWith(failures, ": ")]

  # 22 of 22 used to fail on a front of this kind; the criterion is zero.
  expect_equal(length(failures), 0L,
               info = paste(failures, collapse = "\n"))
})

test_that("the two engines agree on the real fronts", {
  skip_on_cran()
  skip_if_not_installed("JuliaCall")
  skip_if_not(EmpiricalDynamics:::ed_julia_refit_ready(),
              "the Julia refitting engine is not available")

  cases <- fronts_fixture()

  for (case in cases) {
    rows <- sampled_rows(case$equations, how_many = 4L)
    in_r <- suppressWarnings(score_case(case, rows, engine = "r"))
    in_julia <- suppressWarnings(score_case(case, rows, engine = "julia"))

    # Neither engine may lose an equation the other keeps.
    expect_equal(is.na(in_r), is.na(in_julia), info = case$name)

    both <- !is.na(in_r) & !is.na(in_julia)
    if (any(both)) {
      relative <- abs(in_julia[both] - in_r[both]) / pmax(abs(in_r[both]), 1e-8)
      expect_lt(max(relative), 0.05)
    }
  }
})
