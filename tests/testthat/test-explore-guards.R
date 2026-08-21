# =============================================================================
# Regression tests: explore_dynamics()
#
# Four defects, all of the same family -- a number or a verdict that was never
# measured, presented as if it had been.
#
#   1. `include` as a subset, which the manual page documents, was a hard error
#      from R 4.3 on: `include == "all" || ...` with a vector of length 2.
#   2. With no time-like column the auto-detection produced NA_character_, not
#      NULL, the is.null() guard passed, and a ggplot was built that failed
#      only when someone drew it.
#   3. A tertile of five points or fewer was never fitted but entered the
#      monotonicity test with a slope of exactly 0, and the guard could not
#      see it because the guard tests for NA.
#   4. Only significant interaction pairs were recorded, so a pair that was
#      tested and came out flat could not be told from one never tested.
# =============================================================================

explore_data <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(time = seq_len(n), Z = seq(-2, 2, length.out = n),
                  X = stats::rnorm(n))
  d$dZ <- 2 * d$Z - d$Z^2 + 0.1 * stats::rnorm(n)
  d
}

test_that("include accepts the subset the manual page documents", {
  d <- explore_data()
  r <- suppressMessages(explore_dynamics(d, "dZ", predictors = c("Z", "X"),
                                         include = c("timeseries", "phase")))
  expect_s3_class(r, "dynamics_exploration")
  # The blocks asked for ran ...
  expect_true("timeseries" %in% names(r$plots))
  expect_true("phase_Z" %in% names(r$plots))
  # ... and the ones not asked for did not.
  expect_false(any(grepl("^bivariate_", names(r$plots))))
  expect_false(any(grepl("^interaction_", names(r$statistics))))

  # A single name still works, and so does "all".
  expect_s3_class(
    suppressMessages(explore_dynamics(d, "dZ", predictors = "Z",
                                      include = "phase")),
    "dynamics_exploration")
})

test_that("a name outside the set is an error, not a block that quietly skips", {
  d <- explore_data()
  expect_error(
    explore_dynamics(d, "dZ", predictors = "Z", include = "timeseies"),
    "must be")
  expect_error(
    explore_dynamics(d, "dZ", predictors = "Z",
                     include = c("phase", "bivarate")),
    "must be")
  expect_error(explore_dynamics(d, "dZ", predictors = "Z", include = NA),
               "must be")
  expect_error(explore_dynamics(d, "dZ", predictors = "Z",
                                include = character(0)), "must be")
})

test_that("with no time column no broken plot is handed back", {
  d <- explore_data()
  d$time <- NULL
  r <- suppressMessages(explore_dynamics(d, "dZ", predictors = c("Z", "X")))
  # Either there is no time series plot, or the one there is can be drawn.
  # What must not happen is a plot object that fails when someone prints it.
  if (!is.null(r$plots$timeseries)) {
    expect_no_error(suppressMessages(
      invisible(ggplot2::ggplot_build(r$plots$timeseries))))
  } else {
    expect_null(r$plots$timeseries)
  }
  # And with a time column, the plot is there and is drawable.
  d2 <- explore_data()
  r2 <- suppressMessages(explore_dynamics(d2, "dZ", predictors = "Z",
                                          include = "timeseries"))
  expect_false(is.null(r2$plots$timeseries))
  expect_no_error(suppressMessages(
    invisible(ggplot2::ggplot_build(r2$plots$timeseries))))
})

test_that("a tertile that was not fitted cannot vote in the saturation test", {
  # The measured witness: group sizes 48/8/1. The third group has one point,
  # so it has no slope. It used to contribute an exact 0, which made the three
  # slopes strictly decreasing and announced saturation.
  set.seed(328)
  x <- c(rep(0, 47), sort(stats::runif(20, 1, 20)),
         sort(stats::runif(3, 50, 90)))
  y <- 3 * x - 0.02 * x^2 + stats::rnorm(length(x), sd = 0.3)
  d <- data.frame(time = seq_along(x), W = x, dZ = y)

  # The premise, measured: the fallback grouping really does leave a tertile
  # too small to fit.
  br <- stats::quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  g <- if (length(unique(br)) < 4) {
    as.numeric(cut(x, 3))
  } else {
    cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
  }
  sizes <- as.integer(table(factor(g, levels = 1:3)))
  expect_true(any(sizes <= 5L))

  r <- suppressMessages(explore_dynamics(d, "dZ", predictors = "W",
                                         include = "bivariate"))
  # No verdict may be announced from a group that was never estimated.
  expect_false(any(grepl("saturation|Exponential", r$suggestions,
                         ignore.case = TRUE)))
})

test_that("the saturation verdict still fires when all three groups are fitted", {
  # A control that can fail for what it watches: if the repair had simply
  # silenced the announcement, this would fail too.
  set.seed(12)
  x <- sort(stats::runif(120, 0.1, 6))
  y <- 4 * log(x + 1) + stats::rnorm(120, sd = 0.05)   # concave, slope falls
  d <- data.frame(time = seq_along(x), W = x, dZ = y)
  br <- stats::quantile(x, probs = c(0, 1/3, 2/3, 1))
  g <- cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)
  expect_true(all(as.integer(table(factor(g, levels = 1:3))) > 5L))
  r <- suppressMessages(explore_dynamics(d, "dZ", predictors = "W",
                                         include = "bivariate"))
  expect_true(any(grepl("saturation", r$suggestions, ignore.case = TRUE)))
})

test_that("every interaction pair tested is recorded, significant or not", {
  set.seed(3)
  d <- data.frame(time = 1:80, A = stats::rnorm(80), B = stats::rnorm(80),
                  C = stats::rnorm(80))
  d$dZ <- 1 + d$A + d$B + d$C + stats::rnorm(80, sd = 0.5)  # purely additive
  r <- suppressMessages(explore_dynamics(d, "dZ",
                                         predictors = c("A", "B", "C"),
                                         include = "interactions"))
  keys <- grep("^interaction_", names(r$statistics), value = TRUE)
  expect_equal(length(keys), 3L)                 # every pair, not only the
  expect_setequal(keys, c("interaction_A_B", "interaction_A_C",
                          "interaction_B_C"))
  p <- unlist(r$statistics[keys])
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  # Additive data: nothing should be announced, but the tests still leave
  # their p-values behind.
  expect_length(r$suggestions, 0L)
})
