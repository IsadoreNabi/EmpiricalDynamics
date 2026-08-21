# =============================================================================
# Regression tests: sweeping a coefficient of a fit whose predict() does not
# read $coefficients
#
# `analyze_bifurcations()` accepted an `rstanarm` fit -- class
# `stanreg, glm, lm`, so `inherits(equation, "lm")` is TRUE -- and returned a
# well formed table with the SAME fixed point at every parameter value, none
# of them a true one. `analyze_fixed_points()` varies a coefficient by writing
# it into `$coefficients`; for a `stanreg` those are a posterior median, a
# summary of the draws, and `predict()` reads the draws. The write was
# discarded in silence.
#
# The gate that existed could not catch this. It checked the shape and the
# types of the returned table, and the returned table was correctly shaped.
# What was missing is an assertion that the fixed points VARY along the sweep,
# which is the whole content of the word "bifurcation", and it is the first
# test below.
#
# The referent throughout is algebra: the roots of the fitted quadratic with
# the swept coefficient substituted by hand, computed with polyroot().
# =============================================================================

# dZ = a0 + a1*Z + a2*Z^2. Sweeping a2 moves the positive root; the roots are
# in closed form, so the sweep can be checked against something that is not
# the package's own output.
make_sweep_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  Z <- seq(0.2, 6, length.out = n)
  data.frame(Z = Z,
             dZ = 0.05 + 0.67 * Z - 0.1 * Z^2 + stats::rnorm(n, sd = 0.05))
}

# The positive roots of the fitted quadratic once "I(Z^2)" is set to a2, kept
# to those inside the search window.
algebraic_roots <- function(fit, a2, window) {
  cf <- stats::coef(fit)
  r <- Re(polyroot(c(cf[["(Intercept)"]], cf[["I(Z)"]], a2)))
  sort(r[r >= window[1] & r <= window[2]])
}

# a2 in [-0.30, -0.15] keeps the positive root between about 2.3 and 4.6,
# and the search window [0.5, 6] keeps the negative root (near -0.05) out of
# it with room to spare, so the count of fixed points is 1 at every
# parameter value and the posterior cannot change it by wobbling.
SWEEP_RANGE <- c(-0.30, -0.15)
SWEEP_N <- 6L
SWEEP_WINDOW <- c(0.5, 6)

# -----------------------------------------------------------------------------
# 1. The assertion that was missing: a bifurcation sweep must actually vary
# -----------------------------------------------------------------------------

test_that("an lm sweep produces fixed points that MOVE, not n copies of one", {
  d <- make_sweep_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "I(Z^2)",
                              param_range = SWEEP_RANGE, n_param = SWEEP_N,
                              z_range = SWEEP_WINDOW)
  expect_identical(bif$mode, "coefficient")
  expect_gt(nrow(bif$data), 1L)
  # The defect being guarded against is a well formed table holding one
  # number repeated. One distinct value per parameter value that produced a
  # row is the minimum a sweep can mean.
  expect_equal(length(unique(round(bif$data$fixed_point, 8))),
               length(unique(bif$data$parameter_value)))
  expect_gt(stats::sd(bif$data$fixed_point), 0.1)
})

# -----------------------------------------------------------------------------
# 2. Positive control: the lm sweep against the algebra
# -----------------------------------------------------------------------------

test_that("the lm sweep reproduces the roots polyroot() gives", {
  d <- make_sweep_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "I(Z^2)",
                              param_range = SWEEP_RANGE, n_param = SWEEP_N,
                              z_range = SWEEP_WINDOW)
  for (a2 in unique(bif$data$parameter_value)) {
    got <- sort(bif$data$fixed_point[bif$data$parameter_value == a2])
    expect_equal(got, algebraic_roots(fit, a2, SWEEP_WINDOW),
                 tolerance = 1e-4,
                 info = paste("at I(Z^2) =", a2))
  }
  expect_equal(nrow(bif$status), SWEEP_N)
  expect_true(all(bif$status$status == "ok"))
  expect_equal(length(unique(bif$data$parameter_value)), SWEEP_N)
})

# -----------------------------------------------------------------------------
# 3. The refusal: an lm-inheriting object whose predict() ignores the
#    substitution must be rejected by ITS OWN cause, not by some other
#    complaint about the input
# -----------------------------------------------------------------------------

# A fit that inherits from lm and whose predict() reads a private copy of the
# coefficients, so that writing into $coefficients changes nothing. This is
# the property `rstanarm` has, isolated from `rstanarm`: the test does not
# need Stan to be installed to prove the package catches the condition.
fit_deaf <- function(d) {
  m <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  m$private <- stats::coef(m)
  class(m) <- c("deaf_fit", class(m))
  m
}

predict.deaf_fit <- function(object, newdata = NULL, ...) {
  m <- object
  class(m) <- setdiff(class(m), "deaf_fit")
  m$coefficients <- object$private        # the write is discarded, silently
  stats::predict(m, newdata = newdata, ...)
}

test_that("a fit whose predict() ignores $coefficients is refused, by name", {
  registerS3method("predict", "deaf_fit", predict.deaf_fit,
                   envir = environment())
  d <- make_sweep_data()
  m <- fit_deaf(d)

  # The premise of the test, measured rather than assumed: the object really
  # does discard the write, and really does inherit from lm.
  expect_true(inherits(m, "lm"))
  nd <- data.frame(Z = 1.5)
  before <- as.numeric(stats::predict(m, newdata = nd))
  m2 <- m
  m2$coefficients["I(Z^2)"] <- -5
  expect_equal(as.numeric(stats::predict(m2, newdata = nd)), before)

  expect_error(
    analyze_fixed_points(m, variable = "Z", range = SWEEP_WINDOW,
                         coefficient_values = list(`I(Z^2)` = -0.2)),
    class = "ed_substitution_ignored"
  )
  # ... and the sweep refuses too, instead of recording n_param failures.
  expect_error(
    analyze_bifurcations(m, variable = "Z", parameter = "I(Z^2)",
                         param_range = SWEEP_RANGE, n_param = SWEEP_N,
                         z_range = SWEEP_WINDOW),
    class = "ed_substitution_ignored"
  )
  # The message names the cause, not a generic complaint about the argument.
  err <- tryCatch(
    analyze_fixed_points(m, variable = "Z", range = SWEEP_WINDOW,
                         coefficient_values = list(`I(Z^2)` = -0.2)),
    ed_substitution_ignored = function(e) e)
  expect_match(conditionMessage(err), "does not read \\$coefficients")

  # Without a substitution there is nothing to discard: the same object still
  # answers, so the refusal is about the substitution and not about the class.
  expect_s3_class(
    analyze_fixed_points(m, variable = "Z", range = SWEEP_WINDOW),
    "data.frame")
})

test_that("a coefficient whose term is identically zero is its own refusal", {
  # W is held at 0, so the coefficient of W multiplies a column of zeros: no
  # honest object can respond to changing it. That is a different cause from
  # a predict() that ignores the write, and it gets a different condition.
  set.seed(5)
  Z <- rep(seq(0.2, 6, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W,
                  dZ = 0.05 + 0.67 * Z - 0.1 * Z^2 + 0.5 * W +
                    stats::rnorm(180, sd = 0.05))
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + W, data = d)
  expect_error(
    analyze_fixed_points(fit, variable = "Z", range = SWEEP_WINDOW,
                         exogenous_values = list(W = 0),
                         coefficient_values = list(W = 3)),
    class = "ed_substitution_inert"
  )
  # Held anywhere else, the same override is honoured.
  expect_s3_class(
    analyze_fixed_points(fit, variable = "Z", range = SWEEP_WINDOW,
                         exogenous_values = list(W = 1),
                         coefficient_values = list(W = 3)),
    "data.frame")
})

test_that("an override equal to the fitted value is a no-op, not an accusation", {
  d <- make_sweep_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  same <- unname(stats::coef(fit)["I(Z^2)"])
  # A window the fitted root actually falls in, so that "silent" is not
  # satisfied by an absence of fixed points to report.
  win <- c(0.5, 8)
  expect_silent(fps <- analyze_fixed_points(
    fit, variable = "Z", range = win,
    coefficient_values = list(`I(Z^2)` = same)))
  expect_equal(nrow(fps), 1L)
  plain <- analyze_fixed_points(fit, variable = "Z", range = win)
  expect_equal(fps$fixed_point, plain$fixed_point)
})

# -----------------------------------------------------------------------------
# 4. The expansion: rstanarm served through its draws
# -----------------------------------------------------------------------------

stan_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- make_sweep_data()
    fit <- rstanarm::stan_glm(dZ ~ I(Z) + I(Z^2), data = d, chains = 2,
                              iter = 1000, seed = 11, refresh = 0,
                              cores = 1)
    cache <<- list(d = d, stan = fit,
                   lm = stats::lm(dZ ~ I(Z) + I(Z^2), data = d))
    cache
  }
})

test_that("the reconstruction used for the draws matches the fit's own predict()", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  f <- stan_fixture()
  # This is what ed_certify_design() asserts at run time, asserted here from
  # outside: the per draw route, averaged, IS predict(). If it were not, the
  # sweep would be running on a reconstruction that is not the fitted model.
  nd <- data.frame(Z = seq(0.2, 6, length.out = 25))
  B <- ed_draws_of(f$stan)
  expect_false(is.null(B))
  design <- ed_design_of(f$stan, "Z", list())
  X <- design(nd$Z)[, colnames(B), drop = FALSE]
  expect_equal(unname(rowMeans(X %*% t(B))),
               as.numeric(stats::predict(f$stan, newdata = nd,
                                         type = "response")),
               tolerance = 1e-10)
  # And a plain lm must NOT be mistaken for a draws-bearing object.
  expect_null(ed_draws_of(f$lm))
})

test_that("the stanreg sweep varies, and agrees with the lm sweep", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  f <- stan_fixture()
  bif <- analyze_bifurcations(f$stan, variable = "Z", parameter = "I(Z^2)",
                              param_range = SWEEP_RANGE, n_param = SWEEP_N,
                              z_range = SWEEP_WINDOW, n_draws = 60L)
  expect_true(bif$posterior)
  expect_identical(bif$mode, "coefficient")
  expect_true("draw" %in% names(bif$data))
  expect_equal(bif$n_draws, 60L)

  # (a) It varies. This is the assertion the old gate did not have, applied
  #     to the object that broke it.
  expect_gt(nrow(bif$summary), 1L)
  expect_gt(stats::sd(bif$summary$fixed_point), 0.1)
  expect_equal(length(unique(round(bif$summary$fixed_point, 8))),
               nrow(bif$summary))

  # (b) It agrees with the lm sweep within the posterior interval. The
  #     contrast can fail: the lm root is a fixed number computed from a
  #     separate fit by polyroot(), and nothing forces the posterior to
  #     cover it.
  ref <- analyze_bifurcations(f$lm, variable = "Z", parameter = "I(Z^2)",
                              param_range = SWEEP_RANGE, n_param = SWEEP_N,
                              z_range = SWEEP_WINDOW)
  for (a2 in bif$summary$parameter_value) {
    row <- bif$summary[bif$summary$parameter_value == a2, ]
    alg <- algebraic_roots(f$lm, a2, SWEEP_WINDOW)
    expect_equal(length(alg), 1L, info = paste("at I(Z^2) =", a2))
    expect_gte(alg, row$fixed_point_lower)
    expect_lte(alg, row$fixed_point_upper)
    got <- ref$data$fixed_point[ref$data$parameter_value == a2]
    expect_gte(got, row$fixed_point_lower)
    expect_lte(got, row$fixed_point_upper)
  }

  # (c) The interval has width: a posterior sweep that collapsed to a point
  #     would pass (b) and mean nothing.
  expect_true(all(bif$summary$fixed_point_upper >
                    bif$summary$fixed_point_lower))
})

test_that("both modes of the sweep work on a stanreg and on an lm", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  set.seed(7)
  Z <- rep(seq(0.2, 6, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W,
                  dZ = 0.05 + 0.67 * Z - 0.1 * Z^2 + 0.2 * W +
                    stats::rnorm(180, sd = 0.05))
  # Wider window: the exogenous sweep of W slides the upper root from about
  # 6.5 to 7.1, and the count stays at 1 throughout.
  window <- c(0.5, 8)
  fit_lm <- stats::lm(dZ ~ I(Z) + I(Z^2) + W, data = d)
  fit_st <- rstanarm::stan_glm(dZ ~ I(Z) + I(Z^2) + W, data = d, chains = 2,
                               iter = 1000, seed = 7, refresh = 0, cores = 1)

  for (nm in c("lm", "stan")) {
    fit <- if (nm == "lm") fit_lm else fit_st
    coefw <- analyze_bifurcations(fit, variable = "Z", parameter = "I(Z^2)",
                                  param_range = SWEEP_RANGE,
                                  n_param = SWEEP_N, z_range = window,
                                  exogenous_values = list(W = 1),
                                  n_draws = 40L)
    expect_identical(coefw$mode, "coefficient")
    exog <- analyze_bifurcations(fit, variable = "Z", parameter = "W",
                                 param_range = c(-1, 1), n_param = 5,
                                 z_range = window, n_draws = 40L)
    expect_identical(exog$mode, "exogenous")
    # Both modes must vary, in both objects.
    for (b in list(coefw, exog)) {
      tbl <- if (isTRUE(b$posterior)) b$summary else b$data
      expect_gt(nrow(tbl), 1L)
      expect_gt(stats::sd(tbl$fixed_point), 1e-3)
    }
    expect_identical(isTRUE(coefw$posterior), nm == "stan")
    expect_identical(isTRUE(exog$posterior), nm == "stan")
  }
})

test_that("the posterior sweep says so when it is printed", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  f <- stan_fixture()
  bif <- analyze_bifurcations(f$stan, variable = "Z", parameter = "I(Z^2)",
                              param_range = SWEEP_RANGE, n_param = 3L,
                              z_range = SWEEP_WINDOW, n_draws = 30L)
  out <- paste(utils::capture.output(print(bif)), collapse = "\n")
  expect_match(out, "posterior draws")
  expect_match(out, "interval, not a number")
  # The row count must not be reported as a count of fixed points.
  expect_false(grepl(paste0("fixed points found: ", nrow(bif$data)), out,
                     fixed = TRUE))
})

test_that("check_qualitative_behavior counts branches, not draws", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  f <- stan_fixture()
  chk <- check_qualitative_behavior(
    f$stan, data = f$d, variable = "Z",
    expected_features = list(n_fixed_points = 2,
                             stability_pattern = c("unstable", "stable")))
  expect_true(chk$posterior)
  # Two fixed points inside the observed range, whatever the draw count is:
  # the reduction must not report 2 times n_draws rows.
  expect_equal(nrow(chk$checks$fixed_points), 2L)
  expect_true(all(chk$passed))
})

# -----------------------------------------------------------------------------
# 5. The refusals that no fixture reaches by accident
#
# Every branch below is one the package takes only when something is wrong.
# A suite that never visits them proves the happy path and nothing else.
# -----------------------------------------------------------------------------

# An object that inherits from lm, exposes a draws matrix covering its own
# coefficients -- so it is recognised as posterior bearing -- and whose
# predict() returns something that is NOT the average of those draws. The
# certification must refuse it rather than sweep a reconstruction that is not
# the fitted model.
fit_liar <- function(d) {
  m <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  cf <- stats::coef(m)
  m$draws <- matrix(rep(cf, each = 50), nrow = 50,
                    dimnames = list(NULL, names(cf)))
  class(m) <- c("liar_fit", class(m))
  m
}
as.matrix.liar_fit <- function(x, ...) x$draws
predict.liar_fit <- function(object, newdata = NULL, ...) {
  m <- object
  class(m) <- setdiff(class(m), "liar_fit")
  # Off by a constant: the draws and predict() no longer compose.
  stats::predict(m, newdata = newdata, ...) + 17
}

# A fit whose predict() is not reproducible. Whether a substitution "moved"
# anything cannot be established against a moving referent, so the probe must
# say so instead of concluding that the substitution took.
fit_restless <- function(d) {
  m <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  class(m) <- c("restless_fit", class(m))
  m
}
predict.restless_fit <- function(object, newdata = NULL, ...) {
  m <- object
  class(m) <- setdiff(class(m), "restless_fit")
  p <- stats::predict(m, newdata = newdata, ...)
  p + stats::rnorm(length(p), sd = 1)
}

test_that("a draws matrix that does not compose back into predict() is refused", {
  registerS3method("as.matrix", "liar_fit", as.matrix.liar_fit,
                   envir = environment())
  registerS3method("predict", "liar_fit", predict.liar_fit,
                   envir = environment())
  d <- make_sweep_data()
  m <- fit_liar(d)
  # The premise: it really is recognised as posterior bearing.
  expect_false(is.null(ed_draws_of(m)))
  expect_error(
    analyze_fixed_points(m, variable = "Z", range = SWEEP_WINDOW),
    class = "ed_uncertifiable_substitution")
  err <- tryCatch(analyze_fixed_points(m, variable = "Z",
                                       range = SWEEP_WINDOW),
                  ed_uncertifiable_substitution = function(e) e)
  expect_match(conditionMessage(err), "disagrees with its own predict\\(\\)")
})

test_that("a predict() that is not reproducible cannot testify to a substitution", {
  registerS3method("predict", "restless_fit", predict.restless_fit,
                   envir = environment())
  d <- make_sweep_data()
  m <- fit_restless(d)
  expect_error(
    analyze_fixed_points(m, variable = "Z", range = SWEEP_WINDOW,
                         coefficient_values = list(`I(Z^2)` = -0.2)),
    class = "ed_prediction_nondeterministic")
})

test_that("the arguments that govern the posterior are validated", {
  d <- make_sweep_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2), data = d)
  expect_error(analyze_fixed_points(fit, "Z", range = SWEEP_WINDOW,
                                    n_draws = 0), "'n_draws'")
  expect_error(analyze_fixed_points(fit, "Z", range = SWEEP_WINDOW,
                                    n_draws = c(1, 2)), "'n_draws'")
  expect_error(analyze_bifurcations(fit, "Z", "I(Z^2)",
                                    param_range = SWEEP_RANGE, n_param = 3,
                                    z_range = SWEEP_WINDOW, interval = 1),
               "'interval'")
  # The exported reduction refuses a frame that did not come from a posterior
  # sweep: without a draw column every row would be counted into one draw.
  plain <- analyze_fixed_points(fit, "Z", range = c(0.5, 8))
  expect_gt(nrow(plain), 0L)
  expect_error(ed_consensus_fixed_points(plain, 10L), "no 'draw' column")
  expect_error(ed_consensus_fixed_points(plain, 10L, interval = 0), "'interval'")
})

test_that("thinning is even, deterministic, and never invents draws", {
  B <- matrix(seq_len(200), ncol = 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(ed_thin_draws(B, 500), B)     # fewer draws than asked for
  expect_identical(ed_thin_draws(B, Inf), B)
  t1 <- ed_thin_draws(B, 10L)
  expect_equal(nrow(t1), 10L)
  expect_identical(t1, ed_thin_draws(B, 10L))    # no seed, same answer twice
  expect_equal(t1[1, "a"], B[1, "a"])
  expect_equal(t1[10, "a"], B[100, "a"])         # both ends of the posterior
  # Even spacing: the gaps between selected indices differ by at most one.
  gaps <- diff(match(t1[, "a"], B[, "a"]))
  expect_lte(max(gaps) - min(gaps), 1L)
})
