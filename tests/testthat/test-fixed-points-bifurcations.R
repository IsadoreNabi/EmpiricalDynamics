# =============================================================================
# Regression tests: analyze_fixed_points() / analyze_bifurcations()
#
# The referent throughout is algebra, not the package's own output.
# dZ = 2*Z - Z^2 has fixed points 0 (unstable, f' = +2) and 2 (stable,
# f' = -2). The lm path used to evaluate the bare terms of the formula --
# coefficients named "I(Z)" are not symbols of the deparsed RHS -- and
# returned the fixed points of Z + Z^2 without a word. A bifurcation sweep
# over a parameter absent from the equation used to return n_param identical
# copies labelled with the absent parameter's name.
# =============================================================================

make_quadratic_data <- function(n = 200) {
  Z <- seq(-1, 3, length.out = n)
  data.frame(Z = Z, dZ = 2 * Z - Z^2)
}

expect_quadratic_algebra <- function(fps, tol = 1e-4) {
  expect_s3_class(fps, "data.frame")
  expect_equal(nrow(fps), 2L)
  fps <- fps[order(fps$fixed_point), ]
  expect_equal(fps$fixed_point, c(0, 2), tolerance = tol)
  expect_equal(fps$stability, c("unstable", "stable"))
  expect_equal(fps$eigenvalue, c(2, -2), tolerance = 1e-2)
}

# -----------------------------------------------------------------------------
# The lm path must use the fitted coefficients
# -----------------------------------------------------------------------------

test_that("analyze_fixed_points on an lm recovers the algebra", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  # Sanity: the coefficients really are 2 and -1
  expect_equal(unname(stats::coef(fit)), c(2, -1), tolerance = 1e-8)

  fps <- analyze_fixed_points(fit, variable = "Z", range = c(-1, 3),
                              n_grid = 400)
  expect_quadratic_algebra(fps)
})

test_that("the lm path handles coefficients whose names collide with the variable", {
  # Written without I(): the coefficient is literally named "Z". Under the
  # old environment injection the coefficient binding was overwritten by the
  # grid value; under predict() there is nothing to collide.
  d <- make_quadratic_data()
  d$Z2 <- d$Z^2
  fit <- stats::lm(dZ ~ Z + Z2 + 0, data = d)
  # Z2 is a variable of the model, so it must be declared -- but declaring
  # it fixes it, which changes the geometry; instead check the I() form and
  # the plain form agree where both are valid: a 1-regressor model in Z.
  d1 <- data.frame(Z = seq(0, 4, length.out = 100))
  d1$dZ <- 6 - 3 * d1$Z # fixed point at 2, eigenvalue -3
  fit1 <- stats::lm(dZ ~ Z, data = d1)
  fps <- analyze_fixed_points(fit1, variable = "Z", range = c(0, 4))
  expect_equal(nrow(fps), 1L)
  expect_equal(fps$fixed_point, 2, tolerance = 1e-6)
  expect_equal(fps$stability, "stable")
  expect_equal(fps$eigenvalue, -3, tolerance = 1e-2)
})

test_that("a gaussian glm agrees with the lm on the response scale", {
  d <- make_quadratic_data()
  fit_lm <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  fit_glm <- stats::glm(dZ ~ I(Z) + I(Z^2) + 0, data = d,
                        family = stats::gaussian())
  fps_lm <- analyze_fixed_points(fit_lm, variable = "Z", range = c(-1, 3))
  fps_glm <- analyze_fixed_points(fit_glm, variable = "Z", range = c(-1, 3))
  expect_equal(fps_glm$fixed_point, fps_lm$fixed_point, tolerance = 1e-8)
  expect_equal(fps_glm$eigenvalue, fps_lm$eigenvalue, tolerance = 1e-6)
})

test_that("the nls path still recovers the algebra", {
  set.seed(42)
  d <- make_quadratic_data()
  d$dZ <- d$dZ + 1e-3 * stats::rnorm(nrow(d))
  fit <- stats::nls(dZ ~ a * Z + b * Z^2, data = d,
                    start = list(a = 1, b = -0.5))
  fps <- analyze_fixed_points(fit, variable = "Z", range = c(-1, 3),
                              n_grid = 400)
  expect_quadratic_algebra(fps, tol = 1e-2)
})

# -----------------------------------------------------------------------------
# Exogenous variables: honored, declared, and loud when misused
# -----------------------------------------------------------------------------

test_that("exogenous_values shift the fixed points of an lm as the algebra says", {
  # dZ = 2*Z - Z^2 + 0.5*W; at W = 2: z* = 1 +/- sqrt(2), f'(z) = 2 - 2z
  Z <- rep(seq(-1, 3, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W, dZ = 2 * Z - Z^2 + 0.5 * W)
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + W + 0, data = d)

  fps <- analyze_fixed_points(fit, variable = "Z", range = c(-2, 4),
                              n_grid = 400,
                              exogenous_values = list(W = 2))
  fps <- fps[order(fps$fixed_point), ]
  expect_equal(fps$fixed_point, c(1 - sqrt(2), 1 + sqrt(2)),
               tolerance = 1e-4)
  expect_equal(fps$stability, c("unstable", "stable"))
  expect_equal(fps$eigenvalue, c(2 * sqrt(2), -2 * sqrt(2)),
               tolerance = 1e-2)
  # The declaration travels with the result
  expect_identical(attr(fps, "exogenous_values"), list(W = 2))
})

test_that("a model variable without a declared value is an error, not a silent NA sweep", {
  Z <- rep(seq(-1, 3, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W, dZ = 2 * Z - Z^2 + 0.5 * W)
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + W + 0, data = d)
  expect_error(
    analyze_fixed_points(fit, variable = "Z", range = c(-2, 4)),
    "exogenous_values"
  )
})

test_that("the four silent traps are loud now", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)

  # 1. Unnamed list bound nothing
  expect_error(
    analyze_fixed_points(fit, "Z", range = c(-1, 3),
                         exogenous_values = list(2)),
    "named"
  )
  # 2. Declaring the search variable bound nothing
  expect_error(
    analyze_fixed_points(fit, "Z", range = c(-1, 3),
                         exogenous_values = list(Z = 1)),
    "search variable"
  )
  # 3. A name the equation never consults bound nothing
  expect_warning(
    analyze_fixed_points(fit, "Z", range = c(-1, 3),
                         exogenous_values = list(mu = 1)),
    "never consults"
  )
  # 4. A name colliding with a fitted coefficient silently overwrote the
  #    fit (expression path)
  set.seed(42)
  d2 <- make_quadratic_data()
  d2$dZ <- d2$dZ + 1e-3 * stats::rnorm(nrow(d2))
  fit_nls <- stats::nls(dZ ~ a * Z + b * Z^2, data = d2,
                        start = list(a = 1, b = -0.5))
  expect_error(
    analyze_fixed_points(fit_nls, "Z", range = c(-1, 3),
                         exogenous_values = list(a = 5)),
    "coefficient_values"
  )
})

test_that("a free symbol of an nls formula must be declared", {
  set.seed(7)
  Z <- rep(seq(-1, 3, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W,
                  dZ = 2 * Z - Z^2 + 0.5 * W + 1e-3 * stats::rnorm(180))
  fit <- stats::nls(dZ ~ a * Z + b * Z^2 + c0 * W, data = d,
                    start = list(a = 1, b = -0.5, c0 = 0.1))
  expect_error(
    analyze_fixed_points(fit, "Z", range = c(-2, 4)),
    "'W'"
  )
  fps <- analyze_fixed_points(fit, "Z", range = c(-2, 4), n_grid = 400,
                              exogenous_values = list(W = 2))
  fps <- fps[order(fps$fixed_point), ]
  expect_equal(fps$fixed_point, c(1 - sqrt(2), 1 + sqrt(2)),
               tolerance = 1e-2)
})

test_that("coefficient_values overrides a fitted coefficient deliberately", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  # Override the linear coefficient: dZ = 3*Z - Z^2 has fixed points 0 and 3
  fps <- analyze_fixed_points(fit, "Z", range = c(-1, 4), n_grid = 400,
                              coefficient_values = list(`I(Z)` = 3))
  fps <- fps[order(fps$fixed_point), ]
  expect_equal(fps$fixed_point, c(0, 3), tolerance = 1e-4)
  # A name that is not a fitted coefficient is an error
  expect_error(
    analyze_fixed_points(fit, "Z", range = c(-1, 3),
                         coefficient_values = list(mu = 3)),
    "not fitted coefficients"
  )
})

test_that("a variable absent from the equation is an error, not a constant sweep", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  expect_error(
    analyze_fixed_points(fit, variable = "Q", range = c(-1, 3)),
    "does not appear"
  )
})

# -----------------------------------------------------------------------------
# Bifurcations: the parameter must exist, and the grid must be reconstructible
# -----------------------------------------------------------------------------

test_that("a bifurcation parameter absent from the equation is an error", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  expect_error(
    analyze_bifurcations(fit, variable = "Z", parameter = "mu",
                         param_range = c(-1, 1), n_param = 5,
                         z_range = c(-1, 3)),
    "appears nowhere"
  )
  expect_error(
    analyze_bifurcations(fit, variable = "Z", parameter = "Z",
                         param_range = c(-1, 1)),
    "cannot be the variable"
  )
})

test_that("varying a fitted nls coefficient tracks the transcritical bifurcation", {
  set.seed(1)
  d <- data.frame(Z = seq(-1, 3, length.out = 120))
  d$dZ <- 2 * d$Z - d$Z^2 + 1e-3 * stats::rnorm(120)
  fit <- stats::nls(dZ ~ r * Z - Z^2, data = d, start = list(r = 1))

  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "r",
                              param_range = c(-1, 1), n_param = 9,
                              z_range = c(-3, 3))
  expect_s3_class(bif, "bifurcation_analysis")
  expect_identical(bif$mode, "coefficient")
  # Every requested parameter value is accounted for
  expect_equal(nrow(bif$status), 9L)
  expect_true(all(bif$status$status %in% c("ok", "no_fixed_points")))
  # dZ = r*Z - Z^2 has fixed points {0, r}; at r = 1 the nonzero one is
  # stable with eigenvalue -r
  at_r1 <- bif$data[abs(bif$data$parameter_value - 1) < 1e-9, ]
  at_r1 <- at_r1[order(at_r1$fixed_point), ]
  expect_equal(at_r1$fixed_point, c(0, 1), tolerance = 1e-2)
  expect_equal(at_r1$stability, c("unstable", "stable"))
  # ... and at r = -1 the stability is exchanged (transcritical)
  at_rm1 <- bif$data[abs(bif$data$parameter_value + 1) < 1e-9, ]
  at_rm1 <- at_rm1[order(at_rm1$fixed_point), ]
  expect_equal(at_rm1$fixed_point, c(-1, 0), tolerance = 1e-2)
  expect_equal(at_rm1$stability, c("unstable", "stable"))
})

test_that("a parameter value with no fixed points stays visible in the status", {
  # dZ = r - Z^2: two fixed points for r > 0, none for r < 0 (saddle-node)
  set.seed(2)
  d <- data.frame(Z = seq(-2, 2, length.out = 120))
  d$dZ <- 1 - d$Z^2 + 1e-3 * stats::rnorm(120)
  fit <- stats::nls(dZ ~ r - Z^2, data = d, start = list(r = 0.5))

  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "r",
                              param_range = c(-1, 1), n_param = 11,
                              z_range = c(-3, 3))
  expect_equal(nrow(bif$status), 11L)
  neg <- bif$status$parameter_value < -1e-9
  pos <- bif$status$parameter_value > 1e-9
  expect_true(all(bif$status$status[neg] == "no_fixed_points"))
  expect_true(all(bif$status$status[pos] == "ok"))
  expect_true(all(bif$status$n_fixed_points[pos] == 2L))
  # The data only carries real findings, but the requested grid is intact
  expect_true(all(bif$data$parameter_value > -1e-9))
  # At r = 1 the fixed points are +/- 1
  at_r1 <- bif$data[abs(bif$data$parameter_value - 1) < 1e-9, ]
  expect_equal(sort(at_r1$fixed_point), c(-1, 1), tolerance = 1e-2)
})

test_that("varying an lm coefficient by its term label works", {
  d <- make_quadratic_data()
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = d)
  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "I(Z)",
                              param_range = c(1, 3), n_param = 5,
                              z_range = c(-1, 4))
  expect_identical(bif$mode, "coefficient")
  # dZ = c*Z - Z^2: the nonzero fixed point equals the coefficient
  nonzero <- bif$data[abs(bif$data$fixed_point) > 0.1, ]
  expect_equal(nonzero$fixed_point, nonzero$parameter_value,
               tolerance = 1e-4)
})

test_that("an exogenous variable can be the bifurcation parameter", {
  Z <- rep(seq(-1, 3, length.out = 60), times = 3)
  W <- rep(c(-1, 0, 1), each = 60)
  d <- data.frame(Z = Z, W = W, dZ = 2 * Z - Z^2 + 0.5 * W)
  fit <- stats::lm(dZ ~ I(Z) + I(Z^2) + W + 0, data = d)
  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "W",
                              param_range = c(0, 2), n_param = 3,
                              z_range = c(-2, 4))
  expect_identical(bif$mode, "exogenous")
  # At W = 2: z* = 1 +/- sqrt(2)
  at_w2 <- bif$data[abs(bif$data$parameter_value - 2) < 1e-9, ]
  expect_equal(sort(at_w2$fixed_point), c(1 - sqrt(2), 1 + sqrt(2)),
               tolerance = 1e-3)
})

test_that("print.bifurcation_analysis states the grid accounting", {
  set.seed(2)
  d <- data.frame(Z = seq(-2, 2, length.out = 120))
  d$dZ <- 1 - d$Z^2 + 1e-3 * stats::rnorm(120)
  fit <- stats::nls(dZ ~ r - Z^2, data = d, start = list(r = 0.5))
  bif <- analyze_bifurcations(fit, variable = "Z", parameter = "r",
                              param_range = c(-1, 1), n_param = 11,
                              z_range = c(-3, 3))
  out <- paste(utils::capture.output(print(bif)), collapse = "\n")
  expect_match(out, "11 parameter values")
  expect_match(out, "without")
})

# -----------------------------------------------------------------------------
# The symbolic_equation path (needs minpack.lm)
# -----------------------------------------------------------------------------

test_that("the symbolic_equation path recovers the algebra and rejects the traps", {
  testthat::skip_if_not_installed("minpack.lm")
  set.seed(3)
  d <- make_quadratic_data()
  d$dZ <- d$dZ + 1e-2 * stats::rnorm(nrow(d))
  eq <- fit_specified_equation("a * Z + b * Z^2", data = d, response = "dZ")

  fps <- analyze_fixed_points(eq, variable = "Z", range = c(-1, 3),
                              n_grid = 400)
  expect_quadratic_algebra(fps, tol = 5e-2)

  expect_error(
    analyze_fixed_points(eq, "Z", range = c(-1, 3),
                         exogenous_values = list(a = 5)),
    "coefficient_values"
  )
  expect_error(
    analyze_bifurcations(eq, variable = "Z", parameter = "mu",
                         param_range = c(-1, 1), n_param = 3),
    "appears nowhere"
  )
  bif <- analyze_bifurcations(eq, variable = "Z", parameter = "a",
                              param_range = c(1, 3), n_param = 5,
                              z_range = c(-1, 4))
  expect_identical(bif$mode, "coefficient")
  nonzero <- bif$data[abs(bif$data$fixed_point) > 0.1, ]
  # dZ = a*Z + b*Z^2 with b ~ -1: nonzero fixed point at -a/b ~ a
  expect_equal(nonzero$fixed_point, -nonzero$parameter_value /
                 unname(stats::coef(eq)["b"]),
               tolerance = 1e-2)
})
