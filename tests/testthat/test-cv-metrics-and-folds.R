# =============================================================================
# Regression tests: what cross_validate() measures, and on which rows
#
# Four things were wrong at once, and each produced a well formed number that
# meant something other than what it said.
#
#   1. "rolling" was documented as walk-forward and trained on setdiff(1:n,
#      test), so the model saw the whole future of the series -- 41 rows after
#      the test point in a measured n = 50, k = 5 run -- and only 5 of 50 rows
#      were ever held out.
#   2. R-squared was taken against the TEST fold's own mean: a benchmark that
#      can only be computed after the answers are in, a different quantity in
#      every fold, and exactly zero for a single-row fold, which returned -Inf.
#   3. A fold that failed to refit was dropped by na.rm = TRUE, so a mean over
#      two folds printed exactly like a mean over five.
#   4. A glm was refitted with stats::lm(): family and link discarded.
#
# The referents are arithmetic identities, not tolerances.
# =============================================================================

cv_line <- function(n = 60, slope = 1.5, sd = 0.2, seed = 4) {
  set.seed(seed)
  d <- data.frame(time = seq_len(n), x = seq(1, 10, length.out = n))
  d$y <- 2 + slope * d$x + stats::rnorm(n, sd = sd)
  d
}

# -----------------------------------------------------------------------------
# 1. Walk-forward really walks forward
# -----------------------------------------------------------------------------

test_that("rolling trains only on the past, and tiles the tail", {
  n <- 50
  d <- cv_line(n)
  cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 5,
                       method = "rolling", verbose = FALSE)
  expect_equal(length(cv$fold_indices), 5L)
  for (i in seq_len(5)) {
    te <- cv$fold_indices[[i]]
    tr <- cv$train_indices[[i]]
    # The training set is strictly in the past. This is the whole content of
    # the words "walk-forward", and it is what the previous version broke.
    expect_true(max(tr) < min(te))
    expect_true(all(tr >= 1) && all(tr <= n))
    expect_true(all(te >= 1) && all(te <= n))
  }
  # The k test windows tile the tail: contiguous, and ending on the last row.
  all_te <- sort(unlist(cv$fold_indices))
  expect_equal(all_te, (n - 5 + 1):n)
  expect_identical(cv$window, "expanding")
})

test_that("the expanding window grows and the rolling window does not", {
  n <- 60
  d <- cv_line(n)
  exp_cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y",
                           k = 5, method = "rolling", horizon = 2,
                           window = "expanding", verbose = FALSE)
  rol_cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y",
                           k = 5, method = "rolling", horizon = 2,
                           window = "rolling", verbose = FALSE)
  n_exp <- vapply(exp_cv$train_indices, length, integer(1))
  n_rol <- vapply(rol_cv$train_indices, length, integer(1))
  expect_true(all(diff(n_exp) > 0))
  expect_true(all(diff(n_rol) == 0))
  expect_equal(n_exp[1], n_rol[1])
  # "sliding" is the same method under the other literature's name.
  sli_cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y",
                           k = 5, method = "rolling", horizon = 2,
                           window = "sliding", verbose = FALSE)
  expect_identical(sli_cv$window, "rolling")
  expect_identical(sli_cv$train_indices, rol_cv$train_indices)
  expect_equal(sli_cv$rmse, rol_cv$rmse)
})

test_that("fold construction never runs off the end of the data", {
  # create_rolling_folds(50, 5, 20) used to return 50:89 on 50 rows.
  for (n in c(40L, 50L, 97L)) {
    for (k in c(2L, 5L)) {
      for (h in c(1L, 3L, 7L)) {
        if (n - k * h < 2L) next
        f <- create_rolling_folds(n, k, h)
        idx <- unlist(lapply(f, function(z) c(z$train, z$test)))
        expect_true(all(idx >= 1L),
                    info = paste(n, k, h))
        expect_true(all(idx <= n), info = paste(n, k, h))
        expect_true(all(vapply(f, function(z) length(z$test), integer(1)) == h),
                    info = paste(n, k, h))
      }
    }
  }
  # A request that cannot be honoured is refused with the arithmetic that
  # refuses it, instead of producing indices that do not exist.
  expect_error(create_rolling_folds(20, 5, 20), "leaving 0 to train")
})

test_that("block folds partition every row, and an inconsistent size is refused", {
  # The old construction stepped by floor(n/k) from row 1 and truncated to k
  # starts: at n = 50, k = 3 it tested 1:16, 17:32, 33:48 and never held out
  # rows 49 and 50.
  for (n in c(50L, 97L, 100L)) {
    for (k in c(3L, 5L, 7L)) {
      f <- create_block_folds(n, k)
      expect_equal(length(f), k, info = paste(n, k))
      expect_equal(sort(unlist(f)), seq_len(n), info = paste(n, k))
      # Contiguous blocks, in order.
      expect_true(all(vapply(f, function(b) identical(b, min(b):max(b)),
                             logical(1))), info = paste(n, k))
    }
  }
  # block_size = 30 over 50 rows makes 2 blocks, not 5: the loop upstream used
  # to index past the end of the list and die with "subscript out of bounds".
  expect_error(create_block_folds(50, 5, 30), "makes 2 contiguous blocks")
  expect_error(
    cross_validate(stats::lm(y ~ x, data = cv_line(50)), cv_line(50),
                   response = "y", k = 5, method = "block", block_size = 30,
                   verbose = FALSE),
    "makes 2 contiguous blocks")
  expect_silent(f <- create_block_folds(50, 5, 10))
  expect_equal(length(f), 5L)
})

# -----------------------------------------------------------------------------
# 2. The R-squared, anchored on two exact identities
# -----------------------------------------------------------------------------

test_that("a model that predicts the training mean scores exactly zero", {
  d <- cv_line(60)
  # An intercept-only fit predicts the training mean at every held-out row, so
  # ss_res equals ss_tot by construction and R-squared is 0 -- exactly, not
  # approximately. Under the old definition, which benchmarked against the
  # TEST fold's mean, this identity could not even be written.
  cv <- cross_validate(stats::lm(y ~ 1, data = d), d, response = "y", k = 5,
                       method = "block", verbose = FALSE)
  # The identity is exact in exact arithmetic. In doubles the two sums travel
  # by different routes -- lm() reaches the mean through its QR decomposition,
  # mean() by summation -- so they agree to a few ulps rather than bit for
  # bit. The bound is machine epsilon, not a tolerance chosen to make the test
  # pass: 1e-12 is about 4500 times eps and still 10^3 times tighter than any
  # difference that could come from the benchmark being wrong.
  expect_lt(max(abs(cv$r_squared)), 1e-12)
  expect_lt(abs(cv$mean_r_squared), 1e-12)
  expect_lt(abs(cv$r_squared_pooled), 1e-12)
})

test_that("a perfect model scores exactly one", {
  d <- cv_line(60, sd = 0)
  cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 5,
                       method = "block", verbose = FALSE)
  expect_equal(cv$r_squared, rep(1, 5), tolerance = 1e-12)
  expect_equal(cv$r_squared_pooled, 1, tolerance = 1e-12)
})

test_that("a single-row test fold gives a number, not -Inf", {
  d <- cv_line(50)
  cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 5,
                       method = "rolling", horizon = 1, verbose = FALSE)
  expect_true(all(is.finite(cv$r_squared)))
  expect_true(is.finite(cv$mean_r_squared))
  expect_true(is.finite(cv$r_squared_pooled))
})

test_that("the pooled R-squared is a ratio of sums, not a mean of ratios", {
  d <- cv_line(60)
  cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 4,
                       method = "block", verbose = FALSE)
  expect_equal(cv$r_squared_pooled,
               1 - sum(cv$ss_res) / sum(cv$ss_tot))
  # And the per-fold benchmark really is the TRAINING mean.
  for (i in seq_len(4)) {
    tr <- cv$train_indices[[i]]
    te <- cv$fold_indices[[i]]
    expect_equal(cv$ss_tot[i], sum((d$y[te] - mean(d$y[tr]))^2))
  }
})

# -----------------------------------------------------------------------------
# 3. A failed fold cannot leave a number behind that looks complete
# -----------------------------------------------------------------------------

test_that("a fold that does not refit makes the k-fold figures NA", {
  d <- cv_line(60)
  # A formula the refit cannot satisfy on some training sets: an nls whose
  # start is far enough away that a fold or two do not converge is fragile to
  # arrange, so the failure is injected directly through a symbolic equation
  # whose expression is not evaluable.
  eq <- structure(list(expression = "a * x + stop_here(x)"),
                  class = "symbolic_equation")
  cv <- suppressWarnings(cross_validate(eq, d, response = "y", k = 4,
                                        method = "block", verbose = FALSE))
  expect_equal(cv$n_folds_ok, 0L)
  expect_equal(cv$prob_converged, 0)
  expect_true(is.na(cv$mean_rmse))
  expect_true(is.na(cv$mean_r_squared))
  expect_equal(length(cv$failures), 4L)
  expect_named(cv$failures, as.character(1:4))
  # The fallback figure is defined even when nothing refitted: it is the error
  # of predicting the training mean.
  expect_true(is.finite(cv$rmse_with_fallback))
  expect_equal(cv$r_squared_pooled_with_fallback, 0)
})

test_that("with every fold converged, nothing is NA and failures is empty", {
  d <- cv_line(60)
  cv <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 5,
                       method = "block", verbose = FALSE)
  expect_equal(cv$n_folds_ok, 5L)
  expect_equal(cv$prob_converged, 1)
  expect_length(cv$failures, 0L)
  expect_false(is.na(cv$mean_rmse))
  # The conditional figure coincides with the unconditional one exactly when
  # there is nothing to condition on.
  expect_identical(cv$mean_rmse, cv$rmse_given_converged)
  expect_identical(cv$mean_rmse, cv$rmse_with_fallback)
})

test_that("print says why the figures are NA before showing anything else", {
  d <- cv_line(60)
  eq <- structure(list(expression = "a * x + stop_here(x)"),
                  class = "symbolic_equation")
  cv <- suppressWarnings(cross_validate(eq, d, response = "y", k = 3,
                                        method = "block", verbose = FALSE))
  out <- paste(utils::capture.output(print(cv)), collapse = "\n")
  expect_match(out, "did not refit")
  expect_match(out, "not defined")
  expect_match(out, "naive benchmark")
  expect_match(out, "TRAINING rows")
})

# -----------------------------------------------------------------------------
# 4. A glm is refitted as a glm
# -----------------------------------------------------------------------------

test_that("a poisson glm is cross-validated as a poisson, not as an lm", {
  set.seed(6)
  d <- data.frame(x = seq(0.1, 3, length.out = 120))
  d$y <- stats::rpois(120, lambda = exp(0.3 + 0.9 * d$x))
  fit <- stats::glm(y ~ x, family = stats::poisson(), data = d)
  cv <- cross_validate(fit, d, response = "y", k = 4, method = "block",
                       verbose = FALSE)
  # Recompute both candidate refits on the folds the package actually used.
  rmse_lm <- rmse_glm <- numeric(4)
  for (i in 1:4) {
    te <- cv$fold_indices[[i]]
    tr <- cv$train_indices[[i]]
    p_lm <- stats::predict(stats::lm(y ~ x, data = d[tr, ]), newdata = d[te, ])
    p_gl <- stats::predict(stats::glm(y ~ x, family = stats::poisson(),
                                      data = d[tr, ]),
                           newdata = d[te, ], type = "response")
    rmse_lm[i] <- sqrt(mean((d$y[te] - p_lm)^2))
    rmse_glm[i] <- sqrt(mean((d$y[te] - p_gl)^2))
  }
  expect_equal(cv$rmse, rmse_glm)
  # The contrast can fail, and did: before the repair cv$rmse was rmse_lm.
  expect_false(isTRUE(all.equal(rmse_glm, rmse_lm)))
})

test_that("a posterior fit is refused rather than refitted by a different estimator", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")
  set.seed(11)
  d <- data.frame(x = seq(0.2, 6, length.out = 120))
  d$y <- 0.05 + 0.67 * d$x + stats::rnorm(120, sd = 0.3)
  fit <- rstanarm::stan_glm(y ~ x, data = d, chains = 2, iter = 800,
                            seed = 11, refresh = 0, cores = 1)
  expect_error(
    cross_validate(fit, d, response = "y", k = 3, method = "block",
                   verbose = FALSE),
    class = "ed_refit_substitutes_estimator")
})

# -----------------------------------------------------------------------------
# 5. The weights are honoured on the lm and glm paths
# -----------------------------------------------------------------------------

test_that("weights change the lm refit instead of being dropped", {
  d <- cv_line(60)
  w <- c(rep(1, 55), rep(1000, 5))
  a <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 3,
                      method = "block", verbose = FALSE)
  b <- cross_validate(stats::lm(y ~ x, data = d), d, response = "y", k = 3,
                      method = "block", weights = w, verbose = FALSE)
  # Same folds -- "block" is deterministic -- so any difference is the weights.
  expect_identical(a$fold_indices, b$fold_indices)
  expect_false(isTRUE(all.equal(a$rmse, b$rmse)))
  # And they are the weights of each fold's own training rows.
  for (i in seq_len(3)) {
    tr <- a$train_indices[[i]]
    te <- a$fold_indices[[i]]
    m <- stats::lm(y ~ x, data = d[tr, ], weights = w[tr])
    p <- stats::predict(m, newdata = d[te, ])
    expect_equal(b$rmse[i], sqrt(mean((d$y[te] - p)^2)))
  }
})

test_that("weights are honoured on the glm path too", {
  set.seed(8)
  d <- data.frame(x = seq(0.1, 3, length.out = 90))
  d$y <- stats::rpois(90, lambda = exp(0.2 + 0.8 * d$x))
  w <- c(rep(1, 80), rep(50, 10))
  fit <- stats::glm(y ~ x, family = stats::poisson(), data = d)
  a <- cross_validate(fit, d, response = "y", k = 3, method = "block",
                      verbose = FALSE)
  b <- cross_validate(fit, d, response = "y", k = 3, method = "block",
                      weights = w, verbose = FALSE)
  expect_false(isTRUE(all.equal(a$rmse, b$rmse)))
})
