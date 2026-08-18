## The final choice of estimate_sde_iterative(), and the properties it has to
## keep. Which criterion ranks was decided by a pre-registered comparison over
## sixteen real trajectories with a known truth, not by argument; what these
## tests hold is that the machinery around that choice behaves.

make_traj <- function(fitted_list, complexity = NULL) {
  lapply(seq_along(fitted_list), function(i) {
    list(iteration = i,
         drift = list(expression = if (is.null(complexity)) "a * x" else complexity[i]),
         fitted = fitted_list[[i]])
  })
}

test_that("a candidate that predicts non-finite values is excluded and said so", {
  n <- 60L
  x <- seq(1, 5, length.out = n)
  target <- 2 * x
  good <- 2 * x
  bad  <- c(rep(Inf, 5L), 2 * x[-(1:5)])
  h <- make_traj(list(bad, good))
  ch <- select_from_trajectory(target, h, rep(1, n), data = data.frame(x = x))
  expect_identical(ch$selected, 2L)
  expect_identical(ch$excluded_nonfinite, 1L)
})

test_that("when nothing can be scored the choice falls back and declares it", {
  n <- 40L
  target <- rnorm(n)
  h <- make_traj(list(rep(Inf, n), rep(NaN, n)))
  ch <- select_from_trajectory(target, h, rep(1, n))
  expect_identical(ch$rule, "all_non_finite")
  expect_true(all(!is.finite(ch$scores)))
})

test_that("an exact tie goes to the earliest iteration", {
  n <- 50L
  x <- seq(1, 4, length.out = n)
  target <- 3 * x
  same <- 3 * x
  h <- make_traj(list(same, same, same))
  ch <- select_from_trajectory(target, h, rep(1, n), criterion = "deviance")
  expect_identical(ch$selected, 1L)
})

test_that("the weighted deviance is the ruler, and it is weighted", {
  n <- 40L
  target <- c(rep(0, n/2), rep(10, n/2))
  flat <- rep(0, n)
  ## Under flat weights the second half dominates; under weights that switch
  ## it off, the flat prediction is nearly perfect. Same numbers, two rulers.
  w_flat <- rep(1, n)
  w_first <- c(rep(1, n/2), rep(1e-8, n/2))
  expect_gt(weighted_deviance(target, flat, w_flat),
            weighted_deviance(target, flat, w_first))
})

test_that("blocked_cv refits the constants rather than evaluating the literals", {
  skip_if_not_installed("minpack.lm")
  set.seed(4)
  n <- 120L
  x <- seq(1, 6, length.out = n)
  d <- data.frame(x = x)
  target <- 2.5 * x + rnorm(n, sd = 0.05)
  ## The literal is deliberately wrong; refitting on the training blocks should
  ## recover something close to 2.5 and score far better than the literal does.
  score <- blocked_cv_score("9.9 * x", target, d, rep(1, n), n_blocks = 5L)
  literal <- weighted_deviance(target, 9.9 * x, rep(1, n))
  expect_true(is.finite(score))
  expect_lt(score, literal / 10)
})
