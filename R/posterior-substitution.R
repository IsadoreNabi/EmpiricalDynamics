# =============================================================================
# Substituting a coefficient into a fitted object, and proving it took
# =============================================================================
#
# `analyze_fixed_points()` varies a fitted coefficient by writing it into the
# object -- `equation$coefficients[nm] <- v` -- and then evaluating the object
# through its own `predict()`. Evaluating through `predict()` is what made the
# lm and glm paths correct: coefficient names are term labels ("I(Z^2)",
# "(Intercept)", factor contrasts), not symbols of the formula text, so an
# expression rebuilt by hand from the formula does not carry them.
#
# But the assignment and the evaluation are two different assumptions, and only
# the second was ever checked. Writing into `$coefficients` reaches `predict()`
# for an `lm` because for an `lm` the coefficients ARE the model. It does not
# reach `predict()` for every object that inherits from `lm`. An `rstanarm`
# fit has class `stanreg, glm, lm` -- deliberately, so that generic code keeps
# working -- and its `$coefficients` are a posterior median, a SUMMARY of the
# draws. `predict()` reads the draws. The assignment is discarded in silence,
# and a bifurcation sweep returns one well formed table per parameter value
# with the same fixed point in every one of them.
#
# Two things follow, and this file holds both.
#
# 1. The substitution is verified, never assumed. Whatever route was used to
#    perform it, the predictions it was supposed to move must actually move.
#    The check is a property of the object's behaviour and not of its class
#    name, so it catches any future class with the same property without
#    anyone having to add it to a list.
#
# 2. An object that carries posterior draws is served through them rather than
#    refused. The draws route is certified against the object's own
#    `predict()` before it is used: the design matrix is built by the model's
#    own terms and the per draw predictions must average back to exactly what
#    `predict()` returns. `predict()` therefore remains the referent it became
#    when the lm path was repaired -- what changed is that it now certifies the
#    evaluation instead of performing it, which is a stronger use of it, not a
#    weaker one.


# A typed condition: the caller can tell "this object cannot be swept" from
# "this argument is malformed" without matching on message text.
ed_stop <- function(subclass, message) {
  stop(structure(
    class = c(subclass, "ed_error", "error", "condition"),
    list(message = message, call = NULL)
  ))
}


# -----------------------------------------------------------------------------
# Posterior draws, recognised by what the object can do and not by its name
# -----------------------------------------------------------------------------

# Returns an S x p numeric matrix of posterior draws whose columns are exactly
# the fitted coefficients, or NULL when the object carries no such thing.
#
# The test is duck typed on purpose. Nothing here names a package: an object
# qualifies when `as.matrix()` yields a numeric matrix of more than one row
# whose column names cover every fitted coefficient. A plain `lm` does not
# qualify -- `as.matrix()` on it returns a one column list matrix with no
# column names -- so the recognition cannot misfire on the ordinary case.
# Recognition only ever selects a BETTER route; when it declines, the caller
# falls back to the assignment route, which is then verified. Safety never
# rests on this function saying yes.
ed_draws_of <- function(equation) {
  cf <- tryCatch(stats::coef(equation), error = function(e) NULL)
  if (is.null(cf) || !is.numeric(cf) || is.null(names(cf)) ||
      !all(nzchar(names(cf)))) {
    return(NULL)
  }
  dm <- tryCatch(suppressWarnings(as.matrix(equation)),
                 error = function(e) NULL)
  if (is.null(dm) || !is.matrix(dm) || !is.numeric(dm) || nrow(dm) < 2L) {
    return(NULL)
  }
  cn <- colnames(dm)
  if (is.null(cn) || !all(names(cf) %in% cn)) return(NULL)
  dm[, names(cf), drop = FALSE]
}


# Evenly spaced thinning of the draws matrix.
#
# `n_draws` is a budget on root finding, not on the posterior: every draw kept
# is a full sweep of the search grid. 200 is the default because the Monte
# Carlo standard error of a posterior median from S near independent draws is
# about 1.2533 * sigma / sqrt(S), which at S = 200 is 0.089 posterior standard
# deviations, and of a 5% or 95% quantile about 0.15 -- both well under the
# resolution at which anyone reads a bifurcation diagram. Raise it when the
# interval endpoints themselves are the quantity of interest.
#
# The indices are evenly spaced rather than random: no seed is needed, the
# result is reproducible, and chains stay balanced because a chain occupies a
# contiguous block of the draws matrix.
ed_thin_draws <- function(B, n_draws) {
  s <- nrow(B)
  if (!is.finite(n_draws) || n_draws >= s) return(B)
  idx <- unique(round(seq(1, s, length.out = as.integer(n_draws))))
  B[idx, , drop = FALSE]
}


# -----------------------------------------------------------------------------
# The design matrix, built by the model's own terms
# -----------------------------------------------------------------------------

# Returns a closure z -> model matrix (length(z) rows, one column per fitted
# coefficient), or NULL when the object does not expose what it takes to build
# one. `predvars`, `xlevels` and `contrasts` are all consulted, so `I()` terms,
# factors, contrasts and data dependent bases such as `poly()` or `ns()` come
# out exactly as they were fitted. This is the same construction `predict.lm()`
# performs internally; it is not a re-derivation of the formula by hand, which
# is the thing that made the lm path wrong before it was repaired.
ed_design_of <- function(equation, variable, exogenous_values) {
  trm <- tryCatch(stats::delete.response(stats::terms(equation)),
                  error = function(e) NULL)
  if (is.null(trm)) return(NULL)
  xlev <- equation$xlevels
  ctr <- equation$contrasts
  function(z) {
    nd <- stats::setNames(data.frame(z = z), variable)
    for (nm in names(exogenous_values)) nd[[nm]] <- exogenous_values[[nm]]
    mf <- stats::model.frame(trm, nd, xlev = xlev)
    if (is.null(ctr)) {
      stats::model.matrix(trm, mf)
    } else {
      stats::model.matrix(trm, mf, contrasts.arg = ctr)
    }
  }
}


# The response scale transform the draws must pass through, so that a per draw
# prediction is the same object `predict(type = "response")` averages.
ed_linkinv_of <- function(equation, is_glm) {
  if (!is_glm) return(identity)
  fam <- tryCatch(stats::family(equation), error = function(e) NULL)
  if (is.null(fam) || is.null(fam$linkinv)) return(NULL)
  fam$linkinv
}


# Certify the hand built route against the object's own predict().
#
# This is the control that can fail for exactly what it watches. If the design
# matrix were built with the wrong contrasts, if the draws were misaligned with
# the coefficient names, if the link were applied on the wrong side of the
# average -- any of these moves the reconstruction away from `predict()`, and
# the reconstruction is refused instead of used. The identity is exact rather
# than approximate: `predict()` on a draws bearing object returns the posterior
# mean of the per draw prediction, and that is what is recomputed here.
ed_certify_design <- function(f_predict, design, B, linkinv, z, equation,
                              tol = 1e-8) {
  built <- tryCatch(design(z), error = function(e) NULL)
  if (is.null(built) || !is.matrix(built) || is.null(colnames(built))) {
    ed_stop("ed_uncertifiable_substitution", paste0(
      "this ", paste(class(equation), collapse = "/"), " object carries ",
      "posterior draws, but its design matrix could not be rebuilt from its ",
      "own terms, so a per draw sweep cannot be certified against predict()"))
  }
  miss <- setdiff(colnames(B), colnames(built))
  if (length(miss)) {
    ed_stop("ed_uncertifiable_substitution", paste0(
      "the posterior draws carry ", paste0("'", miss, "'", collapse = ", "),
      " which the rebuilt design matrix has no column for (it has: ",
      paste(colnames(built), collapse = ", "),
      "), so a per draw sweep cannot be certified against predict()"))
  }
  want <- tryCatch(as.numeric(f_predict(z)), error = function(e) NULL)
  if (is.null(want) || !any(is.finite(want))) {
    ed_stop("ed_uncertifiable_substitution", paste0(
      "predict() returned no finite value anywhere on the search grid, so ",
      "the per draw reconstruction of this ",
      paste(class(equation), collapse = "/"),
      " object has nothing to be certified against"))
  }
  got <- rowMeans(linkinv(built[, colnames(B), drop = FALSE] %*% t(B)))
  ok <- is.finite(want) & is.finite(got)
  scale <- max(1, max(abs(want[ok])))
  err <- max(abs(got[ok] - want[ok]))
  if (!(err <= tol * scale)) {
    ed_stop("ed_uncertifiable_substitution", paste0(
      "the per draw reconstruction of this ",
      paste(class(equation), collapse = "/"),
      " object disagrees with its own predict() by ", format(err, digits = 3),
      " (tolerance ", format(tol * scale, digits = 3),
      "): the draws, the design matrix and the link do not compose back into ",
      "what the object predicts, so the sweep is refused rather than run on ",
      "a reconstruction that is not the fitted model"))
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# The verification: did the substitution actually move anything?
# -----------------------------------------------------------------------------

# `values` are the requested overrides, `fitted_coefs` the values they replace.
# An override that asks for the value already fitted is a no-op by arithmetic
# and cannot be used to probe anything; it is dropped from the probe rather
# than made into a false accusation.
ed_effective_overrides <- function(values, fitted_coefs) {
  keep <- vapply(names(values), function(nm) {
    a <- as.numeric(values[[nm]])
    b <- unname(fitted_coefs[nm])
    !isTRUE(all.equal(a, b, tolerance = .Machine$double.eps^0.5))
  }, logical(1))
  values[keep]
}


# Is the term's design column identically zero over the points that will be
# swept? Then no object, however honest, can respond to a change in its
# coefficient, and the sweep really would be flat. Returns TRUE, FALSE, or NA
# when the design matrix could not be rebuilt to find out.
ed_column_is_inert <- function(design, nms, z) {
  if (is.null(design)) return(NA)
  x <- tryCatch(design(z), error = function(e) NULL)
  if (is.null(x) || !is.matrix(x)) return(NA)
  present <- intersect(nms, colnames(x))
  if (!length(present)) return(NA)
  all(vapply(present, function(nm) all(x[, nm] == 0 | is.na(x[, nm])),
             logical(1)))
}


# A coefficient whose term is identically zero over the search grid cannot
# move anything, whichever route performs the substitution. The sweep would be
# a constant, which is the failure being repaired here, so it is refused with
# its own cause rather than folded into the one below.
ed_stop_inert <- function(nms) {
  ed_stop("ed_substitution_inert", paste0(
    "the coefficient ", paste0("'", nms, "'", collapse = ", "),
    " multiplies a term that is identically zero everywhere on the search ",
    "grid, so varying it cannot move a single prediction: the sweep would ",
    "return identical copies of the same fixed points labelled with ",
    "different parameter values"))
}


# The postcondition of the assignment route.
#
# Three calls, in this order, because the third only means something if the
# first two agree: predict twice unchanged (a prediction that is not
# reproducible cannot testify about anything), then predict once with the
# substitution in place. If the substitution moved no prediction anywhere on
# the search grid, the sweep about to be run is a constant, and saying so is
# the whole point.
ed_verify_substitution <- function(make_f, before, after, values, design,
                                   grid, equation) {
  nms <- names(values)
  base1 <- make_f(before)(grid)
  base2 <- make_f(before)(grid)
  if (!isTRUE(all.equal(base1, base2))) {
    ed_stop("ed_prediction_nondeterministic", paste0(
      "predict() on this ", paste(class(equation), collapse = "/"),
      " object does not return the same values twice for the same input, ",
      "so it cannot be used to track fixed points: two sweeps of the same ",
      "grid would disagree for reasons that have nothing to do with the ",
      "parameter being varied"))
  }
  moved <- make_f(after)(grid)
  ok <- suppressWarnings(any(is.finite(base1) & is.finite(moved) &
                               base1 != moved))
  if (isTRUE(ok)) return(invisible(TRUE))
  if (all(!is.finite(base1))) {
    ed_stop("ed_uncertifiable_substitution", paste0(
      "predict() returned no finite value anywhere on the search grid, so ",
      "whether the substitution of ", paste0("'", nms, "'", collapse = ", "),
      " took effect cannot be established"))
  }
  inert <- ed_column_is_inert(design, nms, grid)
  if (isTRUE(inert)) ed_stop_inert(nms)
  ed_stop("ed_substitution_ignored", paste0(
    "writing ", paste0("'", nms, "'", collapse = ", "),
    " into the coefficients of this ",
    paste(class(equation), collapse = "/"),
    " object left every prediction on the search grid unchanged: its ",
    "predict() does not read $coefficients, so the requested override was ",
    "discarded in silence and the sweep would return identical copies of ",
    "the same fixed points labelled with different parameter values",
    if (identical(inert, NA)) paste0(
      " (the design matrix could not be rebuilt to rule out a term that is ",
      "identically zero over this range)") else ""))
}


# -----------------------------------------------------------------------------
# The root finding core, shared by the point path and every posterior draw
# -----------------------------------------------------------------------------

# Grid search for sign changes, refinement by uniroot, and the central
# difference stability classification, given the function and its values on the
# grid. Extracted so that one draw of a posterior is found exactly the way a
# point estimate is: same grid rule, same de-duplication at half a grid step,
# same marginal threshold. There is no second root finder to drift from this
# one.
#
# The reason for a code is returned as the "ed_reason" attribute instead of
# being announced: a caller sweeping several hundred draws must be able to
# report once rather than once per draw.
ed_fixed_points_from <- function(f, f_vals, grid, range, n_grid) {
  empty <- function(reason) {
    out <- data.frame(fixed_point = numeric(0), stability = character(0),
                      eigenvalue = numeric(0))
    attr(out, "ed_reason") <- reason
    out
  }
  valid_idx <- which(!is.na(f_vals))
  if (length(valid_idx) < 2) return(empty("unevaluable"))

  # A grid point that lands exactly on a root has sign 0. Counting it through
  # diff(sign(...)) != 0 counted the same crossing twice -- once into the zero
  # and once out of it -- and uniroot then refined the same root from both
  # adjacent brackets, so it appeared twice in the output. An exact zero is a
  # root in its own right; a crossing is a strict sign flip between
  # consecutive evaluable grid points.
  s <- sign(f_vals[valid_idx])
  fixed_points <- grid[valid_idx][s == 0]

  strict_flips <- which(s[-length(s)] * s[-1] < 0)
  for (i in strict_flips) {
    idx1 <- valid_idx[i]
    idx2 <- valid_idx[i + 1]
    fp <- tryCatch(stats::uniroot(f, c(grid[idx1], grid[idx2]))$root,
                   error = function(e) NA)
    if (!is.na(fp)) fixed_points <- c(fixed_points, fp)
  }

  # Two detections closer than half a grid step are the same root seen twice:
  # the grid cannot resolve distinct roots below its own spacing.
  if (length(fixed_points) > 1) {
    fixed_points <- sort(fixed_points)
    step <- (range[2] - range[1]) / (n_grid - 1)
    fixed_points <- fixed_points[c(TRUE, diff(fixed_points) > step / 2)]
  }
  if (length(fixed_points) == 0) return(empty("none"))

  eps <- 1e-6
  stability <- character(length(fixed_points))
  eigenvalues <- numeric(length(fixed_points))
  for (i in seq_along(fixed_points)) {
    fp <- fixed_points[i]
    f_plus <- f(fp + eps)
    f_minus <- f(fp - eps)
    if (is.na(f_plus) || is.na(f_minus)) {
      eigenvalues[i] <- NA
      stability[i] <- "unknown"
    } else {
      df_dz <- (f_plus - f_minus) / (2 * eps)
      eigenvalues[i] <- df_dz
      stability[i] <- if (df_dz < -eps) {
        "stable"
      } else if (df_dz > eps) {
        "unstable"
      } else {
        "marginal"
      }
    }
  }
  out <- data.frame(fixed_point = fixed_points, stability = stability,
                    eigenvalue = eigenvalues)
  attr(out, "ed_reason") <- "ok"
  out
}


# -----------------------------------------------------------------------------
# Reducing a posterior of fixed points to branches
# -----------------------------------------------------------------------------

#' Reduce a Posterior of Fixed Points to Branches
#'
#' Turns the per draw output of \code{\link{analyze_fixed_points}} on a
#' posterior-bearing fit into one row per branch, with an interval instead of
#' a number.
#'
#' A posterior sweep produces a set of fixed points per draw, and the size of
#' that set can differ between draws. That is not noise: it is the posterior
#' disagreeing about how many fixed points there are, which near a saddle-node
#' is the interesting fact rather than a nuisance to be averaged away.
#'
#' The reduction is therefore explicitly conditional. The modal count \code{k}
#' is taken -- ties broken toward the smaller count, so the summary never
#' claims more structure than the posterior supports -- the draws attaining
#' \code{k} are kept, their roots are ordered along the line, and the j-th root
#' of each is one branch. Ordering along the line is well defined wherever the
#' count is constant, and stays continuous through a branch crossing because
#' the roots coincide there. \code{prob_n_fixed_points} reports the fraction of
#' draws the summary was computed from, so the interval is never read as
#' unconditional.
#'
#' @param fps Data frame returned by \code{\link{analyze_fixed_points}} whose
#'   \code{"posterior"} attribute is \code{TRUE}; it must carry a \code{draw}
#'   column.
#' @param n_draws Number of draws that were swept, including those that found
#'   no fixed point at all and therefore contribute no row to \code{fps}.
#' @param interval Width of the reported posterior interval.
#'
#' @return Data frame with one row per branch: \code{fixed_point} (the
#'   posterior median), \code{stability}, \code{eigenvalue},
#'   \code{branch}, \code{fixed_point_lower}, \code{fixed_point_upper},
#'   \code{prob_stability}, \code{n_fixed_points},
#'   \code{prob_n_fixed_points} and \code{n_draws}. Zero rows when most draws
#'   found no fixed point; the fraction that did is then in the
#'   \code{"prob_n_fixed_points"} attribute.
#' @export
ed_consensus_fixed_points <- function(fps, n_draws, interval = 0.9) {
  if (!is.data.frame(fps)) {
    stop("'fps' must be the data frame analyze_fixed_points() returned",
         call. = FALSE)
  }
  if (!is.numeric(n_draws) || length(n_draws) != 1L || is.na(n_draws) ||
      n_draws < 1) {
    stop("'n_draws' must be the number of draws that were swept, at least 1",
         call. = FALSE)
  }
  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval) ||
      interval <= 0 || interval >= 1) {
    stop("'interval' must be a single number strictly between 0 and 1",
         call. = FALSE)
  }
  if (nrow(fps) > 0L && is.null(fps$draw)) {
    # Without the draw column every row would be counted into one draw, and
    # the reduction would return a confident summary of an arithmetic error.
    stop("'fps' has no 'draw' column: it did not come from a sweep over ",
         "posterior draws (its \"posterior\" attribute would be TRUE)",
         call. = FALSE)
  }
  n_draws <- as.integer(n_draws)
  counts <- integer(n_draws)
  if (nrow(fps) > 0L) {
    tb <- table(factor(fps$draw, levels = seq_len(n_draws)))
    counts <- as.integer(tb)
  }
  tk <- table(counts)
  k <- as.integer(names(tk)[which.max(as.integer(tk))])
  prob_k <- mean(counts == k)

  empty <- data.frame(
    fixed_point = numeric(0), stability = character(0),
    eigenvalue = numeric(0), branch = integer(0),
    fixed_point_lower = numeric(0), fixed_point_upper = numeric(0),
    prob_stability = numeric(0), n_fixed_points = integer(0),
    prob_n_fixed_points = numeric(0), n_draws = integer(0),
    stringsAsFactors = FALSE)
  if (k == 0L) {
    attr(empty, "prob_n_fixed_points") <- prob_k
    attr(empty, "n_draws") <- n_draws
    return(empty)
  }

  keep <- which(counts == k)
  sub <- fps[fps$draw %in% keep, , drop = FALSE]
  sub <- sub[order(sub$draw, sub$fixed_point), , drop = FALSE]
  m_fp <- matrix(sub$fixed_point, nrow = k)
  m_st <- matrix(sub$stability, nrow = k)
  m_ev <- matrix(sub$eigenvalue, nrow = k)

  lo <- (1 - interval) / 2
  qs <- c(lo, 0.5, 1 - lo)
  out <- do.call(rbind, lapply(seq_len(k), function(j) {
    q <- unname(stats::quantile(m_fp[j, ], probs = qs, na.rm = TRUE))
    tt <- table(m_st[j, ])
    # Deterministic tie break: most frequent label, then alphabetical.
    st <- names(tt)[order(-as.integer(tt), names(tt))][1]
    data.frame(
      fixed_point = q[2], stability = st,
      eigenvalue = stats::median(m_ev[j, ], na.rm = TRUE),
      branch = j,
      fixed_point_lower = q[1], fixed_point_upper = q[3],
      prob_stability = mean(m_st[j, ] == st),
      n_fixed_points = k, prob_n_fixed_points = prob_k,
      n_draws = n_draws, stringsAsFactors = FALSE)
  }))
  attr(out, "prob_n_fixed_points") <- prob_k
  attr(out, "n_draws") <- n_draws
  attr(out, "interval") <- interval
  out
}
