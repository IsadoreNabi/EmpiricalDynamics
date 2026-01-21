#' @title Preprocessing Functions for Time Series Data
#' @description Functions for data preparation, variable specification, and 
#'   numerical differentiation including Total Variation Regularized (TVR)
#'   differentiation for noisy economic data.
#' @name preprocessing
NULL


#' Specify Variable Types for Dynamical Analysis
#'
#' Classifies variables in a dataset according to their role in the dynamical
#' system being studied.
#'
#' @param data A data.frame containing the time series data.
#' @param endogenous Character vector of endogenous state variable names.
#' @param endogenous_coupled Character vector of coupled endogenous variables.
#' @param coupled Alias for endogenous_coupled (for compatibility).
#' @param exogenous Character vector of exogenous forcing variable names.
#' @param slow_parameter Character vector of slowly-varying parameter names.
#' @param time Name of the time column (auto-detected if NULL).
#' @param time_col Alias for time (for compatibility).
#'
#' @return The input data.frame with variable specifications added as the 
#'   "var_spec" attribute.
#'
#' @details
#' Variable types:
#' \itemize{
#'   \item \strong{endogenous}: Variables whose dynamics are modeled (appear as dZ/dt).
#'   \item \strong{endogenous_coupled}: Variables that co-evolve with endogenous vars.
#'   \item \strong{exogenous}: Variables that influence the system but are not modeled.
#'   \item \strong{slow_parameter}: Variables that change on much longer timescales.
#' }
#'
#' @examples
#' data <- data.frame(
#'   time = 1:10,
#'   profit_rate = runif(10),
#'   capital_stock = runif(10),
#'   interest_rate = runif(10)
#' )
#' 
#' data <- specify_variables(data,
#'   endogenous = "profit_rate",
#'   endogenous_coupled = "capital_stock",
#'   exogenous = "interest_rate"
#' )
#' attr(data, "var_spec")
#'
#' @export
specify_variables <- function(data, 
                              endogenous = NULL,
                              endogenous_coupled = NULL,
                              coupled = NULL,
                              exogenous = NULL,
                              slow_parameter = NULL,
                              time = NULL,
                              time_col = NULL) {
  
  # Alias for compatibility
  if (is.null(endogenous_coupled) && !is.null(coupled)) {
    endogenous_coupled <- coupled
  }
  if (is.null(time) && !is.null(time_col)) {
    time <- time_col
  }
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame", call. = FALSE)
  }
  
  all_vars <- c(endogenous, endogenous_coupled, exogenous, slow_parameter)
  # Filter out NULLs before checking names
  all_vars <- all_vars[!is.null(all_vars)]
  
  missing <- setdiff(all_vars, names(data))
  
  if (length(missing) > 0) {
    stop("Variables not found in data: ", paste(missing, collapse = ", "), 
         call. = FALSE)
  }
  
  # Auto-detect time column
  if (is.null(time)) {
    time_candidates <- c("time", "t", "date", "year", "period", "Time", "Date", "Year")
    time <- intersect(time_candidates, names(data))[1]
  }
  
  # Create variable specification
  var_spec <- list(
    endogenous = endogenous,
    endogenous_coupled = endogenous_coupled,
    coupled = endogenous_coupled,
    exogenous = exogenous,
    slow_parameter = slow_parameter,
    time = time,
    time_col = time,
    all_modeled = c(endogenous, endogenous_coupled),
    all_predictors = c(endogenous, endogenous_coupled, exogenous, slow_parameter)
  )
  
  attr(data, "var_spec") <- var_spec
  class(data) <- c("specified_data", class(data))
  
  message("Variable specification:")
  if (!is.null(endogenous)) message("  Endogenous: ", paste(endogenous, collapse = ", "))
  if (!is.null(endogenous_coupled)) message("  Coupled: ", paste(endogenous_coupled, collapse = ", "))
  if (!is.null(exogenous)) message("  Exogenous: ", paste(exogenous, collapse = ", "))
  if (!is.null(time)) message("  Time: ", time)
  
  data
}


#' Compute Derivative of a Time Series
#'
#' Main dispatcher function for numerical differentiation. Supports multiple
#' methods appropriate for different data characteristics.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points (NULL assumes dt=1).
#' @param method Differentiation method: "tvr", "savgol", "spline", 
#'   "finite_diff", "fd", or "spectral".
#' @param ... Additional arguments passed to the specific method.
#'
#' @return Numeric vector of estimated derivatives with diagnostic attributes.
#'
#' @details
#' Available methods:
#' \itemize{
#'   \item \strong{tvr}: Total Variation Regularized differentiation (recommended 
#'     for economic data with trends and shocks).
#'   \item \strong{savgol}: Savitzky-Golay filter (moderate noise, preserves peaks).
#'   \item \strong{spline}: Smoothing spline (high noise, prioritizes trend).
#'   \item \strong{finite_diff} or \strong{fd}: Centered finite differences (low noise).
#'   \item \strong{spectral}: FFT-based (periodic data only).
#' }
#'
#' @seealso \code{\link{compute_derivative_tvr}}, \code{\link{suggest_differentiation_method}}
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' Z <- sin(t) + rnorm(length(t), sd = 0.1)
#' 
#' # Finite differences (fast, no dependencies)
#' dZ_fd <- compute_derivative(Z, t, method = "finite_diff")
#' 
#' # Access the derivative vector for plotting
#' plot(t, dZ_fd$derivative, type = "l", main = "Derivative Comparison")
#' lines(t, cos(t), col = "red", lty = 2) # True derivative
#' 
#' \donttest{
#' # TVR (requires CVXR)
#' if (requireNamespace("CVXR", quietly = TRUE)) {
#'   dZ_tvr <- compute_derivative(Z, t, method = "tvr")
#' }
#' }
#'
#' @export
compute_derivative <- function(Z, t = NULL, 
                               method = c("tvr", "savgol", "spline", 
                                          "finite_diff", "spectral"),
                               ...) {
  # Alias for compatibility: "fd" -> "finite_diff"
  if (!missing(method) && length(method) == 1 && identical(method, "fd")) {
    method <- "finite_diff"
  }
  
  method <- match.arg(method)
  
  # Validate inputs
  if (exists("validate_timeseries", mode = "function")) {
    validate_timeseries(Z, "Z", allow_na = FALSE, min_length = 5)
  } else {
    if (length(Z) < 5) stop("'Z' must have at least 5 observations.")
    if (any(is.na(Z))) stop("'Z' contains NA values.")
  }
  
  if (is.null(t)) {
    t <- seq_along(Z)
  } else if (length(t) != length(Z)) {
    stop("'t' and 'Z' must have the same length", call. = FALSE)
  }
  
  # Dispatch to appropriate method
  result <- switch(method,
                   "tvr"         = compute_derivative_tvr(Z, t, ...),
                   "savgol"      = compute_derivative_savgol(Z, t, ...),
                   "spline"      = compute_derivative_spline(Z, t, ...),
                   "finite_diff" = compute_derivative_fd(Z, t, ...),
                   "spectral"    = compute_derivative_spectral(Z, t, ...)
  )
  
  # Add common attributes
  attr(result, "method") <- method
  attr(result, "n") <- length(Z)
  attr(result, "t") <- t
  
  # Ensure result is a list with derivative element for compatibility
  if (!is.list(result)) {
    result <- list(derivative = as.numeric(result))
    attr(result$derivative, "method") <- method
  }
  
  result
}


#' Total Variation Regularized Differentiation
#'
#' Computes derivatives by solving a convex optimization problem that balances
#' fidelity to the data against smoothness of the derivative.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points (NULL assumes dt=1).
#' @param lambda Regularization parameter ("auto" for cross-validation selection).
#' @param solver Optimization backend: "osqp", "ecos", or "scs".
#' @param ... Additional arguments (ignored).
#'
#' @return Object of class "tvr_derivative" (also a list with $derivative).
#'
#' @export
compute_derivative_tvr <- function(Z, t = NULL, lambda = "auto", 
                                   solver = "osqp", ...) {
  
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for TVR differentiation. Please install it.")
  }
  
  n <- length(Z)
  if (is.null(t)) t <- 1:n
  dt <- diff(t)
  
  if (!exists("build_integration_matrix", mode = "function")) {
    stop("Internal function 'build_integration_matrix' not found.")
  }
  
  A <- build_integration_matrix(n, dt)
  D <- build_difference_matrix(n)
  
  # Define optimization problem
  dZ <- CVXR::Variable(n)
  Z_centered <- Z - Z[1]
  
  fidelity <- CVXR::sum_squares(Z_centered - A %*% dZ)
  total_variation <- CVXR::norm1(D %*% dZ)
  
  # Auto-select lambda if requested
  if (identical(lambda, "auto")) {
    lambda <- select_lambda_cv_tvr(Z, t, A, D, solver)
    message(sprintf("Lambda selected automatically: %.4e", lambda))
  }
  
  # Solve optimization problem
  objective <- CVXR::Minimize(fidelity + lambda * total_variation)
  problem <- CVXR::Problem(objective)
  
  result <- tryCatch(
    CVXR::solve(problem, solver = toupper(solver)),
    error = function(e) {
      # Try alternative solver
      message("Primary solver failed, trying SCS...")
      tryCatch(
        CVXR::solve(problem, solver = "SCS"),
        error = function(e2) list(status = "error")
      )
    }
  )
  
  if (is.null(result$status) || !result$status %in% c("optimal", "optimal_inaccurate")) {
    warning(paste("Solver did not reach optimum. Status:", result$status))
  }
  
  dZ_hat <- as.vector(result$getValue(dZ))
  
  # Compute reconstruction for diagnostics
  Z_reconstructed <- Z[1] + c(0, cumsum(dZ_hat[-n] * dt))
  reconstruction_error <- sqrt(mean((Z - Z_reconstructed)^2))
  
  # Return with attributes
  out <- structure(
    dZ_hat,
    lambda = lambda,
    solver_status = result$status,
    fidelity_term = as.numeric(result$getValue(fidelity)),
    tv_term = as.numeric(result$getValue(total_variation)),
    reconstruction_rmse = reconstruction_error,
    Z_reconstructed = Z_reconstructed,
    class = c("tvr_derivative", "numeric")
  )
  
  # Also return as list for compatibility, AND assign S3 class to list for tests
  res_list <- list(
    derivative = as.numeric(out),
    lambda = lambda,
    solver_status = result$status,
    reconstruction_rmse = reconstruction_error,
    Z_reconstructed = Z_reconstructed,
    raw = out
  )
  class(res_list) <- c("tvr_derivative", "list")
  
  res_list
}


#' Cross-Validation Selection of Lambda for TVR
#'
#' Selects the regularization parameter lambda using leave-one-out-like
#' cross-validation.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param A Integration matrix.
#' @param D Difference matrix.
#' @param solver Optimization backend.
#' @param lambda_seq Sequence of lambda values to evaluate.
#' @param verbose Print progress?
#'
#' @return Selected lambda value.
#'
#' @export
select_lambda_cv_tvr <- function(Z, t, A = NULL, D = NULL, solver = "osqp",
                                 lambda_seq = 10^seq(-4, 2, length.out = 30),
                                 verbose = FALSE) {
  
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("CVXR required")
  }
  
  n <- length(Z)
  dt <- diff(t)
  Z_centered <- Z - Z[1]
  
  if (is.null(A)) A <- build_integration_matrix(n, dt)
  if (is.null(D)) D <- build_difference_matrix(n)
  
  cv_errors <- sapply(lambda_seq, function(lam) {
    if (verbose) cat(".")
    
    dZ <- CVXR::Variable(n)
    fidelity <- CVXR::sum_squares(Z_centered - A %*% dZ)
    tv <- CVXR::norm1(D %*% dZ)
    prob <- CVXR::Problem(CVXR::Minimize(fidelity + lam * tv))
    
    res <- tryCatch(
      CVXR::solve(prob, solver = toupper(solver)),
      error = function(e) list(status = "error")
    )
    
    if (!is.null(res$status) && res$status %in% c("optimal", "optimal_inaccurate")) {
      dZ_hat <- as.vector(res$getValue(dZ))
      Z_reconstructed <- Z[1] + c(0, cumsum(dZ_hat[-n] * dt))
      mean((Z - Z_reconstructed)^2)
    } else {
      Inf
    }
  })
  
  if (verbose) cat("\n")
  
  # One-standard-error rule
  valid_idx <- is.finite(cv_errors)
  if (sum(valid_idx) < 3) {
    warning("Few valid lambdas. Using minimum directly.")
    return(lambda_seq[which.min(cv_errors)])
  }
  
  min_idx <- which.min(cv_errors)
  se <- stats::sd(cv_errors[valid_idx]) / sqrt(sum(valid_idx))
  threshold <- cv_errors[min_idx] + se
  
  # Largest lambda within threshold (more regularization = simpler model)
  candidates <- which(cv_errors <= threshold & seq_along(cv_errors) >= min_idx)
  lambda_seq[max(candidates)]
}


#' Savitzky-Golay Derivative
#'
#' Computes derivatives using the Savitzky-Golay filter.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param p Polynomial order (default 3).
#' @param n Filter length (must be odd, default auto-selected).
#' @param m Derivative order (default 1).
#' @param ... Additional arguments (ignored).
#'
#' @return List with derivative vector.
#'
#' @export
compute_derivative_savgol <- function(Z, t = NULL, p = 3, n = NULL, m = 1, ...) {
  
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop("Package 'signal' is required for Savitzky-Golay differentiation.")
  }
  
  N <- length(Z)
  if (is.null(t)) t <- 1:N
  dt <- mean(diff(t))
  
  # Auto-select filter length
  if (is.null(n)) {
    n <- min(2 * floor(N / 10) + 1, 21)
    n <- max(n, p + 2)
    if (n %% 2 == 0) n <- n + 1
  }
  
  # Ensure valid parameters
  if (n > N) n <- N - (N %% 2 == 0)
  if (p >= n) p <- n - 1
  
  # Apply filter
  dZ <- signal::sgolayfilt(Z, p = p, n = n, m = m) / (dt^m)
  
  list(
    derivative = as.numeric(dZ),
    filter_order = p,
    filter_length = n,
    dt = dt
  )
}


#' Smoothing Spline Derivative
#'
#' Computes derivatives by fitting a smoothing spline and differentiating it.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param spar Smoothing parameter (NULL for automatic selection).
#' @param df Degrees of freedom (alternative to spar).
#' @param ... Additional arguments (ignored).
#'
#' @return List with derivative vector.
#'
#' @export
compute_derivative_spline <- function(Z, t = NULL, spar = NULL, df = NULL, ...) {
  
  N <- length(Z)
  if (is.null(t)) t <- 1:N
  
  # Fit smoothing spline
  if (!is.null(df)) {
    sp <- stats::smooth.spline(t, Z, df = df)
  } else if (!is.null(spar)) {
    sp <- stats::smooth.spline(t, Z, spar = spar)
  } else {
    sp <- stats::smooth.spline(t, Z)
  }
  
  # Predict derivative
  dZ <- stats::predict(sp, t, deriv = 1)$y
  
  list(
    derivative = as.numeric(dZ),
    spar = sp$spar,
    df = sp$df,
    cv_error = sp$cv.crit
  )
}


#' Centered Finite Differences
#'
#' Computes derivatives using centered finite differences.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param ... Additional arguments (ignored).
#'
#' @return List with derivative vector.
#'
#' @export
compute_derivative_fd <- function(Z, t = NULL, ...) {
  
  N <- length(Z)
  if (is.null(t)) t <- 1:N
  
  dZ <- numeric(N)
  
  # Forward difference at first point
  dZ[1] <- (Z[2] - Z[1]) / (t[2] - t[1])
  
  # Centered differences in interior
  for (i in 2:(N-1)) {
    dZ[i] <- (Z[i+1] - Z[i-1]) / (t[i+1] - t[i-1])
  }
  
  # Backward difference at last point
  dZ[N] <- (Z[N] - Z[N-1]) / (t[N] - t[N-1])
  
  list(derivative = as.numeric(dZ))
}


#' Spectral (FFT) Differentiation
#'
#' Computes derivatives using the Fourier transform.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param ... Additional arguments (ignored).
#'
#' @return List with derivative vector.
#'
#' @section Warning:
#' This method assumes the signal is periodic.
#'
#' @export
compute_derivative_spectral <- function(Z, t = NULL, ...) {
  
  warning("Spectral differentiation assumes periodic data. ",
          "For economic time series, consider using 'tvr' method instead.")
  
  N <- length(Z)
  if (is.null(t)) t <- 1:N
  
  dt <- mean(diff(t))
  
  # Compute frequency vector
  if (N %% 2 == 0) {
    k <- c(0:(N/2-1), 0, (-(N/2-1)):-1)
  } else {
    k <- c(0:((N-1)/2), (-((N-1)/2)):-1)
  }
  
  omega <- 2 * pi * k / (N * dt)
  
  # Differentiate in frequency domain
  Z_fft <- stats::fft(Z)
  dZ_fft <- 1i * omega * Z_fft
  dZ <- Re(stats::fft(dZ_fft, inverse = TRUE)) / N
  
  list(
    derivative = as.numeric(dZ),
    dt = dt
  )
}


#' Suggest Differentiation Method Based on Data Characteristics
#'
#' Analyzes the time series to recommend the most appropriate differentiation
#' method based on detected features like trend, periodicity, shocks, and noise.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#'
#' @return List with suggested method and diagnostic information.
#'
#' @examples
#' t <- 1:100
#' Z <- 0.1 * t + rnorm(100) # Trend with noise
#' result <- suggest_differentiation_method(Z, t)
#' print(result$suggested_method)
#'
#' @export
suggest_differentiation_method <- function(Z, t = NULL) {
  
  if (is.null(t)) t <- seq_along(Z)
  n <- length(Z)
  
  diagnostics <- list()
  
  # 1. Detect trend (simplified Mann-Kendall)
  dZ_simple <- diff(Z)
  trend_sign <- sign(outer(Z, Z, "-"))
  S <- sum(trend_sign[lower.tri(trend_sign)])
  diagnostics$trend_statistic <- S / (n * (n-1) / 2)
  diagnostics$has_trend <- abs(diagnostics$trend_statistic) > 0.3
  
  # 2. Detect periodicity
  max_lag <- min(floor(n/3), 50)
  if (max_lag > 1) {
    acf_vals <- stats::acf(Z, lag.max = max_lag, plot = FALSE)$acf[-1]
    diagnostics$max_acf <- max(abs(acf_vals))
    diagnostics$is_periodic <- any(abs(acf_vals[min(10, length(acf_vals)):length(acf_vals)]) > 0.5)
  } else {
    diagnostics$is_periodic <- FALSE
    diagnostics$max_acf <- 0
  }
  
  # 3. Detect shocks/discontinuities
  mad_dZ <- stats::mad(dZ_simple)
  if (mad_dZ > 0) {
    outliers <- abs(dZ_simple - stats::median(dZ_simple)) > 4 * mad_dZ
    diagnostics$n_shocks <- sum(outliers)
  } else {
    diagnostics$n_shocks <- 0
  }
  diagnostics$has_shocks <- diagnostics$n_shocks > 0
  
  # 4. Estimate noise level
  d2Z <- diff(dZ_simple)
  sd_dZ <- stats::sd(dZ_simple)
  if (sd_dZ > 0) {
    diagnostics$noise_level <- stats::sd(d2Z) / sd_dZ
  } else {
    diagnostics$noise_level <- 0
  }
  diagnostics$is_noisy <- diagnostics$noise_level > 0.5
  
  # Decision tree
  if (diagnostics$has_shocks || diagnostics$has_trend) {
    method <- "tvr"
    reason <- "Series has trend and/or shocks detected"
  } else if (diagnostics$is_periodic && !diagnostics$is_noisy) {
    method <- "savgol"
    reason <- "Periodic series with moderate noise"
  } else if (diagnostics$is_noisy) {
    method <- "spline"
    reason <- "High noise level detected"
  } else {
    method <- "finite_diff"
    reason <- "Smooth series without special characteristics"
  }
  
  diagnostics$suggested_method <- method
  diagnostics$recommended <- method
  diagnostics$reason <- reason
  
  message(sprintf("Suggested method: %s\nReason: %s", method, reason))
  
  diagnostics
}


#' Compare Differentiation Methods
#'
#' Applies multiple differentiation methods to the same data and produces
#' a comparison plot.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#' @param methods Character vector of methods to compare.
#' @param plot Produce comparison plot?
#'
#' @return List of derivative vectors from each method.
#'
#' @export
compare_differentiation_methods <- function(Z, t = NULL, 
                                            methods = c("tvr", "savgol", "spline", "finite_diff"),
                                            plot = TRUE) {
  
  if (is.null(t)) t <- seq_along(Z)
  
  results <- lapply(methods, function(m) {
    tryCatch({
      res <- compute_derivative(Z, t, method = m)
      if (is.list(res) && "derivative" %in% names(res)) {
        res$derivative
      } else {
        as.numeric(res)
      }
    }, error = function(e) {
      warning(sprintf("Method '%s' failed: %s", m, e$message))
      rep(NA, length(Z))
    })
  })
  names(results) <- methods
  
  if (plot) {
    colors <- c("red", "blue", "green3", "orange", "purple")
    y_range <- range(unlist(results), na.rm = TRUE)
    
    graphics::plot(t, results[[1]], type = "n", 
                   ylim = y_range,
                   main = "Comparison of Differentiation Methods",
                   xlab = "t", ylab = expression(dot(Z)))
    
    for (i in seq_along(methods)) {
      if (!all(is.na(results[[i]]))) {
        graphics::lines(t, results[[i]], col = colors[i], lwd = 1.5)
      }
    }
    
    graphics::legend("topright", methods, col = colors[seq_along(methods)], 
                     lwd = 1.5, cex = 0.8, bg = "white")
    graphics::grid()
  }
  
  invisible(results)
}


#' Diagnostic Plot for TVR Differentiation
#'
#' Produces a four-panel diagnostic plot showing the original series,
#' estimated derivative, reconstruction comparison, and residuals.
#'
#' @param Z Numeric vector of original observations.
#' @param t Numeric vector of time points.
#' @param dZ_tvr TVR derivative object (from \code{compute_derivative_tvr}).
#'
#' @return Invisibly returns a list of diagnostic values.
#'
#' @export
plot_tvr_diagnostic <- function(Z, t = NULL, dZ_tvr) {
  
  if (is.null(t)) t <- seq_along(Z)
  
  # Handle both list and raw vector formats
  if (is.list(dZ_tvr)) {
    dZ_vals <- dZ_tvr$derivative
    Z_reconstructed <- dZ_tvr$Z_reconstructed
    lambda <- dZ_tvr$lambda
    rmse <- dZ_tvr$reconstruction_rmse
  } else {
    dZ_vals <- as.numeric(dZ_tvr)
    Z_reconstructed <- attr(dZ_tvr, "Z_reconstructed")
    lambda <- attr(dZ_tvr, "lambda")
    rmse <- attr(dZ_tvr, "reconstruction_rmse")
  }
  
  if (is.null(Z_reconstructed)) {
    dt <- diff(t)
    n <- length(Z)
    Z_reconstructed <- Z[1] + c(0, cumsum(dZ_vals[-n] * dt))
  }
  
  # Set up plot
  old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par))
  
  # 1. Original series
  graphics::plot(t, Z, type = "l", main = "Original Series Z(t)", 
                 xlab = "t", ylab = "Z", col = "darkblue", lwd = 1.5)
  graphics::grid()
  
  # 2. Estimated derivative
  graphics::plot(t, dZ_vals, type = "l", main = "Estimated Derivative (TVR)",
                 xlab = "t", ylab = expression(dot(Z)), col = "darkred", lwd = 1.5)
  graphics::abline(h = 0, col = "gray", lty = 2)
  graphics::grid()
  if (!is.null(lambda)) {
    graphics::mtext(sprintf("lambda = %.2e", lambda), side = 3, line = -1.5, 
                    adj = 0.95, cex = 0.8)
  }
  
  # 3. Reconstruction vs original
  graphics::plot(t, Z, type = "l", main = "Reconstruction vs Original", 
                 xlab = "t", ylab = "Z", col = "darkblue", lwd = 1.5)
  graphics::lines(t, Z_reconstructed, col = "red", lty = 2, lwd = 1.5)
  graphics::legend("topleft", c("Original", "Reconstructed"), 
                   col = c("darkblue", "red"), lty = c(1, 2), lwd = 1.5, cex = 0.8)
  graphics::grid()
  if (!is.null(rmse)) {
    graphics::mtext(sprintf("RMSE = %.4f", rmse), side = 3, line = -1.5, 
                    adj = 0.95, cex = 0.8)
  }
  
  # 4. Reconstruction residual
  residual <- Z - Z_reconstructed
  graphics::plot(t, residual, type = "l", 
                 main = "Reconstruction Residual", xlab = "t", ylab = "Residual",
                 col = "darkgreen", lwd = 1.5)
  graphics::abline(h = 0, col = "gray", lty = 2)
  graphics::grid()
  
  invisible(list(
    Z_reconstructed = Z_reconstructed,
    residuals = residual,
    lambda = lambda,
    rmse = rmse
  ))
}


#' Print Method for TVR Derivative
#'
#' @param x A tvr_derivative object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
print.tvr_derivative <- function(x, ...) {
  cat("TVR Derivative Estimate\n")
  cat("=======================\n")
  cat(sprintf("  Length: %d\n", length(x$derivative)))
  cat(sprintf("  Lambda: %.4e\n", x$lambda))
  cat(sprintf("  Solver status: %s\n", x$solver_status))
  cat(sprintf("  Reconstruction RMSE: %.6f\n", x$reconstruction_rmse))
  cat(sprintf("  Range: [%.4f, %.4f]\n", min(x$derivative), max(x$derivative)))
  invisible(x)
}


#' Plot Method for TVR Derivative
#'
#' @param x A tvr_derivative object.
#' @param t Time vector (uses attribute if NULL).
#' @param ... Additional plot arguments.
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
plot.tvr_derivative <- function(x, t = NULL, ...) {
  if (is.null(t)) {
    t <- attr(x, "t")
    if (is.null(t)) t <- seq_along(x$derivative)
  }
  
  graphics::plot(t, x$derivative, type = "l", 
                 main = "TVR Derivative Estimate",
                 xlab = "t", ylab = expression(dot(Z)),
                 col = "darkred", lwd = 1.5, ...)
  graphics::abline(h = 0, col = "gray", lty = 2)
  graphics::grid()
  
  if (!is.null(x$lambda)) {
    graphics::mtext(sprintf("lambda = %.2e", x$lambda), 
                    side = 3, line = -1.5, adj = 0.95, cex = 0.8)
  }
  
  invisible(x)
}


#' Diagnose Sampling Frequency
#'
#' Evaluates whether the sampling frequency is appropriate for capturing
#' the dynamics of the phenomenon.
#'
#' @param Z Numeric vector of observations.
#' @param t Numeric vector of time points.
#'
#' @return List with diagnostic information and recommendations.
#'
#' @export
diagnose_sampling_frequency <- function(Z, t = NULL) {
  
  if (is.null(t)) t <- seq_along(Z)
  n <- length(Z)
  
  dZ <- diff(Z)
  dt <- diff(t)
  
  diagnostics <- list()
  
  # Check for uniform sampling
  dt_var <- stats::var(dt) / mean(dt)^2
  diagnostics$uniform_sampling <- dt_var < 0.01
  diagnostics$mean_dt <- mean(dt)
  diagnostics$dt_mean <- mean(dt)  # Alias
  diagnostics$dt_cv <- sqrt(dt_var)
  
  # Autocorrelation of differences
  acf_dZ <- stats::acf(dZ, lag.max = min(20, floor(n/4)), plot = FALSE)$acf[-1]
  diagnostics$dZ_autocorr <- acf_dZ
  
  # Slow decay suggests oversampling
  if (length(acf_dZ) >= 5) {
    decay_rate <- abs(acf_dZ[5] / acf_dZ[1])
    diagnostics$oversampled <- decay_rate > 0.8 && acf_dZ[1] > 0.5
  } else {
    diagnostics$oversampled <- FALSE
  }
  
  # Detect oscillations and check Nyquist
  spec <- stats::spectrum(Z, plot = FALSE)
  peak_freq <- spec$freq[which.max(spec$spec)]
  diagnostics$dominant_frequency <- peak_freq
  diagnostics$samples_per_cycle <- 1 / (peak_freq * mean(dt))
  diagnostics$nyquist_ok <- diagnostics$samples_per_cycle >= 2
  diagnostics$aliasing_risk <- diagnostics$samples_per_cycle < 4
  
  # Generate recommendation
  diagnostics$recommendation <- if (diagnostics$oversampled) {
    "Consider subsampling data"
  } else if (diagnostics$aliasing_risk) {
    "Sample more frequently to avoid aliasing"
  } else {
    "Sampling frequency appears adequate"
  }
  diagnostics$message <- diagnostics$recommendation  # Alias
  
  # Generate warnings
  if (diagnostics$oversampled) {
    warning("Data may be oversampled: autocorrelation of differences decays slowly")
  }
  
  if (diagnostics$aliasing_risk) {
    warning(sprintf("Aliasing risk: only %.1f samples per dominant cycle (recommend >= 4)",
                    diagnostics$samples_per_cycle))
  }
  
  diagnostics
}


#' Compute Derivatives for Specified Variables
#'
#' Convenience function to compute derivatives for all endogenous variables
#' in a specified dataset.
#'
#' @param data Data frame with variable specifications (from \code{specify_variables}).
#' @param method Differentiation method.
#' @param prefix Prefix for derivative column names (default "d_").
#' @param ... Additional arguments passed to \code{compute_derivative}.
#'
#' @return Data frame with derivative columns added.
#'
#' @export
compute_derivatives <- function(data, method = "tvr", prefix = "d_", ...) {
  
  var_spec <- attr(data, "var_spec")
  
  if (is.null(var_spec)) {
    stop("Data must be processed with specify_variables() first", call. = FALSE)
  }
  
  # Get time vector
  t <- if (!is.null(var_spec$time)) data[[var_spec$time]] else NULL
  
  # Compute derivatives for all modeled variables
  for (varname in var_spec$all_modeled) {
    deriv_name <- paste0(prefix, varname)
    message(sprintf("Computing derivative: %s -> %s", varname, deriv_name))
    
    result <- compute_derivative(data[[varname]], t, method = method, ...)
    if (is.list(result) && "derivative" %in% names(result)) {
      data[[deriv_name]] <- result$derivative
    } else {
      data[[deriv_name]] <- as.numeric(result)
    }
  }
  
  data
}


#' Create Candidate Transformations
#'
#' Generates derived variables based on specified transformations that may
#' be theoretically relevant.
#'
#' @param data Data frame.
#' @param transformations List of formulas specifying transformations.
#' @param variables Character vector of variable names (alternative interface).
#' @param cols Alias for variables.
#'
#' @return Data frame with transformation columns added.
#'
#' @examples
#' data <- data.frame(X = 1:10, Y = 10:1)
#' 
#' # Simple interface
#' data <- create_transformations(data, variables = c("X", "Y"))
#' 
#' # Formula interface
#' data <- create_transformations(data, transformations = list(
#'   ratios = ~ X/Y
#' ))
#'
#' @export
create_transformations <- function(data, transformations = NULL, 
                                   variables = NULL, cols = NULL) {
  
  # If variables/cols are passed instead of transformations formula
  if (is.null(transformations) && (!is.null(variables) || !is.null(cols))) {
    vars <- if (!is.null(variables)) variables else cols
    result <- data
    
    for (v in vars) {
      if (v %in% names(data)) {
        result[[paste0("log_", v)]] <- log(abs(data[[v]]) + 1e-10)
        result[[paste0("sqrt_", v)]] <- sqrt(abs(data[[v]]))
        result[[paste0("sq_", v)]] <- data[[v]]^2
      }
    }
    return(result)
  }
  
  # Formula interface
  if (is.null(transformations)) {
    return(data)
  }
  
  for (name in names(transformations)) {
    formula <- transformations[[name]]
    terms <- all.vars(formula)
    
    # Check all variables exist
    missing <- setdiff(terms, names(data))
    if (length(missing) > 0) {
      warning(sprintf("Skipping '%s': missing variables %s", 
                      name, paste(missing, collapse = ", ")))
      next
    }
    
    # Parse and evaluate each term
    expr_text <- as.character(formula)[2]
    term_exprs <- strsplit(expr_text, "\\+")[[1]]
    term_exprs <- trimws(term_exprs)
    
    for (i in seq_along(term_exprs)) {
      expr <- term_exprs[i]
      col_name <- if (length(term_exprs) == 1) name else paste0(name, "_", i)
      
      tryCatch({
        data[[col_name]] <- eval(parse(text = expr), envir = data)
        message(sprintf("Created: %s = %s", col_name, expr))
      }, error = function(e) {
        warning(sprintf("Failed to create '%s': %s", col_name, e$message))
      })
    }
  }
  
  data
}