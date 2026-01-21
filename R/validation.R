# =============================================================================
# EmpiricalDynamics: validation.R
# Step E: Validation - Cross-validation, trajectory simulation, qualitative tests
# =============================================================================

#' @title Validation of Discovered Equations
#' @description Functions for validating discovered differential equations through
#'   cross-validation, trajectory simulation, and qualitative behavior analysis.
#' @name validation
NULL

# =============================================================================
# CROSS-VALIDATION
# =============================================================================

#' Cross-Validate Discovered Equation
#'
#' Performs k-fold or block cross-validation to assess out-of-sample predictive
#' performance of the discovered equation.
#'
#' @param equation Fitted equation object from \code{fit_specified_equation} or
#'   \code{symbolic_search}. Can also be an object of class \code{lm} or \code{nls}.
#' @param data Data frame containing all variables.
#' @param response Name of the response column (derivative).
#' @param derivative_col Alias for response (for compatibility).
#' @param k Number of folds for cross-validation.
#' @param method CV method: "random", "block", "rolling".
#' @param block_size For block methods, size of contiguous blocks.
#' @param horizon For rolling CV, forecast horizon.
#' @param refit_derivative Logical; whether to recompute derivatives for each fold (currently unused).
#' @param diff_method Differentiation method if refitting (currently unused).
#' @param verbose Print progress.
#'
#' @return Object of class "cv_result" containing:
#'   \item{rmse}{Root mean squared error per fold}
#'   \item{mae}{Mean absolute error per fold}
#'   \item{r_squared}{R-squared per fold}
#'   \item{mean_rmse}{Average RMSE across folds}
#'   \item{sd_rmse}{Standard deviation of RMSE}
#'   \item{predictions}{List of predicted vs actual per fold}
#'   \item{fold_indices}{Indices used for each fold}
#'
#' @examples
#' \donttest{
#' # Toy example using lm
#' data <- data.frame(
#'   time = 1:50,
#'   y = seq(1, 10, length.out = 50) + stats::rnorm(50, sd = 0.1)
#' )
#' # Simple linear model as a proxy for a discovered equation
#' model <- stats::lm(y ~ time, data = data)
#'
#' # Run cross-validation
#' cv_res <- cross_validate(
#'   equation = model,
#'   data = data,
#'   response = "y",
#'   k = 3,
#'   method = "random"
#' )
#' print(cv_res)
#' }
#' @export
cross_validate <- function(equation, data, response = NULL,
                           derivative_col = NULL,
                           k = 5,
                           method = c("block", "random", "rolling"),
                           block_size = NULL,
                           horizon = 1,
                           refit_derivative = FALSE,
                           diff_method = "tvr",
                           verbose = TRUE) {
  
  # Compatibility: derivative_col is alias for response
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  if (is.null(response)) {
    stop("Must specify 'response' or 'derivative_col'", call. = FALSE)
  }
  
  method <- match.arg(method)
  n <- nrow(data)
  
  # Generate fold indices
  fold_indices <- switch(method,
                         "random" = create_random_folds(n, k),
                         "block" = create_block_folds(n, k, block_size),
                         "rolling" = create_rolling_folds(n, k, horizon)
  )
  
  # Extract formula from equation
  if (inherits(equation, "symbolic_equation")) {
    eq_expr <- equation$expression
    eq_type <- "symbolic"
  } else if (inherits(equation, "nls")) {
    form <- stats::formula(equation)
    eq_type <- "nls"
  } else if (inherits(equation, "lm")) {
    form <- stats::formula(equation)
    eq_type <- "lm"
  } else {
    stop("Unknown equation type. Expected symbolic_equation, nls, or lm object.")
  }
  
  # Storage for results
  results <- list(
    rmse = numeric(k),
    mae = numeric(k),
    r_squared = numeric(k),
    predictions = vector("list", k),
    fold_indices = fold_indices
  )
  
  for (i in seq_len(k)) {
    if (verbose) message("Fold ", i, "/", k)
    
    test_idx <- fold_indices[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)
    
    train_data <- data[train_idx, , drop = FALSE]
    test_data <- data[test_idx, , drop = FALSE]
    
    # Refit on training data
    fitted_model <- tryCatch({
      if (eq_type == "lm") {
        stats::lm(form, data = train_data)
      } else if (eq_type == "nls") {
        if (!requireNamespace("minpack.lm", quietly = TRUE)) {
          stop("Package 'minpack.lm' is required for NLS fitting.")
        }
        # For NLS, need starting values
        start_vals <- as.list(stats::coef(equation))
        minpack.lm::nlsLM(form, data = train_data, start = start_vals,
                          control = minpack.lm::nls.lm.control(maxiter = 200))
      } else {
        # Symbolic equation - use fit_specified_equation
        fit_specified_equation(
          eq_expr,
          data = train_data,
          response = response,
          start = as.list(stats::coef(equation))
        )
      }
    }, error = function(e) {
      warning("Fold ", i, " fitting failed: ", e$message)
      NULL
    })
    
    if (is.null(fitted_model)) {
      results$rmse[i] <- NA
      results$mae[i] <- NA
      results$r_squared[i] <- NA
      next
    }
    
    # Predict on test data
    pred <- stats::predict(fitted_model, newdata = test_data)
    actual <- test_data[[response]]
    
    # Calculate metrics
    residuals <- actual - pred
    results$rmse[i] <- sqrt(mean(residuals^2, na.rm = TRUE))
    results$mae[i] <- mean(abs(residuals), na.rm = TRUE)
    
    ss_res <- sum(residuals^2, na.rm = TRUE)
    ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
    results$r_squared[i] <- 1 - ss_res / ss_tot
    
    results$predictions[[i]] <- data.frame(
      index = test_idx,
      actual = actual,
      predicted = pred,
      residual = residuals
    )
  }
  
  # Summary statistics
  results$mean_rmse <- mean(results$rmse, na.rm = TRUE)
  results$sd_rmse <- stats::sd(results$rmse, na.rm = TRUE)
  results$mean_mae <- mean(results$mae, na.rm = TRUE)
  results$mean_r_squared <- mean(results$r_squared, na.rm = TRUE)
  results$method <- method
  results$k <- k
  
  class(results) <- "cv_result"
  return(results)
}

#' Create Random Folds
#' @keywords internal
create_random_folds <- function(n, k) {
  folds <- sample(rep(1:k, length.out = n))
  lapply(1:k, function(i) which(folds == i))
}

#' Create Block Folds (for time series)
#' @keywords internal
create_block_folds <- function(n, k, block_size = NULL) {
  if (is.null(block_size)) {
    block_size <- floor(n / k)
  }
  
  # Create contiguous blocks
  starts <- seq(1, n, by = block_size)
  if (length(starts) > k) starts <- starts[1:k]
  
  lapply(seq_along(starts), function(i) {
    start <- starts[i]
    end <- min(start + block_size - 1, n)
    start:end
  })
}

#' Create Rolling Folds (walk-forward validation)
#' @keywords internal
create_rolling_folds <- function(n, k, horizon = 1) {
  # Expanding window with fixed test horizon
  min_train <- floor(n / (k + 1))
  
  lapply(1:k, function(i) {
    train_end <- min_train + (i - 1) * horizon
    test_start <- train_end + 1
    test_end <- min(test_start + horizon - 1, n)
    test_start:test_end
  })
}

#' Print CV Results
#'
#' @param x Object of class cv_result.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
print.cv_result <- function(x, ...) {
  cat("Cross-Validation Results\n")
  cat("========================\n")
  cat("Method:", x$method, "with", x$k, "folds\n\n")
  cat("RMSE:     ", sprintf("%.4f (SD: %.4f)", x$mean_rmse, x$sd_rmse), "\n")
  cat("MAE:      ", sprintf("%.4f", x$mean_mae), "\n")
  cat("R-squared:", sprintf("%.4f", x$mean_r_squared), "\n\n")
  cat("Per-fold RMSE:", paste(sprintf("%.4f", x$rmse), collapse = ", "), "\n")
  invisible(x)
}

#' Plot CV Results
#'
#' @param x Object of class cv_result.
#' @param type Type of plot: "predictions", "folds", or "both".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object or a list of ggplot objects.
#' @export
plot.cv_result <- function(x, type = c("predictions", "folds", "both"), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  type <- match.arg(type)
  
  plots <- list()
  
  if (type %in% c("predictions", "both")) {
    # Combine all predictions
    all_pred <- do.call(rbind, x$predictions)
    
    plots$pred <- ggplot2::ggplot(all_pred, ggplot2::aes(x = actual, y = predicted)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      ggplot2::labs(
        title = "Cross-Validation: Predicted vs Actual",
        subtitle = sprintf("Mean RMSE: %.4f", x$mean_rmse),
        x = "Actual dZ/dt",
        y = "Predicted dZ/dt"
      ) +
      ed_theme()
  }
  
  if (type %in% c("folds", "both")) {
    fold_df <- data.frame(
      fold = factor(1:x$k),
      rmse = x$rmse,
      mae = x$mae,
      r_squared = x$r_squared
    )
    
    plots$folds <- ggplot2::ggplot(fold_df, ggplot2::aes(x = fold, y = rmse)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::geom_hline(yintercept = x$mean_rmse, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "RMSE by Fold",
        x = "Fold",
        y = "RMSE"
      ) +
      ed_theme()
  }
  
  if (length(plots) == 1) {
    print(plots[[1]])
    return(invisible(plots[[1]]))
  } else if (length(plots) > 1) {
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(grobs = plots, ncol = 2)
    } else {
      warning("Package 'gridExtra' needed to plot both panels simultaneously.")
      print(plots$pred)
      print(plots$folds)
    }
  }
  
  invisible(plots)
}


# =============================================================================
# TRAJECTORY SIMULATION
# =============================================================================

#' Simulate Trajectory from SDE
#'
#' Simulates trajectories using the discovered SDE to assess whether the model
#' can reproduce observed dynamics.
#'
#' @param sde SDE object from \code{construct_sde} or \code{estimate_sde_iterative}.
#' @param initial_conditions Named vector of initial values for all variables.
#' @param times Numeric vector of time points.
#' @param n_sims Number of Monte Carlo simulations (for stochastic models).
#' @param method Integration method: "euler", "milstein", "rk4" (deterministic only).
#' @param exogenous_data Data frame with exogenous variable trajectories (if any).
#' @param seed Random seed for reproducibility.
#'
#' @return Object of class "trajectory_simulation" containing:
#'   \item{trajectories}{Array of simulated trajectories (time x variable x simulation)}
#'   \item{times}{Time points}
#'   \item{summary}{Summary statistics (mean, quantiles) at each time}
#'
#' @examples
#' \donttest{
#' # Toy example: dX = 0.5 * X
#' # Mock SDE object structure
#' sde <- list(
#'   drift = list(expression = "0.5 * X"),
#'   diffusion = list(expression = "0.1"), # Add noise
#'   variable = "X"
#' )
#' class(sde) <- "sde_model"
#'
#' # Simulation
#' sim <- simulate_trajectory(
#'   sde = sde,
#'   initial_conditions = c(X = 1),
#'   times = seq(0, 1, by = 0.1),
#'   n_sims = 10,
#'   seed = 123
#' )
#' print(sim$summary$mean)
#' }
#' @export
simulate_trajectory <- function(sde, initial_conditions, times,
                                n_sims = 100,
                                method = c("euler", "milstein", "rk4"),
                                exogenous_data = NULL,
                                seed = NULL) {
  
  method <- match.arg(method)
  
  # Set seed if provided (without modifying .GlobalEnv per CRAN policy)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Parse SDE components
  if (inherits(sde, "sde_model")) {
    drift_expr <- if (!is.null(sde$drift$expression)) {
      sde$drift$expression
    } else {
      sde$drift$string
    }
    diffusion_expr <- if (!is.null(sde$diffusion)) {
      if (!is.null(sde$diffusion$expression)) {
        sde$diffusion$expression
      } else {
        "0"  # No diffusion
      }
    } else {
      "0"
    }
    drift_coefs <- stats::coef(sde$drift)
    diffusion_coefs <- if (!is.null(sde$diffusion) && inherits(sde$diffusion$fit, "lm")) {
      stats::coef(sde$diffusion$fit)
    } else {
      NULL
    }
    var_name <- sde$variable
    if (is.null(var_name)) {
      var_name <- names(initial_conditions)[1]
    }
  } else {
    stop("sde must be an object from construct_sde or estimate_sde_iterative")
  }
  
  n_times <- length(times)
  dt <- diff(times)
  
  # Storage for trajectories
  trajectories <- array(NA, dim = c(n_times, length(initial_conditions), n_sims),
                        dimnames = list(NULL, names(initial_conditions), NULL))
  
  # Initialize
  for (s in 1:n_sims) {
    trajectories[1, , s] <- initial_conditions
  }
  
  # Create evaluation environment
  eval_env <- new.env()
  if (!is.null(drift_coefs)) {
    for (nm in names(drift_coefs)) {
      eval_env[[nm]] <- drift_coefs[nm]
    }
  }
  if (!is.null(diffusion_coefs)) {
    for (nm in names(diffusion_coefs)) {
      eval_env[[nm]] <- diffusion_coefs[nm]
    }
  }
  
  # Simulation loop
  for (s in 1:n_sims) {
    for (t in 1:(n_times - 1)) {
      current_dt <- dt[t]
      
      # Set current state in environment
      for (nm in names(initial_conditions)) {
        eval_env[[nm]] <- trajectories[t, nm, s]
      }
      
      # Add exogenous variables if provided
      if (!is.null(exogenous_data)) {
        for (nm in setdiff(names(exogenous_data), "time")) {
          # Interpolate exogenous data to current time
          eval_env[[nm]] <- stats::approx(exogenous_data$time, exogenous_data[[nm]], 
                                          times[t], rule = 2)$y
        }
      }
      
      # Evaluate drift
      drift_val <- tryCatch(
        eval(parse(text = drift_expr), envir = eval_env),
        error = function(e) NA
      )
      
      # Evaluate diffusion
      diffusion_val <- tryCatch(
        eval(parse(text = diffusion_expr), envir = eval_env),
        error = function(e) 0
      )
      
      if (is.na(drift_val)) {
        trajectories[(t+1):n_times, , s] <- NA
        break
      }
      
      # Generate noise
      dW <- stats::rnorm(1, mean = 0, sd = sqrt(current_dt))
      
      # Integration step
      if (method == "euler") {
        # Euler-Maruyama
        new_val <- trajectories[t, var_name, s] + 
          drift_val * current_dt + 
          diffusion_val * dW
        
      } else if (method == "milstein") {
        # Milstein scheme (for scalar diffusion)
        # Requires derivative of diffusion - approximate numerically
        eps <- 1e-6
        eval_env[[var_name]] <- trajectories[t, var_name, s] + eps
        diffusion_plus <- tryCatch(
          eval(parse(text = diffusion_expr), envir = eval_env),
          error = function(e) diffusion_val
        )
        diffusion_deriv <- (diffusion_plus - diffusion_val) / eps
        
        new_val <- trajectories[t, var_name, s] + 
          drift_val * current_dt + 
          diffusion_val * dW +
          0.5 * diffusion_val * diffusion_deriv * (dW^2 - current_dt)
        
      } else if (method == "rk4") {
        # RK4 (deterministic part only)
        eval_env[[var_name]] <- trajectories[t, var_name, s]
        k1 <- eval(parse(text = drift_expr), envir = eval_env)
        
        eval_env[[var_name]] <- trajectories[t, var_name, s] + 0.5 * current_dt * k1
        k2 <- eval(parse(text = drift_expr), envir = eval_env)
        
        eval_env[[var_name]] <- trajectories[t, var_name, s] + 0.5 * current_dt * k2
        k3 <- eval(parse(text = drift_expr), envir = eval_env)
        
        eval_env[[var_name]] <- trajectories[t, var_name, s] + current_dt * k3
        k4 <- eval(parse(text = drift_expr), envir = eval_env)
        
        new_val <- trajectories[t, var_name, s] + 
          (current_dt / 6) * (k1 + 2*k2 + 2*k3 + k4) +
          diffusion_val * dW
      }
      
      trajectories[t + 1, var_name, s] <- new_val
    }
  }
  
  # Compute summary statistics
  summary_stats <- list(
    mean = apply(trajectories[, var_name, , drop = FALSE], 1, mean, na.rm = TRUE),
    sd = apply(trajectories[, var_name, , drop = FALSE], 1, stats::sd, na.rm = TRUE),
    q05 = apply(trajectories[, var_name, , drop = FALSE], 1, stats::quantile, 0.05, na.rm = TRUE),
    q25 = apply(trajectories[, var_name, , drop = FALSE], 1, stats::quantile, 0.25, na.rm = TRUE),
    q50 = apply(trajectories[, var_name, , drop = FALSE], 1, stats::quantile, 0.50, na.rm = TRUE),
    q75 = apply(trajectories[, var_name, , drop = FALSE], 1, stats::quantile, 0.75, na.rm = TRUE),
    q95 = apply(trajectories[, var_name, , drop = FALSE], 1, stats::quantile, 0.95, na.rm = TRUE)
  )
  
  result <- list(
    trajectories = trajectories,
    times = times,
    summary = summary_stats,
    n_sims = n_sims,
    method = method,
    variable = var_name,
    initial_conditions = initial_conditions
  )
  
  class(result) <- "trajectory_simulation"
  return(result)
}

#' Plot Simulated Trajectories
#'
#' @param x Object of class trajectory_simulation.
#' @param observed_data Optional observed data to overlay.
#' @param show_trajectories Show individual trajectories?
#' @param n_show Number of trajectories to show.
#' @param alpha_traj Transparency for trajectories.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#' @export
plot.trajectory_simulation <- function(x, observed_data = NULL, 
                                       show_trajectories = TRUE,
                                       n_show = 20,
                                       alpha_traj = 0.2, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  var_name <- x$variable
  
  # Create base data frame
  df_summary <- data.frame(
    time = x$times,
    mean = x$summary$mean,
    q05 = x$summary$q05,
    q25 = x$summary$q25,
    q50 = x$summary$q50,
    q75 = x$summary$q75,
    q95 = x$summary$q95
  )
  
  p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = time))
  
  # Add individual trajectories
  if (show_trajectories && x$n_sims > 0) {
    n_to_show <- min(n_show, x$n_sims)
    for (s in 1:n_to_show) {
      traj_df <- data.frame(
        time = x$times,
        value = x$trajectories[, var_name, s]
      )
      p <- p + ggplot2::geom_line(
        data = traj_df,
        ggplot2::aes(y = value),
        alpha = alpha_traj,
        color = "gray50"
      )
    }
  }
  
  # Add confidence ribbon
  p <- p + 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95), 
                         fill = "steelblue", alpha = 0.2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), 
                         fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = q50), color = "steelblue", linewidth = 1)
  
  # Add observed data if provided
  if (!is.null(observed_data)) {
    if ("time" %in% names(observed_data) && var_name %in% names(observed_data)) {
      p <- p + ggplot2::geom_point(
        data = observed_data,
        ggplot2::aes(x = time, y = .data[[var_name]]),
        color = "red", size = 2
      )
    }
  }
  
  p <- p + 
    ggplot2::labs(
      title = paste("Simulated Trajectories:", var_name),
      subtitle = sprintf("N = %d simulations, Method: %s", x$n_sims, x$method),
      x = "Time",
      y = var_name
    ) +
    ed_theme()
  
  print(p)
  invisible(p)
}

#' Compare Simulated and Observed Trajectories
#'
#' Computes metrics comparing simulated trajectories to observed data.
#'
#' @param simulation Trajectory simulation object.
#' @param observed_data Data frame with observed values.
#' @param time_col Name of time column.
#' @param var_col Name of variable column to compare.
#'
#' @return Data frame with comparison metrics.
#' @export
compare_trajectories <- function(simulation, observed_data, 
                                 time_col = "time", var_col = NULL) {
  
  if (is.null(var_col)) var_col <- simulation$variable
  
  # Interpolate simulated summary to observed times
  obs_times <- observed_data[[time_col]]
  obs_values <- observed_data[[var_col]]
  
  sim_mean <- stats::approx(simulation$times, simulation$summary$mean, obs_times)$y
  sim_median <- stats::approx(simulation$times, simulation$summary$q50, obs_times)$y
  sim_q05 <- stats::approx(simulation$times, simulation$summary$q05, obs_times)$y
  sim_q95 <- stats::approx(simulation$times, simulation$summary$q95, obs_times)$y
  
  # Coverage: proportion of observed values within 90% CI
  coverage_90 <- mean(obs_values >= sim_q05 & obs_values <= sim_q95, na.rm = TRUE)
  
  # RMSE of mean trajectory
  rmse_mean <- sqrt(mean((obs_values - sim_mean)^2, na.rm = TRUE))
  rmse_median <- sqrt(mean((obs_values - sim_median)^2, na.rm = TRUE))
  
  # Correlation
  cor_mean <- stats::cor(obs_values, sim_mean, use = "complete.obs")
  
  # Mean absolute deviation
  mad_mean <- mean(abs(obs_values - sim_mean), na.rm = TRUE)
  
  metrics <- data.frame(
    metric = c("RMSE (mean)", "RMSE (median)", "MAD (mean)", 
               "Correlation", "Coverage (90% CI)"),
    value = c(rmse_mean, rmse_median, mad_mean, cor_mean, coverage_90)
  )
  
  return(metrics)
}


# =============================================================================
# QUALITATIVE VALIDATION
# =============================================================================

#' Analyze Fixed Points
#'
#' Finds and characterizes fixed points of the discovered equation.
#'
#' @param equation Fitted equation object.
#' @param variable Name of the main variable.
#' @param range Numeric vector of length 2 specifying search range.
#' @param n_grid Number of grid points for initial search.
#' @param exogenous_values Named list of fixed values for exogenous variables.
#'
#' @return Data frame of fixed points with stability classification.
#' @examples
#' \donttest{
#' # Toy example: dZ = 2*Z - Z^2 (Logistic growth)
#' data <- data.frame(Z = seq(0.1, 3, length.out=50))
#' data$dZ <- 2*data$Z - data$Z^2
#' model <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = data)
#'
#' # Analyze (note: linear models on dZ aren't direct ODEs, but this demonstrates structure)
#' # For correct usage, 'equation' should be from fit_specified_equation
#' fp <- analyze_fixed_points(model, variable = "Z", range = c(0, 3))
#' }
#' @export
analyze_fixed_points <- function(equation, variable, 
                                 range = c(-10, 10),
                                 n_grid = 100,
                                 exogenous_values = list()) {
  
  # Extract expression
  if (inherits(equation, "symbolic_equation")) {
    expr_str <- if (!is.null(equation$expression)) {
      equation$expression
    } else {
      equation$string
    }
    coefs <- stats::coef(equation)
  } else if (inherits(equation, "nls") || inherits(equation, "lm")) {
    expr_str <- deparse(stats::formula(equation)[[3]])
    coefs <- stats::coef(equation)
  } else {
    stop("Unknown equation type")
  }
  
  # Create evaluation function
  eval_env <- new.env()
  if (!is.null(coefs)) {
    for (nm in names(coefs)) {
      eval_env[[nm]] <- coefs[nm]
    }
  }
  for (nm in names(exogenous_values)) {
    eval_env[[nm]] <- exogenous_values[[nm]]
  }
  
  f <- function(z) {
    eval_env[[variable]] <- z
    tryCatch(
      eval(parse(text = expr_str), envir = eval_env),
      error = function(e) NA
    )
  }
  
  # Grid search for sign changes
  grid <- seq(range[1], range[2], length.out = n_grid)
  f_vals <- sapply(grid, f)
  
  # Find sign changes
  valid_idx <- which(!is.na(f_vals))
  if (length(valid_idx) < 2) {
    message("Could not evaluate function on grid")
    return(data.frame(
      fixed_point = numeric(0),
      stability = character(0),
      eigenvalue = numeric(0)
    ))
  }
  
  sign_changes <- which(diff(sign(f_vals[valid_idx])) != 0)
  
  fixed_points <- numeric(0)
  
  for (i in sign_changes) {
    idx1 <- valid_idx[i]
    idx2 <- valid_idx[i + 1]
    # Refine with uniroot
    fp <- tryCatch({
      stats::uniroot(f, c(grid[idx1], grid[idx2]))$root
    }, error = function(e) NA)
    
    if (!is.na(fp)) {
      fixed_points <- c(fixed_points, fp)
    }
  }
  
  if (length(fixed_points) == 0) {
    message("No fixed points found in specified range")
    return(data.frame(
      fixed_point = numeric(0),
      stability = character(0),
      eigenvalue = numeric(0)
    ))
  }
  
  # Classify stability (compute derivative at fixed point)
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
      
      if (df_dz < -eps) {
        stability[i] <- "stable"
      } else if (df_dz > eps) {
        stability[i] <- "unstable"
      } else {
        stability[i] <- "marginal"
      }
    }
  }
  
  result <- data.frame(
    fixed_point = fixed_points,
    stability = stability,
    eigenvalue = eigenvalues
  )
  
  return(result)
}

#' Analyze Bifurcations
#'
#' Examines how fixed points change as a parameter varies.
#'
#' @param equation Fitted equation object.
#' @param variable Name of the main variable.
#' @param parameter Name of the parameter to vary.
#' @param param_range Range for parameter values.
#' @param n_param Number of parameter values to test.
#' @param z_range Range for searching fixed points.
#' @param exogenous_values Fixed values for other variables.
#'
#' @return Object of class "bifurcation_analysis".
#' @export
analyze_bifurcations <- function(equation, variable, parameter,
                                 param_range = c(-5, 5),
                                 n_param = 50,
                                 z_range = c(-10, 10),
                                 exogenous_values = list()) {
  
  param_vals <- seq(param_range[1], param_range[2], length.out = n_param)
  
  all_fps <- list()
  
  for (i in seq_along(param_vals)) {
    # Set parameter value
    exog <- exogenous_values
    exog[[parameter]] <- param_vals[i]
    
    fps <- tryCatch({
      analyze_fixed_points(equation, variable, range = z_range,
                           exogenous_values = exog)
    }, error = function(e) {
      data.frame(fixed_point = numeric(0), stability = character(0), eigenvalue = numeric(0))
    })
    
    if (nrow(fps) > 0) {
      fps$parameter_value <- param_vals[i]
      fps$parameter_name <- parameter
      all_fps[[i]] <- fps
    }
  }
  
  bifurc_data <- do.call(rbind, all_fps)
  
  result <- list(
    data = bifurc_data,
    parameter = parameter,
    variable = variable,
    param_range = param_range
  )
  
  class(result) <- "bifurcation_analysis"
  return(result)
}

#' Plot Bifurcation Diagram
#'
#' @param x Object of class bifurcation_analysis.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#' @export
plot.bifurcation_analysis <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  if (is.null(x$data) || nrow(x$data) == 0) {
    message("No bifurcation data to plot")
    return(invisible(NULL))
  }
  
  p <- ggplot2::ggplot(x$data, 
                       ggplot2::aes(x = parameter_value, y = fixed_point, 
                                    color = stability)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(
      values = c("stable" = "blue", "unstable" = "red", "marginal" = "gray50", "unknown" = "gray80")
    ) +
    ggplot2::labs(
      title = "Bifurcation Diagram",
      x = x$parameter,
      y = paste(x$variable, "fixed point"),
      color = "Stability"
    ) +
    ed_theme()
  
  print(p)
  invisible(p)
}

#' Check Qualitative Behavior
#'
#' Comprehensive check of whether the discovered equation exhibits
#' expected qualitative features.
#'
#' @param equation Fitted equation object.
#' @param data Original data.
#' @param variable Main variable name.
#' @param expected_features List of expected qualitative features:
#'   \itemize{
#'     \item n_fixed_points: Expected number of fixed points
#'     \item stability_pattern: e.g., c("stable", "unstable", "stable")
#'     \item monotonicity: Expected sign of derivative ("positive", "negative", "none")
#'     \item bounded: Whether dynamics should be bounded
#'   }
#'
#' @return Object of class "qualitative_check".
#' @export
check_qualitative_behavior <- function(equation, data, variable,
                                       expected_features = list()) {
  
  results <- list(
    checks = list(),
    passed = logical(0),
    messages = character(0)
  )
  
  # Determine data range for the variable
  if (variable %in% names(data)) {
    var_range <- range(data[[variable]], na.rm = TRUE)
    var_range <- var_range + c(-1, 1) * 0.2 * diff(var_range)  # Expand 20%
  } else {
    var_range <- c(-10, 10)
  }
  
  # Check 1: Number of fixed points
  fps <- analyze_fixed_points(equation, variable, range = var_range)
  results$checks$fixed_points <- fps
  
  if (!is.null(expected_features$n_fixed_points)) {
    passed <- nrow(fps) == expected_features$n_fixed_points
    results$passed <- c(results$passed, fp_count = passed)
    results$messages <- c(results$messages,
                          sprintf("Fixed points: found %d, expected %d - %s",
                                  nrow(fps), expected_features$n_fixed_points,
                                  if(passed) "PASS" else "FAIL"))
  }
  
  # Check 2: Stability pattern
  if (!is.null(expected_features$stability_pattern) && nrow(fps) > 0) {
    actual_pattern <- fps$stability[order(fps$fixed_point)]
    expected_pattern <- expected_features$stability_pattern
    
    passed <- length(actual_pattern) == length(expected_pattern) &&
      all(actual_pattern == expected_pattern)
    
    results$passed <- c(results$passed, stability = passed)
    results$messages <- c(results$messages,
                          sprintf("Stability pattern: %s vs expected %s - %s",
                                  paste(actual_pattern, collapse = ", "),
                                  paste(expected_pattern, collapse = ", "),
                                  if(passed) "PASS" else "FAIL"))
  }
  
  results$n_passed <- sum(results$passed)
  results$n_total <- length(results$passed)
  
  class(results) <- "qualitative_check"
  return(results)
}

#' Print Qualitative Check Results
#'
#' @param x Object of class qualitative_check.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
print.qualitative_check <- function(x, ...) {
  cat("Qualitative Behavior Check\n")
  cat("==========================\n\n")
  
  for (msg in x$messages) {
    cat(msg, "\n")
  }
  
  cat("\n")
  cat(sprintf("Summary: %d/%d checks passed\n", x$n_passed, x$n_total))
  
  if (x$n_passed == x$n_total && x$n_total > 0) {
    cat("All qualitative checks PASSED\n")
  }
  
  invisible(x)
}


# =============================================================================
# COMPREHENSIVE VALIDATION
# =============================================================================

#' Comprehensive Model Validation
#'
#' Runs a battery of validation tests on the discovered equation.
#'
#' @param equation Fitted equation object.
#' @param sde SDE object (optional, for trajectory validation).
#' @param data Original data frame.
#' @param response Name of response column.
#' @param derivative_col Alias for response (for compatibility).
#' @param variable Main variable name.
#' @param time_col Time column name.
#' @param cv_folds Number of CV folds.
#' @param n_sims Number of trajectory simulations.
#' @param expected_features List of expected qualitative features.
#' @param verbose Print progress.
#'
#' @return Object of class "validation_result".
#' @export
validate_model <- function(equation, sde = NULL, data, response = NULL,
                           derivative_col = NULL, variable,
                           time_col = "time",
                           cv_folds = 5,
                           n_sims = 50,
                           expected_features = list(),
                           verbose = TRUE) {
  
  # Compatibility
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  if (is.null(response)) {
    # Try to auto-detect
    d_cols <- grep("^d_|^d[A-Z]", names(data), value = TRUE)
    if (length(d_cols) == 1) {
      response <- d_cols[1]
    } else {
      stop("Must specify 'response' or 'derivative_col'", call. = FALSE)
    }
  }
  
  results <- list()
  
  # 1. Cross-validation
  if (verbose) message("Running cross-validation...")
  results$cv <- tryCatch({
    cross_validate(equation, data, response = response, k = cv_folds, 
                   method = "block", verbose = FALSE)
  }, error = function(e) {
    warning("Cross-validation failed: ", e$message)
    NULL
  })
  
  # 2. In-sample fit statistics
  if (verbose) message("Computing fit statistics...")
  results$fit_stats <- tryCatch({
    pred <- stats::predict(equation, newdata = data)
    actual <- data[[response]]
    residuals <- actual - pred
    
    list(
      r_squared = 1 - sum(residuals^2, na.rm = TRUE) / sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE),
      adj_r_squared = NA,
      rmse = sqrt(mean(residuals^2, na.rm = TRUE)),
      mae = mean(abs(residuals), na.rm = TRUE),
      aic = tryCatch(stats::AIC(equation), error = function(e) NA),
      bic = tryCatch(stats::BIC(equation), error = function(e) NA)
    )
  }, error = function(e) {
    list(r_squared = NA, rmse = NA, mae = NA, aic = NA, bic = NA)
  })
  
  # 3. Residual diagnostics
  if (verbose) message("Running residual diagnostics...")
  results$residuals <- tryCatch({
    if (inherits(equation, "symbolic_equation") && !is.null(equation$fit)) {
      resid <- stats::residuals(equation$fit)
    } else if (inherits(equation, "lm") || inherits(equation, "nls")) {
      resid <- stats::residuals(equation)
    } else {
      pred <- stats::predict(equation, newdata = data)
      resid <- data[[response]] - pred
    }
    residual_diagnostics(resid, data, plot = FALSE)
  }, error = function(e) NULL)
  
  # 4. Trajectory simulation (if SDE provided)
  if (!is.null(sde) && time_col %in% names(data)) {
    if (verbose) message("Simulating trajectories...")
    results$trajectory <- tryCatch({
      times <- data[[time_col]]
      init <- setNames(data[[variable]][1], variable)
      
      sim <- simulate_trajectory(sde, init, times, n_sims = n_sims)
      
      # Compare with observed
      metrics <- compare_trajectories(sim, data, time_col, variable)
      
      list(simulation = sim, comparison = metrics)
    }, error = function(e) {
      warning("Trajectory simulation failed: ", e$message)
      NULL
    })
  }
  
  # 5. Qualitative checks
  if (verbose) message("Checking qualitative behavior...")
  results$qualitative <- tryCatch({
    check_qualitative_behavior(equation, data, variable, expected_features)
  }, error = function(e) NULL)
  
  # Overall assessment
  results$summary <- list(
    cv_rmse = if (!is.null(results$cv)) results$cv$mean_rmse else NA,
    in_sample_r2 = results$fit_stats$r_squared,
    residual_tests_passed = if (!is.null(results$residuals) && !is.null(results$residuals$tests)) {
      sum(results$residuals$tests$p_value > 0.05, na.rm = TRUE)
    } else {
      NA
    },
    qualitative_passed = if (!is.null(results$qualitative)) {
      results$qualitative$n_passed
    } else {
      NA
    }
  )
  
  class(results) <- "validation_result"
  
  if (verbose) {
    message("\nValidation complete.")
    print(results)
  }
  
  return(results)
}

#' Print Validation Results
#'
#' @param x Object of class validation_result.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
print.validation_result <- function(x, ...) {
  cat("\n")
  cat("Model Validation Summary\n")
  cat("========================\n\n")
  
  cat("Predictive Performance:\n")
  cat(sprintf("  Cross-validation RMSE: %.4f\n", x$summary$cv_rmse))
  cat(sprintf("  In-sample R-squared:   %.4f\n", x$summary$in_sample_r2))
  
  if (!is.null(x$fit_stats)) {
    if (!is.na(x$fit_stats$aic)) cat(sprintf("  AIC: %.2f\n", x$fit_stats$aic))
    if (!is.na(x$fit_stats$bic)) cat(sprintf("  BIC: %.2f\n", x$fit_stats$bic))
  }
  
  cat("\nResidual Diagnostics:\n")
  if (!is.null(x$residuals) && !is.null(x$residuals$tests)) {
    n_tests <- nrow(x$residuals$tests)
    n_passed <- sum(x$residuals$tests$p_value > 0.05, na.rm = TRUE)
    cat(sprintf("  Tests passed (p > 0.05): %d/%d\n", n_passed, n_tests))
  } else {
    cat("  Not available\n")
  }
  
  if (!is.null(x$trajectory)) {
    cat("\nTrajectory Comparison:\n")
    print(x$trajectory$comparison)
  }
  
  if (!is.null(x$qualitative)) {
    cat("\nQualitative Checks:\n")
    cat(sprintf("  Passed: %d/%d\n", 
                x$qualitative$n_passed, 
                x$qualitative$n_total))
  }
  
  invisible(x)
}

#' Plot Validation Results
#'
#' @param x Object of class validation_result.
#' @param ... Additional arguments (ignored).
#'
#' @return A list of ggplot objects (invisible).
#' @export
plot.validation_result <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  plots <- list()
  
  # CV results
  if (!is.null(x$cv)) {
    plots$cv <- tryCatch(plot(x$cv, type = "predictions"), error = function(e) NULL)
  }
  
  # Residual diagnostics
  if (!is.null(x$residuals)) {
    # Note: plot_residual_diagnostics_panel must handle its own imports (hist, curve, etc.)
    tryCatch(plot_residual_diagnostics_panel(x$residuals), error = function(e) NULL)
  }
  
  # Trajectory comparison
  if (!is.null(x$trajectory) && !is.null(x$trajectory$simulation)) {
    plots$traj <- tryCatch(plot(x$trajectory$simulation), error = function(e) NULL)
  }
  
  invisible(plots)
}


# =============================================================================
# SENSITIVITY ANALYSIS
# =============================================================================

#' Parameter Sensitivity Analysis
#'
#' Examines how sensitive the model predictions are to parameter perturbations.
#'
#' @param equation Fitted equation object.
#' @param data Data for evaluation.
#' @param response Name of response column.
#' @param derivative_col Alias for response.
#' @param perturbation_pct Percentage perturbation (default 10%).
#' @param n_bootstrap Number of bootstrap samples for uncertainty.
#'
#' @return Data frame with sensitivity metrics for each parameter.
#' @export
sensitivity_analysis <- function(equation, data, response = NULL,
                                 derivative_col = NULL,
                                 perturbation_pct = 10,
                                 n_bootstrap = 100) {
  
  # Compatibility
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  
  coefs <- stats::coef(equation)
  if (is.null(coefs) || length(coefs) == 0) {
    warning("No coefficients found in equation")
    return(data.frame(parameter = character(0), estimate = numeric(0), 
                      sensitivity = numeric(0), elasticity = numeric(0)))
  }
  
  n_coef <- length(coefs)
  baseline_pred <- stats::predict(equation, newdata = data)
  baseline_rmse <- sqrt(mean((data[[response]] - baseline_pred)^2, na.rm = TRUE))
  
  results <- data.frame(
    parameter = names(coefs),
    estimate = coefs,
    sensitivity = NA,
    elasticity = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n_coef)) {
    # Perturb parameter up and down
    delta <- abs(coefs[i]) * perturbation_pct / 100
    if (delta < 1e-10) delta <- perturbation_pct / 100
    
    # Create modified coefficient vectors
    coefs_up <- coefs
    coefs_up[i] <- coefs[i] + delta
    
    coefs_down <- coefs
    coefs_down[i] <- coefs[i] - delta
    
    # Compute predictions with perturbed parameters
    pred_up <- eval_with_coefs(equation, data, coefs_up)
    pred_down <- eval_with_coefs(equation, data, coefs_down)
    
    if (!is.null(pred_up) && !is.null(pred_down) && 
        !all(is.na(pred_up)) && !all(is.na(pred_down))) {
      # Sensitivity: change in prediction per unit change in parameter
      sensitivity <- mean(abs(pred_up - pred_down) / (2 * delta), na.rm = TRUE)
      
      # Elasticity: % change in prediction per % change in parameter
      elasticity <- mean(abs((pred_up - pred_down) / baseline_pred) / 
                           (2 * perturbation_pct / 100), na.rm = TRUE)
      
      results$sensitivity[i] <- sensitivity
      results$elasticity[i] <- elasticity
    }
  }
  
  # Rank parameters by sensitivity
  results <- results[order(-results$sensitivity, na.last = TRUE), ]
  rownames(results) <- NULL
  
  return(results)
}

#' Evaluate equation with modified coefficients
#' @keywords internal
eval_with_coefs <- function(equation, data, new_coefs) {
  tryCatch({
    if (inherits(equation, "symbolic_equation")) {
      expr_str <- if (!is.null(equation$expression)) {
        equation$expression
      } else {
        equation$string
      }
      eval_env <- new.env()
      for (nm in names(new_coefs)) {
        eval_env[[nm]] <- new_coefs[nm]
      }
      
      sapply(1:nrow(data), function(i) {
        for (nm in names(data)) {
          eval_env[[nm]] <- data[[nm]][i]
        }
        eval(parse(text = expr_str), envir = eval_env)
      })
    } else {
      # For nls/lm, harder to modify - would need to refit
      NULL
    }
  }, error = function(e) NULL)
}


#' Bootstrap Confidence Intervals for Parameters
#'
#' Computes bootstrap confidence intervals for equation parameters.
#'
#' @param equation Fitted equation object.
#' @param data Original data.
#' @param response Name of response column.
#' @param derivative_col Alias for response.
#' @param n_boot Number of bootstrap samples.
#' @param conf_level Confidence level (default 0.95).
#' @param block_size Block size for block bootstrap (time series).
#'
#' @return Data frame with parameter estimates and confidence intervals.
#' @export
bootstrap_parameters <- function(equation, data, response = NULL,
                                 derivative_col = NULL,
                                 n_boot = 500,
                                 conf_level = 0.95,
                                 block_size = NULL) {
  
  # Compatibility
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  
  n <- nrow(data)
  if (is.null(block_size)) {
    block_size <- max(1, floor(sqrt(n)))
  }
  
  original_coefs <- stats::coef(equation)
  if (is.null(original_coefs) || length(original_coefs) == 0) {
    warning("No coefficients found in equation")
    return(data.frame(parameter = character(0), estimate = numeric(0),
                      se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)))
  }
  
  n_coef <- length(original_coefs)
  
  # Storage for bootstrap estimates
  boot_coefs <- matrix(NA, nrow = n_boot, ncol = n_coef)
  colnames(boot_coefs) <- names(original_coefs)
  
  # Get expression for refitting
  if (inherits(equation, "symbolic_equation")) {
    eq_expr <- if (!is.null(equation$expression)) equation$expression else equation$string
  } else if (inherits(equation, "nls")) {
    eq_form <- stats::formula(equation)
  } else if (inherits(equation, "lm")) {
    eq_form <- stats::formula(equation)
  }
  
  for (b in 1:n_boot) {
    # Block bootstrap indices
    boot_idx <- block_bootstrap_indices(n, block_size)
    boot_data <- data[boot_idx, , drop = FALSE]
    
    # Refit on bootstrap sample
    tryCatch({
      if (inherits(equation, "lm")) {
        fit <- stats::lm(eq_form, data = boot_data)
        boot_coefs[b, ] <- stats::coef(fit)[names(original_coefs)]
      } else if (inherits(equation, "nls")) {
        if (!requireNamespace("minpack.lm", quietly = TRUE)) {
          stop("minpack.lm required")
        }
        fit <- minpack.lm::nlsLM(eq_form, data = boot_data, 
                                 start = as.list(original_coefs),
                                 control = minpack.lm::nls.lm.control(maxiter = 100))
        boot_coefs[b, ] <- stats::coef(fit)[names(original_coefs)]
      } else if (inherits(equation, "symbolic_equation")) {
        fit <- fit_specified_equation(eq_expr, data = boot_data,
                                      response = response,
                                      start = as.list(original_coefs))
        boot_coefs[b, ] <- stats::coef(fit)[names(original_coefs)]
      }
    }, error = function(e) {
      # Skip failed fits
    })
  }
  
  # Compute confidence intervals
  alpha <- 1 - conf_level
  results <- data.frame(
    parameter = names(original_coefs),
    estimate = original_coefs,
    se = apply(boot_coefs, 2, stats::sd, na.rm = TRUE),
    ci_lower = apply(boot_coefs, 2, stats::quantile, alpha/2, na.rm = TRUE),
    ci_upper = apply(boot_coefs, 2, stats::quantile, 1 - alpha/2, na.rm = TRUE),
    n_successful = colSums(!is.na(boot_coefs)),
    stringsAsFactors = FALSE
  )
  
  rownames(results) <- NULL
  
  return(results)
}


#' Block Bootstrap Indices
#' @keywords internal
block_bootstrap_indices <- function(n, block_size) {
  n_blocks <- ceiling(n / block_size)
  block_starts <- sample(1:(n - block_size + 1), n_blocks, replace = TRUE)
  
  indices <- unlist(lapply(block_starts, function(s) s:(s + block_size - 1)))
  indices <- indices[indices <= n]
  
  if (length(indices) < n) {
    indices <- c(indices, sample(1:n, n - length(indices), replace = TRUE))
  } else if (length(indices) > n) {
    indices <- indices[1:n]
  }
  
  indices
}


#' Coefficient Change Between Equations
#' @keywords internal
coefficient_change <- function(eq1, eq2) {
  c1 <- stats::coef(eq1)
  c2 <- stats::coef(eq2)
  
  if (is.null(c1) || is.null(c2)) return(Inf)
  
  common <- intersect(names(c1), names(c2))
  if (length(common) == 0) return(Inf)
  
  sqrt(mean((c1[common] - c2[common])^2))
}