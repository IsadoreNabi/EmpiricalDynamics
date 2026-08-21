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
#' @param weights Optional vector of observation weights, one per row of
#'   \code{data}. Each fold is refitted with the weights of its own training
#'   rows, so that a search run under weights is validated under them too.
#' @param refit_engine Which engine re-estimates the constants on each training
#'   fold: \code{"r"} (the default) uses this package's nonlinear least squares;
#'   \code{"julia"} hands the equation to SymbolicRegression.jl's own constant
#'   optimisation, which is the optimiser and the loss the search itself used,
#'   and requires Julia and the \code{JuliaCall} package. The default is
#'   \code{"r"} so that the figure a given equation gets does not depend on
#'   what happens to be installed on the machine.
#' @param verbose Print progress.
#'
#' @details A discovered equation arrives with its constants already fitted, as
#'   numeric literals, and there is nothing in it left to estimate. Before the
#'   folds begin, such an equation is passed through
#'   \code{\link{parameterize_equation}}: its literals become free parameters
#'   starting at the values the search found, and each fold re-estimates them on
#'   its own training rows. An equation specified by the user, which already has
#'   named free parameters, keeps taking its starting values from its own
#'   coefficients. An equation with no constants at all is not refitted -- there
#'   is nothing to refit -- and is scored by evaluating it on the held-out rows;
#'   the returned object says so in \code{refitted}.
#'
#' @return Object of class "cv_result" containing:
#'   \item{rmse}{Root mean squared error per fold}
#'   \item{mae}{Mean absolute error per fold}
#'   \item{r_squared}{R-squared per fold}
#'   \item{mean_rmse}{Average RMSE across folds}
#'   \item{sd_rmse}{Standard deviation of RMSE}
#'   \item{predictions}{List of predicted vs actual per fold}
#'   \item{fold_indices}{Indices used for each fold}
#'   \item{refitted}{Whether the equation was re-estimated on each fold}
#'   \item{parameterization}{The parameterised expression that was refitted, when
#'     the equation was a discovered one, and \code{NULL} otherwise}
#'   \item{refit_engine}{The engine that did the refitting}
#'   \item{stalled_folds}{Folds in which the Julia optimiser did not move off
#'     its starting constants, so the error reported for them is that of the
#'     equation as the search left it. Empty for the R engine.}
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
                           weights = NULL,
                           refit_engine = c("r", "julia"),
                           verbose = TRUE) {

  # Compatibility: derivative_col is alias for response
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  if (is.null(response)) {
    stop("Must specify 'response' or 'derivative_col'", call. = FALSE)
  }

  method <- match.arg(method)
  refit_engine <- match.arg(refit_engine)
  equation <- as_symbolic_equation(equation)
  n <- nrow(data)

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != n) {
      stop("'weights' must be a numeric vector with one value per row of 'data'.",
           call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must be finite and non-negative.", call. = FALSE)
    }
  }
  
  # Generate fold indices
  fold_indices <- switch(method,
                         "random" = create_random_folds(n, k),
                         "block" = create_block_folds(n, k, block_size),
                         "rolling" = create_rolling_folds(n, k, horizon)
  )
  
  # Extract formula from equation
  if (inherits(equation, "symbolic_equation")) {
    eq_expr <- if (!is.null(equation$expression)) {
      equation$expression
    } else {
      equation$string
    }
    if (is.null(eq_expr)) {
      stop("The equation carries neither an expression nor a string.", call. = FALSE)
    }
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
  
  # How the equation will be re-estimated, decided once rather than per fold.
  plan <- if (eq_type == "symbolic") {
    cv_refit_plan(eq_expr, equation, data)
  } else {
    list(kind = eq_type, expression = NULL, start = NULL,
         n_parameters = NA_integer_)
  }

  if (identical(refit_engine, "julia") && eq_type != "symbolic") {
    stop("'refit_engine = \"julia\"' applies to symbolic equations only.",
         call. = FALSE)
  }

  # Storage for results
  results <- list(
    rmse = numeric(k),
    mae = numeric(k),
    r_squared = numeric(k),
    predictions = vector("list", k),
    fold_indices = fold_indices
  )

  # Folds where the Julia optimiser did not move off its starting point.
  stalled <- integer(0)

  for (i in seq_len(k)) {
    if (verbose) message("Fold ", i, "/", k)

    test_idx <- fold_indices[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)

    train_data <- data[train_idx, , drop = FALSE]
    test_data <- data[test_idx, , drop = FALSE]
    train_weights <- if (is.null(weights)) NULL else weights[train_idx]

    # Refit on training data and predict on the held-out rows
    pred <- tryCatch({
      if (eq_type == "lm") {
        stats::predict(stats::lm(form, data = train_data), newdata = test_data)

      } else if (eq_type == "nls") {
        if (!requireNamespace("minpack.lm", quietly = TRUE)) {
          stop("Package 'minpack.lm' is required for NLS fitting.")
        }
        # For NLS, need starting values
        start_vals <- as.list(stats::coef(equation))
        args <- list(formula = form, data = train_data, start = start_vals,
                     control = minpack.lm::nls.lm.control(maxiter = 200))
        if (!is.null(train_weights)) args$weights <- train_weights
        stats::predict(do.call(minpack.lm::nlsLM, args), newdata = test_data)

      } else if (identical(refit_engine, "julia")) {
        got <- ed_julia_refit(eq_expr, train_data = train_data,
                              test_data = test_data, response = response,
                              weights = train_weights)
        if (!isTRUE(got$ok)) stop(got$message, call. = FALSE)
        if (identical(plan$kind, "discovered") && !isTRUE(got$moved)) {
          stalled <- c(stalled, i)
        }
        got$predictions

      } else if (identical(plan$kind, "fixed")) {
        # Nothing to estimate: the equation is scored as it stands.
        rep_len(as.numeric(eval(parse(text = plan$expression),
                                envir = create_eval_env(test_data))),
                nrow(test_data))

      } else {
        fitted_model <- fit_specified_equation(
          plan$expression,
          data = train_data,
          response = response,
          start = plan$start,
          method = cv_fit_method(plan$expression, train_data),
          weights = train_weights
        )
        stats::predict(fitted_model, newdata = test_data)
      }
    }, error = function(e) {
      warning("Fold ", i, " fitting failed: ", conditionMessage(e))
      NULL
    })

    if (is.null(pred)) {
      results$rmse[i] <- NA
      results$mae[i] <- NA
      results$r_squared[i] <- NA
      next
    }

    pred <- as.numeric(pred)
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
  results$refitted <- !identical(plan$kind, "fixed")
  results$parameterization <- if (identical(plan$kind, "discovered")) {
    plan$expression
  } else {
    NULL
  }
  results$refit_engine <- if (eq_type == "symbolic") refit_engine else "r"

  # Reported rather than hidden: the search's constant optimisation perturbs
  # its restarts multiplicatively, so a constant the search left at zero cannot
  # move, and the error below is the error of the equation as it stands, not of
  # the equation refitted.
  results$stalled_folds <- stalled
  if (length(stalled) > 0L) {
    warning("The Julia optimiser did not move off its starting constants in ",
            length(stalled), " of ", k, " folds; the error reported for those ",
            "folds is that of the equation as the search left it")
  }

  class(results) <- "cv_result"
  return(results)
}


#' Accept an equation in any of the shapes this package hands out
#'
#' \code{symbolic_search()} returns its candidates in \code{all_equations} as
#' plain lists -- a string, a complexity, an error -- while
#' \code{get_pareto_set()} and \code{select_equation()} wrap the same content in
#' a \code{symbolic_equation}. Validation used to insist on the second and
#' refuse the first with "Unknown equation type", so the records the search
#' itself produces could not be handed to the package's own cross-validation.
#' Anything carrying an expression is now accepted and classed here.
#'
#' @keywords internal
#' @noRd
as_symbolic_equation <- function(equation) {
  if (inherits(equation, "symbolic_equation")) return(equation)
  if (!is.list(equation) || is.data.frame(equation)) return(equation)

  expression <- equation$expression
  if (is.null(expression)) expression <- equation$string
  if (is.null(expression) || !is.character(expression) ||
      length(expression) != 1L || !nzchar(expression)) {
    return(equation)
  }

  equation$expression <- expression
  if (is.null(equation$string)) equation$string <- expression
  class(equation) <- c("symbolic_equation", class(equation))
  equation
}


#' Decide how an equation is to be re-estimated on each fold
#'
#' Three cases. An equation whose expression mentions names the data does not
#' carry has free parameters already and is refitted as it stands. A discovered
#' equation mentions only data columns, because its constants are literals: it
#' is parameterised, and the literals become the starting values. An equation
#' that mentions only data columns and carries no literal has nothing to
#' estimate.
#'
#' @keywords internal
#' @noRd
cv_refit_plan <- function(eq_expr, equation, data) {
  free <- setdiff(all.vars(parse(text = eq_expr)), names(data))

  if (length(free) > 0L) {
    coefs <- stats::coef(equation)
    start <- if (!is.null(coefs) && length(coefs) > 0L) as.list(coefs) else NULL
    return(list(kind = "specified", expression = eq_expr, start = start,
                n_parameters = length(free)))
  }

  parameterised <- parameterize_equation(eq_expr, data = data)
  if (parameterised$n_parameters == 0L) {
    return(list(kind = "fixed", expression = eq_expr, start = list(),
                n_parameters = 0L))
  }

  list(kind = "discovered", expression = parameterised$expression,
       start = parameterised$start, n_parameters = parameterised$n_parameters)
}


#' Which optimiser to refit an expression with
#'
#' An equation that mentions no variable at all -- a bare constant, which the
#' search does propose at complexity one -- has a one-by-one gradient, and the
#' nonlinear least squares routines recycle it against the full weight vector.
#' Such an equation is fitted by general optimisation instead, where a scalar
#' prediction is expected.
#'
#' @keywords internal
#' @noRd
cv_fit_method <- function(expression, data) {
  mentioned <- intersect(all.vars(parse(text = expression)), names(data))
  if (length(mentioned) == 0L) "optim" else "LM"
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
  cat("Method:", x$method, "with", x$k, "folds\n")
  if (!is.null(x$refitted)) {
    if (isTRUE(x$refitted)) {
      cat("Refitted per fold with the", x$refit_engine, "engine")
      if (!is.null(x$parameterization)) {
        cat(" as:", x$parameterization)
      }
      cat("\n")
    } else {
      cat("Not refitted: the equation carries no constant to estimate\n")
    }
  }
  cat("\n")
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

# Extract the pieces every qualitative routine needs from a fitted equation:
# how it must be evaluated, which symbols it consults, and which coefficients
# it carries. "expression" equations (symbolic_equation, nls) evaluate their
# formula text in an environment; "model" equations (lm, and glm through
# inheritance) must be evaluated by their own predict() machinery, because
# their coefficient names -- "I(Z)", "(Intercept)", factor contrasts -- are
# term labels, not symbols of the formula text, and injecting them into an
# evaluation environment binds nothing the expression ever consults.
ed_equation_pieces <- function(equation) {
  if (inherits(equation, "symbolic_equation")) {
    expr_str <- if (!is.null(equation$expression)) {
      equation$expression
    } else {
      equation$string
    }
    list(kind = "expression", expr_str = expr_str,
         symbols = all.vars(str2lang(expr_str)),
         coefs = stats::coef(equation))
  } else if (inherits(equation, "nls")) {
    # nls parameter names are genuine symbols of its formula, so the
    # expression path is exact for it.
    expr_str <- paste(deparse(stats::formula(equation)[[3]]), collapse = " ")
    list(kind = "expression", expr_str = expr_str,
         symbols = all.vars(str2lang(expr_str)),
         coefs = stats::coef(equation))
  } else if (inherits(equation, "lm")) {
    trm <- stats::delete.response(stats::terms(equation))
    list(kind = "model", expr_str = NULL,
         symbols = all.vars(trm),
         coefs = stats::coef(equation))
  } else {
    stop("Unknown equation type: expected a 'symbolic_equation', 'nls', ",
         "'lm' or 'glm' object, got ",
         paste(class(equation), collapse = "/"), call. = FALSE)
  }
}

# A declared value that carries no name binds nothing: names(x) is NULL and
# the injection loop runs zero times, silently. Refuse that up front.
ed_check_named_list <- function(x, what) {
  if (length(x) == 0L) return(invisible(NULL))
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("every element of '", what, "' must be named: an unnamed value ",
         "cannot bind to any symbol of the equation and would be ignored ",
         "silently", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("'", what, "' declares the name '",
         nms[anyDuplicated(nms)], "' more than once", call. = FALSE)
  }
  invisible(NULL)
}

#' Analyze Fixed Points
#'
#' Finds and characterizes fixed points of the discovered equation.
#'
#' For \code{lm} and \code{glm} equations the right-hand side is evaluated
#' through the model's own \code{predict()} machinery (on the response scale
#' for \code{glm}), so the fitted coefficients enter exactly as fitted --
#' \code{I()} terms, factors and link functions included. For
#' \code{symbolic_equation} and \code{nls} equations the expression is
#' evaluated with the fitted coefficients bound by name.
#'
#' A fitted object that carries posterior draws -- an \code{rstanarm} fit, for
#' instance, whose class is \code{stanreg, glm, lm} -- is swept through its
#' draws instead of through its point estimate, and the result then carries
#' one block of rows per draw with a \code{draw} column identifying them; see
#' the \code{n_draws} argument and the \code{"posterior"} attribute of the
#' return value. Objects are recognised by whether they can produce a draws
#' matrix covering their own coefficients, not by class name.
#'
#' When \code{coefficient_values} is used, the substitution is verified rather
#' than assumed: the requested override must actually move the predictions
#' somewhere on the search grid. An object whose \code{predict()} does not read
#' \code{$coefficients} -- which is true of every fit that stores its estimate
#' elsewhere and merely summarises it into \code{$coefficients} -- is refused
#' with an \code{"ed_substitution_ignored"} condition instead of being swept
#' into a table of identical copies. A coefficient whose term is identically
#' zero over the search range is refused as \code{"ed_substitution_inert"}.
#'
#' Every name consulted by the equation must be accounted for: the search
#' \code{variable}, a fitted coefficient, or a value declared in
#' \code{exogenous_values}. A free symbol without a value, an unnamed
#' declaration, a declaration that collides with the search variable or with
#' a fitted coefficient -- each of these is reported instead of being
#' silently ignored or silently overwriting the fit. To vary a fitted
#' coefficient deliberately, use \code{coefficient_values}.
#'
#' @param equation Fitted equation object (\code{symbolic_equation},
#'   \code{nls}, \code{lm} or \code{glm}), or any object inheriting from
#'   \code{lm} that either honours a coefficient substitution or exposes
#'   posterior draws covering its own coefficients.
#' @param variable Name of the main variable.
#' @param range Numeric vector of length 2 specifying search range.
#' @param n_grid Number of grid points for initial search.
#' @param exogenous_values Named list of fixed values for exogenous variables.
#' @param coefficient_values Named list of deliberate overrides for fitted
#'   coefficients; every name must match a fitted coefficient.
#' @param n_draws Number of posterior draws to sweep when \code{equation}
#'   carries a posterior; ignored otherwise. Draws are thinned at even
#'   spacing, so the selection is reproducible without a seed. The default of
#'   200 places the Monte Carlo error of a posterior median near 0.09
#'   posterior standard deviations and of a 5% or 95% quantile near 0.15.
#'
#' @return Data frame of fixed points with stability classification. The
#'   values any exogenous variables were held at travel with the result as
#'   its \code{"exogenous_values"} attribute, because a fixed point found
#'   under fixed exogenous values is conditional on that fixing; the
#'   deliberate coefficient overrides travel as \code{"coefficient_values"}.
#'   Both attributes are attached whether or not any fixed point was found.
#'   The \code{"posterior"} attribute says whether the sweep was run over
#'   posterior draws; when it is \code{TRUE} the frame has a \code{draw}
#'   column, holds one block of rows per draw, and \code{"n_draws"} gives how
#'   many draws were swept. \code{ed_consensus_fixed_points()} reduces such a
#'   frame to one row per branch.
#' @examples
#' # dZ = 2*Z - Z^2: fixed points at 0 (unstable, f' = +2)
#' #                            and at 2 (stable,   f' = -2)
#' data <- data.frame(Z = seq(0.1, 3, length.out = 50))
#' data$dZ <- 2 * data$Z - data$Z^2
#' model <- stats::lm(dZ ~ I(Z) + I(Z^2) + 0, data = data)
#' analyze_fixed_points(model, variable = "Z", range = c(-1, 3))
#' @export
analyze_fixed_points <- function(equation, variable,
                                 range = c(-10, 10),
                                 n_grid = 100,
                                 exogenous_values = list(),
                                 coefficient_values = list(),
                                 n_draws = 200L) {

  if (!is.character(variable) || length(variable) != 1L || !nzchar(variable)) {
    stop("'variable' must be a single variable name", call. = FALSE)
  }
  if (!is.numeric(range) || length(range) != 2L || !all(is.finite(range)) ||
      range[1] >= range[2]) {
    stop("'range' must be two finite numbers in increasing order",
         call. = FALSE)
  }
  if (!is.numeric(n_draws) || length(n_draws) != 1L || is.na(n_draws) ||
      n_draws < 1) {
    stop("'n_draws' must be a single number of posterior draws, at least 1",
         call. = FALSE)
  }
  ed_check_named_list(exogenous_values, "exogenous_values")
  ed_check_named_list(coefficient_values, "coefficient_values")

  pieces <- ed_equation_pieces(equation)
  coefs <- pieces$coefs

  if (variable %in% names(exogenous_values)) {
    stop("'exogenous_values' declares the search variable '", variable,
         "': the declared value would be overwritten at every grid point ",
         "and bind nothing", call. = FALSE)
  }
  if (!(variable %in% pieces$symbols)) {
    stop("variable '", variable, "' does not appear in the equation ",
         "(its variables are: ", paste(pieces$symbols, collapse = ", "),
         "); the equation would be constant in it", call. = FALSE)
  }
  if (length(coefficient_values)) {
    unknown <- setdiff(names(coefficient_values), names(coefs))
    if (length(unknown)) {
      stop("'coefficient_values' names ",
           paste0("'", unknown, "'", collapse = ", "),
           " which are not fitted coefficients of this equation ",
           "(fitted: ", paste(names(coefs), collapse = ", "), ")",
           call. = FALSE)
    }
  }
  unused <- setdiff(names(exogenous_values), pieces$symbols)
  if (length(unused)) {
    warning("'exogenous_values' declares ",
            paste0("'", unused, "'", collapse = ", "),
            ", which the equation never consults: the declaration has ",
            "no effect", call. = FALSE)
  }

  grid <- seq(range[1], range[2], length.out = n_grid)
  draw_f <- NULL
  vals_grid <- NULL
  b_draws <- NULL

  if (pieces$kind == "model") {
    # lm / glm: evaluate through the model's own prediction machinery.
    needed <- setdiff(pieces$symbols, c(variable, names(exogenous_values)))
    if (length(needed)) {
      stop("the model consults ",
           paste0("'", needed, "'", collapse = ", "),
           " but no value was declared for ",
           if (length(needed) == 1L) "it" else "them",
           ": pass fixed values through 'exogenous_values'", call. = FALSE)
    }
    is_glm <- inherits(equation, "glm")
    # One evaluator, built from whichever object it is handed: the fit as it
    # came, or the fit with the deliberate overrides written into it. The
    # verification below compares the two through this same route, so what it
    # certifies is the route the sweep actually uses.
    make_f <- function(m) {
      function(z) {
        nd <- stats::setNames(data.frame(z = z), variable)
        for (nm in names(exogenous_values)) {
          nd[[nm]] <- exogenous_values[[nm]]
        }
        tryCatch({
          if (is_glm) {
            as.numeric(stats::predict(m, newdata = nd, type = "response"))
          } else {
            as.numeric(stats::predict(m, newdata = nd))
          }
        }, error = function(e) rep(NA_real_, length(z)))
      }
    }
    design <- ed_design_of(equation, variable, exogenous_values)
    b_draws <- ed_draws_of(equation)
    # An override that asks for the value already fitted moves nothing by
    # arithmetic and must not be mistaken for an override that was ignored.
    effective <- ed_effective_overrides(coefficient_values, coefs)

    if (is.null(b_draws)) {
      substituted <- equation
      if (length(coefficient_values)) {
        for (nm in names(coefficient_values)) {
          substituted$coefficients[nm] <- as.numeric(coefficient_values[[nm]])
        }
      }
      if (length(effective)) {
        ed_verify_substitution(make_f, equation, substituted, effective,
                               design, grid, equation)
      }
      f_vec <- make_f(substituted)
      f <- function(z) f_vec(z)
    } else {
      linkinv <- ed_linkinv_of(equation, is_glm)
      if (is.null(linkinv)) {
        ed_stop("ed_uncertifiable_substitution", paste0(
          "this ", paste(class(equation), collapse = "/"), " object carries ",
          "posterior draws but exposes no family, so a per draw prediction ",
          "cannot be put on the same scale as its own predict()"))
      }
      ed_certify_design(make_f(equation), design, b_draws, linkinv, grid,
                        equation)
      if (length(effective) &&
          isTRUE(ed_column_is_inert(design, names(effective), grid))) {
        ed_stop_inert(names(effective))
      }
      b_draws <- ed_thin_draws(b_draws, n_draws)
      for (nm in names(coefficient_values)) {
        b_draws[, nm] <- as.numeric(coefficient_values[[nm]])
      }
      cn <- colnames(b_draws)
      vals_grid <- linkinv(design(grid)[, cn, drop = FALSE] %*% t(b_draws))
      draw_f <- function(s) {
        b <- b_draws[s, ]
        function(z) {
          x <- tryCatch(design(z), error = function(e) NULL)
          if (is.null(x)) return(rep(NA_real_, length(z)))
          as.numeric(linkinv(x[, cn, drop = FALSE] %*% b))
        }
      }
    }
  } else {
    expr_str <- pieces$expr_str
    clash <- intersect(names(exogenous_values), names(coefs))
    if (length(clash)) {
      stop("'exogenous_values' declares ",
           paste0("'", clash, "'", collapse = ", "),
           ", which ", if (length(clash) == 1L) "is a" else "are",
           " fitted coefficient", if (length(clash) == 1L) "" else "s",
           " of the equation: the declaration would silently overwrite ",
           "the fit. To vary a fitted coefficient deliberately, use ",
           "'coefficient_values'", call. = FALSE)
    }
    if (length(coefficient_values)) {
      for (nm in names(coefficient_values)) {
        coefs[nm] <- as.numeric(coefficient_values[[nm]])
      }
    }
    known <- c(variable, names(coefs), names(exogenous_values))
    free <- setdiff(pieces$symbols, known)
    free <- free[!vapply(free, exists, logical(1), envir = baseenv())]
    if (length(free)) {
      stop("the equation consults ",
           paste0("'", free, "'", collapse = ", "),
           " but no value was declared for ",
           if (length(free) == 1L) "it" else "them",
           ": pass fixed values through 'exogenous_values'", call. = FALSE)
    }

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
    f_vec <- function(z) sapply(z, f)
  }

  if (is.null(draw_f)) {
    found <- ed_fixed_points_from(f, f_vec(grid), grid, range, n_grid)
    reason <- attr(found, "ed_reason")
    result <- data.frame(fixed_point = found$fixed_point,
                         stability = found$stability,
                         eigenvalue = found$eigenvalue)
    n_used <- NA_integer_
  } else {
    # One sweep per posterior draw, each found by the same core the point
    # estimate goes through. The reasons are collected and reported once:
    # a caller sweeping hundreds of draws must not be told hundreds of times.
    n_used <- nrow(b_draws)
    per <- vector("list", n_used)
    reasons <- character(n_used)
    for (s in seq_len(n_used)) {
      one <- ed_fixed_points_from(draw_f(s), vals_grid[, s], grid, range,
                                  n_grid)
      reasons[s] <- attr(one, "ed_reason")
      if (nrow(one)) {
        one$draw <- s
        per[[s]] <- one
      }
    }
    result <- do.call(rbind, per)
    if (is.null(result)) {
      result <- data.frame(fixed_point = numeric(0), stability = character(0),
                           eigenvalue = numeric(0), draw = integer(0))
    }
    reason <- if (all(reasons == "unevaluable")) {
      "unevaluable"
    } else if (!any(reasons == "ok")) {
      "none"
    } else {
      "ok"
    }
  }
  if (identical(reason, "unevaluable")) {
    message("Could not evaluate function on grid")
  } else if (identical(reason, "none")) {
    message("No fixed points found in specified range")
  }

  # A fixed point found with exogenous variables held fixed is conditional
  # on that fixing: the declaration travels with the result -- and it travels
  # with an empty result too, because "no fixed point under W = 2" is as
  # conditional a statement as any other.
  rownames(result) <- NULL
  attr(result, "variable") <- variable
  attr(result, "exogenous_values") <- exogenous_values
  attr(result, "coefficient_values") <- coefficient_values
  attr(result, "posterior") <- !is.null(draw_f)
  attr(result, "n_draws") <- n_used
  return(result)
}

#' Analyze Bifurcations
#'
#' Examines how fixed points change as a parameter varies.
#'
#' The parameter must actually appear in the equation, and is varied
#' according to what it is: a fitted coefficient is overridden through
#' \code{coefficient_values} of \code{\link{analyze_fixed_points}}, and a
#' free variable of the expression is bound through
#' \code{exogenous_values}. A name that is neither is an error -- sweeping
#' it would evaluate the same equation \code{n_param} times and label the
#' identical copies with the parameter's name.
#'
#' When one name could denote both: for \code{lm}/\code{glm} equations a
#' plain term's coefficient label coincides with the variable's name, and
#' the name is taken as the variable (labels that are never variables, like
#' \code{"I(Z)"}, denote their coefficient); for \code{nls} and
#' \code{symbolic_equation} equations a coefficient name is a symbol of the
#' formula itself and is taken as the coefficient.
#'
#' When \code{equation} carries posterior draws, the sweep is run once per
#' draw and the answer is a distribution of fixed points at every parameter
#' value rather than a point. This is the interventional reading of the
#' sweep -- the coordinate is SET to each value in every draw, the remaining
#' uncertainty left as estimated -- and not the conditional one, which would
#' keep only the draws that happen to agree with the value and would empty out
#' at exactly the extreme parameter values a bifurcation diagram exists to
#' visit. \code{$data} then carries a \code{draw} column and \code{$summary}
#' holds the branch-by-branch reduction.
#'
#' @param equation Fitted equation object: a \code{symbolic_equation}, an
#'   \code{nls}, an \code{lm} or a \code{glm}. An object inheriting from
#'   \code{lm} is accepted when it either honours a coefficient substitution
#'   or exposes posterior draws covering its own coefficients -- an
#'   \code{rstanarm} fit does the latter. One that does neither is refused
#'   rather than swept into identical copies; see
#'   \code{\link{analyze_fixed_points}}.
#' @param variable Name of the main variable.
#' @param parameter Name of the parameter to vary. Must be a fitted
#'   coefficient or a variable of the equation, and must differ from
#'   \code{variable}.
#' @param param_range Range for parameter values.
#' @param n_param Number of parameter values to test.
#' @param z_range Range for searching fixed points.
#' @param exogenous_values Fixed values for other variables.
#' @param n_draws Number of posterior draws to sweep at each parameter value
#'   when \code{equation} carries a posterior; ignored otherwise. See
#'   \code{\link{analyze_fixed_points}}.
#' @param interval Width of the posterior interval reported in
#'   \code{$summary}; ignored when \code{equation} carries no posterior. The
#'   default of 0.9 is the width \code{rstanarm} itself reports, chosen there
#'   because the 2.5% and 97.5% quantiles need several times as many draws
#'   for the same Monte Carlo error as the 5% and 95% ones.
#'
#' @return Object of class "bifurcation_analysis". \code{$data} holds the
#'   fixed points found; \code{$status} holds one row per requested
#'   parameter value -- \code{"ok"}, \code{"no_fixed_points"} or
#'   \code{"error"} with the error's message -- so the requested grid can
#'   always be reconstructed and a failure is never mistaken for an
#'   absence; \code{$mode} records whether the parameter was varied as a
#'   \code{"coefficient"} or as an \code{"exogenous"} variable.
#'   \code{$posterior} says whether the sweep carries uncertainty. When it
#'   does, \code{$data} gains a \code{draw} column and holds one block of rows
#'   per draw per parameter value; \code{$status} gains \code{n_draws} and
#'   \code{prob_n_fixed_points}, and its \code{n_fixed_points} is the count
#'   most draws agreed on rather than the number of rows; and \code{$summary}
#'   holds one row per parameter value per branch, with
#'   \code{fixed_point_lower} and \code{fixed_point_upper} bounding a
#'   \code{$interval} posterior interval. \code{$summary} is \code{NULL} for a
#'   point estimate, where \code{$data} already is the answer.
#' @examples
#' # dZ = r*Z - Z^2: transcritical bifurcation at r = 0
#' set.seed(1)
#' data <- data.frame(Z = seq(-1, 3, length.out = 80))
#' data$dZ <- 2 * data$Z - data$Z^2 + 0.01 * stats::rnorm(80)
#' fit <- stats::nls(dZ ~ r * Z - Z^2, data = data, start = list(r = 1))
#' bif <- analyze_bifurcations(fit, variable = "Z", parameter = "r",
#'                             param_range = c(-1, 1), n_param = 11,
#'                             z_range = c(-3, 3))
#' print(bif)
#' @export
analyze_bifurcations <- function(equation, variable, parameter,
                                 param_range = c(-5, 5),
                                 n_param = 50,
                                 z_range = c(-10, 10),
                                 exogenous_values = list(),
                                 n_draws = 200L,
                                 interval = 0.9) {

  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval) ||
      interval <= 0 || interval >= 1) {
    stop("'interval' must be a single number strictly between 0 and 1",
         call. = FALSE)
  }
  if (!is.character(parameter) || length(parameter) != 1L ||
      !nzchar(parameter)) {
    stop("'parameter' must be a single name", call. = FALSE)
  }
  if (identical(parameter, variable)) {
    stop("'parameter' and 'variable' are both '", variable,
         "': the bifurcation parameter cannot be the variable whose ",
         "fixed points are being tracked", call. = FALSE)
  }
  ed_check_named_list(exogenous_values, "exogenous_values")
  if (parameter %in% names(exogenous_values)) {
    stop("'exogenous_values' already fixes '", parameter,
         "', the parameter being varied", call. = FALSE)
  }

  pieces <- ed_equation_pieces(equation)
  # For an lm/glm a plain term's coefficient label coincides with the
  # variable's own name; there the name denotes the variable, and labels
  # that are never variables -- "I(Z)", "(Intercept)" -- denote their
  # coefficient. For an expression (nls, symbolic_equation) a coefficient
  # name is a genuine symbol of the formula and denotes the coefficient.
  candidates <- if (pieces$kind == "model") {
    list(exogenous = pieces$symbols, coefficient = names(pieces$coefs))
  } else {
    list(coefficient = names(pieces$coefs), exogenous = pieces$symbols)
  }
  mode <- if (parameter %in% candidates[[1L]]) {
    names(candidates)[1L]
  } else if (parameter %in% candidates[[2L]]) {
    names(candidates)[2L]
  } else {
    stop("parameter '", parameter, "' appears nowhere in the equation: ",
         "it is neither a fitted coefficient (",
         paste(names(pieces$coefs), collapse = ", "),
         ") nor a variable (",
         paste(pieces$symbols, collapse = ", "),
         "). Sweeping it would return n_param identical copies of the ",
         "same fixed points", call. = FALSE)
  }

  # Names the equation never consults are diagnosed once, here, instead of
  # once per parameter value inside the loop.
  unused <- setdiff(names(exogenous_values), pieces$symbols)
  if (length(unused)) {
    warning("'exogenous_values' declares ",
            paste0("'", unused, "'", collapse = ", "),
            ", which the equation never consults: dropped", call. = FALSE)
    exogenous_values <- exogenous_values[setdiff(names(exogenous_values),
                                                 unused)]
  }

  param_vals <- seq(param_range[1], param_range[2], length.out = n_param)

  all_fps <- vector("list", length(param_vals))
  all_sum <- vector("list", length(param_vals))
  status <- data.frame(
    parameter_value = param_vals,
    n_fixed_points = NA_integer_,
    status = NA_character_,
    message = NA_character_,
    n_draws = NA_integer_,
    prob_n_fixed_points = NA_real_
  )
  posterior <- FALSE
  n_used <- NA_integer_

  for (i in seq_along(param_vals)) {
    fps <- tryCatch({
      if (mode == "coefficient") {
        suppressMessages(analyze_fixed_points(
          equation, variable, range = z_range,
          exogenous_values = exogenous_values,
          coefficient_values = stats::setNames(list(param_vals[i]),
                                               parameter),
          n_draws = n_draws
        ))
      } else {
        exog <- exogenous_values
        exog[[parameter]] <- param_vals[i]
        suppressMessages(analyze_fixed_points(
          equation, variable, range = z_range,
          exogenous_values = exog, n_draws = n_draws
        ))
      }
    }, error = function(e) e)

    if (inherits(fps, "ed_error")) {
      # Whether this object honours a coefficient substitution at all, and
      # whether the swept term can move anything, are properties of the
      # object and the range -- not of the parameter value. Recording them
      # per value would bury the one fact that matters under n_param copies
      # of it, and hand back an object whose $status is entirely failures.
      stop(fps)
    }
    if (inherits(fps, "error")) {
      # A failure is recorded as a failure: it must never be mistaken for
      # a parameter value that simply has no fixed points.
      status$status[i] <- "error"
      status$message[i] <- conditionMessage(fps)
      next
    }

    if (isTRUE(attr(fps, "posterior"))) {
      posterior <- TRUE
      n_used <- attr(fps, "n_draws")
      cons <- ed_consensus_fixed_points(fps, n_used, interval)
      status$n_fixed_points[i] <- if (nrow(cons)) {
        cons$n_fixed_points[1]
      } else {
        0L
      }
      status$n_draws[i] <- n_used
      status$prob_n_fixed_points[i] <- attr(cons, "prob_n_fixed_points")
      if (nrow(cons) > 0) {
        cons$parameter_value <- param_vals[i]
        cons$parameter_name <- parameter
        all_sum[[i]] <- cons
      }
    } else {
      status$n_fixed_points[i] <- nrow(fps)
    }
    status$status[i] <- if (nrow(fps) > 0) "ok" else "no_fixed_points"
    if (nrow(fps) > 0) {
      fps$parameter_value <- param_vals[i]
      fps$parameter_name <- parameter
      all_fps[[i]] <- fps
    }
  }

  # The two posterior columns are carried through the loop so that a value
  # can be written at any index, and dropped again when no draw was ever
  # swept: a point estimate's $status must keep the shape it always had.
  if (!posterior) {
    status$n_draws <- NULL
    status$prob_n_fixed_points <- NULL
  }

  bifurc_data <- do.call(rbind, all_fps)
  if (is.null(bifurc_data)) {
    bifurc_data <- data.frame(
      fixed_point = numeric(0), stability = character(0),
      eigenvalue = numeric(0), parameter_value = numeric(0),
      parameter_name = character(0)
    )
    if (posterior) bifurc_data$draw <- integer(0)
  }
  rownames(bifurc_data) <- NULL

  bifurc_summary <- NULL
  if (posterior) {
    bifurc_summary <- do.call(rbind, all_sum)
    if (is.null(bifurc_summary)) {
      bifurc_summary <- ed_consensus_fixed_points(
        bifurc_data[0, , drop = FALSE], 1L, interval)[0, , drop = FALSE]
      bifurc_summary$parameter_value <- numeric(0)
      bifurc_summary$parameter_name <- character(0)
    }
    rownames(bifurc_summary) <- NULL
  }

  result <- list(
    data = bifurc_data,
    summary = bifurc_summary,
    status = status,
    mode = mode,
    parameter = parameter,
    variable = variable,
    param_range = param_range,
    exogenous_values = exogenous_values,
    posterior = posterior,
    n_draws = n_used,
    interval = if (posterior) interval else NA_real_
  )

  class(result) <- "bifurcation_analysis"
  return(result)
}

#' Print a Bifurcation Analysis
#'
#' @param x Object of class bifurcation_analysis.
#' @param ... Additional arguments (ignored).
#'
#' @return \code{x}, invisibly.
#' @export
print.bifurcation_analysis <- function(x, ...) {
  cat("Bifurcation analysis: fixed points of '", x$variable,
      "' as '", x$parameter, "' varies (as a ",
      if (is.null(x$mode)) "parameter" else x$mode, ")\n", sep = "")
  if (!is.null(x$status)) {
    n_ok <- sum(x$status$status == "ok")
    n_empty <- sum(x$status$status == "no_fixed_points")
    n_err <- sum(x$status$status == "error")
    cat(sprintf(
      "  %d parameter values in [%g, %g]: %d with fixed points, %d without, %d failed\n",
      nrow(x$status), x$param_range[1], x$param_range[2],
      n_ok, n_empty, n_err))
    if (n_err > 0) {
      first_err <- which(x$status$status == "error")[1]
      cat("  first failure, at ", x$parameter, " = ",
          format(x$status$parameter_value[first_err]), ": ",
          x$status$message[first_err], "\n", sep = "")
    }
  }
  if (length(x$exogenous_values)) {
    cat("  exogenous variables held fixed: ",
        paste(names(x$exogenous_values),
              vapply(x$exogenous_values, format, character(1)),
              sep = " = ", collapse = ", "), "\n", sep = "")
  }
  if (isTRUE(x$posterior)) {
    # The count of rows is not the count of fixed points here: it is the
    # count of fixed points times the count of draws. Saying "fixed points
    # found: 8400" would be arithmetically true and read as a finding.
    cat("  swept over ", x$n_draws,
        " posterior draws per parameter value: each fixed point is a ",
        format(100 * x$interval), "% interval, not a number\n", sep = "")
    if (!is.null(x$summary) && nrow(x$summary) > 0) {
      cat("  branches summarised: ", nrow(x$summary), " (",
          paste(names(table(x$summary$stability)),
                table(x$summary$stability),
                sep = ": ", collapse = ", "), ")\n", sep = "")
      worst <- min(x$summary$prob_n_fixed_points)
      cat("  draws agreeing on the number of fixed points: ",
          format(100 * min(1, worst), digits = 3), "% at worst\n", sep = "")
    }
  } else if (!is.null(x$data) && nrow(x$data) > 0) {
    cat("  fixed points found: ", nrow(x$data), " (",
        paste(names(table(x$data$stability)),
              table(x$data$stability),
              sep = ": ", collapse = ", "), ")\n", sep = "")
  }
  invisible(x)
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
  # A posterior sweep returns one block of rows per draw. Counting those rows
  # would report the number of draws times the number of fixed points and
  # compare it to an expectation about fixed points. The posterior is reduced
  # to the count most draws agreed on, and how much of the posterior that was
  # is reported rather than assumed away.
  results$posterior <- isTRUE(attr(fps, "posterior"))
  if (results$posterior) {
    fps <- ed_consensus_fixed_points(fps, attr(fps, "n_draws"))
    results$messages <- c(results$messages, sprintf(
      "Posterior fit: checks are against the modal count over %d draws, held by %s of them",
      attr(fps, "n_draws"),
      paste0(format(100 * attr(fps, "prob_n_fixed_points"), digits = 3), "%")))
  }
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

  equation <- as_symbolic_equation(equation)

  # A discovered equation has its constants written in as literals and reports
  # no coefficients; the sensitivity being asked about is precisely the
  # sensitivity to those constants, so they are made into parameters first.
  if (inherits(equation, "symbolic_equation") &&
      length(stats::coef(equation)) == 0L) {
    eq_expr <- if (!is.null(equation$expression)) equation$expression else equation$string
    if (!is.null(eq_expr)) {
      plan <- cv_refit_plan(eq_expr, equation, data)
      if (identical(plan$kind, "discovered")) {
        equation$expression <- plan$expression
        equation$string <- plan$expression
        equation$coefficients <- unlist(plan$start)
      }
    }
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
#' @param weights Optional vector of observation weights, one per row of
#'   \code{data}.
#'
#' @details A discovered equation carries its constants as literals and has no
#'   coefficients to resample. It is therefore parameterised first, exactly as
#'   in \code{\link{cross_validate}}, and what is resampled are the parameters
#'   that stand in for those literals.
#'
#' @return Data frame with parameter estimates and confidence intervals.
#' @export
bootstrap_parameters <- function(equation, data, response = NULL,
                                 derivative_col = NULL,
                                 n_boot = 500,
                                 conf_level = 0.95,
                                 block_size = NULL,
                                 weights = NULL) {

  # Compatibility
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }

  equation <- as_symbolic_equation(equation)
  n <- nrow(data)
  if (is.null(block_size)) {
    block_size <- max(1, floor(sqrt(n)))
  }

  if (!is.null(weights) &&
      (!is.numeric(weights) || length(weights) != n)) {
    stop("'weights' must be a numeric vector with one value per row of 'data'.",
         call. = FALSE)
  }

  # Get expression for refitting, and the parameters to resample
  eq_expr <- NULL
  if (inherits(equation, "symbolic_equation")) {
    eq_expr <- if (!is.null(equation$expression)) equation$expression else equation$string
    plan <- cv_refit_plan(eq_expr, equation, data)
    if (identical(plan$kind, "fixed")) {
      warning("The equation carries no constant to resample")
      return(data.frame(parameter = character(0), estimate = numeric(0),
                        se = numeric(0), ci_lower = numeric(0),
                        ci_upper = numeric(0)))
    }
    eq_expr <- plan$expression
    original_coefs <- if (identical(plan$kind, "discovered")) {
      # The literal the search wrote down is a starting value, not the estimate
      # this resampling is about: the estimate is what the constant comes to on
      # the whole sample, which is what each replicate re-estimates. Taking the
      # literal instead would report an interval that need not even contain the
      # point it is drawn around.
      on_full <- tryCatch(
        stats::coef(fit_specified_equation(
          eq_expr, data = data, response = response,
          start = as.list(unlist(plan$start)),
          method = cv_fit_method(eq_expr, data), weights = weights)),
        error = function(e) NULL)
      if (is.null(on_full)) {
        warning("The equation could not be fitted on the full sample; ",
                "the constants found by the search are reported instead")
        unlist(plan$start)
      } else {
        on_full
      }
    } else {
      stats::coef(equation)
    }
  } else if (inherits(equation, "nls") || inherits(equation, "lm")) {
    eq_form <- stats::formula(equation)
    original_coefs <- stats::coef(equation)
  } else {
    original_coefs <- stats::coef(equation)
  }

  if (is.null(original_coefs) || length(original_coefs) == 0) {
    warning("No coefficients found in equation")
    return(data.frame(parameter = character(0), estimate = numeric(0),
                      se = numeric(0), ci_lower = numeric(0), ci_upper = numeric(0)))
  }

  n_coef <- length(original_coefs)

  # Storage for bootstrap estimates
  boot_coefs <- matrix(NA, nrow = n_boot, ncol = n_coef)
  colnames(boot_coefs) <- names(original_coefs)

  for (b in 1:n_boot) {
    # Block bootstrap indices
    boot_idx <- block_bootstrap_indices(n, block_size)
    boot_data <- data[boot_idx, , drop = FALSE]
    boot_weights <- if (is.null(weights)) NULL else weights[boot_idx]

    # Refit on bootstrap sample
    tryCatch({
      if (inherits(equation, "lm")) {
        args <- list(formula = eq_form, data = boot_data)
        if (!is.null(boot_weights)) args$weights <- boot_weights
        fit <- do.call(stats::lm, args)
        boot_coefs[b, ] <- stats::coef(fit)[names(original_coefs)]
      } else if (inherits(equation, "nls")) {
        if (!requireNamespace("minpack.lm", quietly = TRUE)) {
          stop("minpack.lm required")
        }
        args <- list(formula = eq_form, data = boot_data,
                     start = as.list(original_coefs),
                     control = minpack.lm::nls.lm.control(maxiter = 100))
        if (!is.null(boot_weights)) args$weights <- boot_weights
        fit <- do.call(minpack.lm::nlsLM, args)
        boot_coefs[b, ] <- stats::coef(fit)[names(original_coefs)]
      } else if (inherits(equation, "symbolic_equation")) {
        fit <- fit_specified_equation(eq_expr, data = boot_data,
                                      response = response,
                                      start = as.list(original_coefs),
                                      method = cv_fit_method(eq_expr, boot_data),
                                      weights = boot_weights)
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