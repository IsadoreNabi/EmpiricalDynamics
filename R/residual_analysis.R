#' @title Residual Analysis and Stochastic Differential Equations
#' @description Functions for analyzing residual structure, modeling conditional
#'   variance, constructing stochastic differential equations (SDEs), and
#'   iterative GLS estimation for heteroscedastic systems.
#' @name residual_analysis
NULL


#' Compute Residuals from Symbolic Equation
#'
#' Calculates residuals from a fitted symbolic equation or SDE model.
#'
#' @param model A symbolic_equation or sde_model object
#' @param data Data frame containing the variables
#' @param target Name of target variable (auto-detected if NULL)
#'
#' @return Numeric vector of residuals
#'
#' @export
compute_residuals <- function(model, data, target = NULL) {
  
  # Get target variable
  if (is.null(target)) {
    if (inherits(model, "sde_model")) {
      target <- model$drift$target
    } else if (inherits(model, "symbolic_equation")) {
      target <- model$target
    }
  }
  
  if (is.null(target)) {
    # Try to find derivative column
    d_cols <- grep("^d_|^d[A-Z]", names(data), value = TRUE)
    if (length(d_cols) == 1) {
      target <- d_cols[1]
    } else {
      stop("Cannot auto-detect target variable. Please specify 'target'.", 
           call. = FALSE)
    }
  }
  
  y <- data[[target]]
  
  # Get predictions
  if (inherits(model, "sde_model")) {
    pred <- predict(model$drift, newdata = data)
  } else if (inherits(model, "symbolic_equation")) {
    pred <- predict(model, newdata = data)
  } else if (inherits(model, "nls") || inherits(model, "lm")) {
    pred <- predict(model, newdata = data)
  } else {
    stop("Unsupported model type", call. = FALSE)
  }
  
  if (is.null(pred)) {
    stop("Could not compute predictions from model", call. = FALSE)
  }
  
  y - pred
}


#' Comprehensive Residual Diagnostics
#'
#' Performs a battery of statistical tests on model residuals to check
#' for autocorrelation, heteroscedasticity, and normality.
#'
#' @param residuals Numeric vector of residuals (or model object)
#' @param data Optional data frame for conditional tests
#' @param predictors Variable names for heteroscedasticity tests
#' @param max_lag Maximum lag for autocorrelation tests
#' @param plot Produce diagnostic plots?
#'
#' @return A list of test results with class "residual_diagnostics"
#'
#' @export
residual_diagnostics <- function(residuals, data = NULL, predictors = NULL,
                                 max_lag = 10, plot = TRUE) {
  
  # Extract residuals if model object passed
  if (inherits(residuals, c("symbolic_equation", "sde_model", "lm", "nls"))) {
    if (is.null(data)) {
      stop("'data' required when passing model object", call. = FALSE)
    }
    residuals <- compute_residuals(residuals, data)
  }
  
  n <- length(residuals)
  results <- list()
  
  # 1. Basic statistics
  results$basic <- list(
    n = n,
    mean = mean(residuals, na.rm = TRUE),
    sd = sd(residuals, na.rm = TRUE),
    min = min(residuals, na.rm = TRUE),
    max = max(residuals, na.rm = TRUE),
    skewness = compute_skewness(residuals),
    kurtosis = compute_kurtosis(residuals)
  )
  
  # 2. Ljung-Box test for autocorrelation
  results$ljung_box <- tryCatch({
    if (n >= 10) {
      Box.test(residuals, lag = min(max_lag, floor(n/5)), type = "Ljung-Box")
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))
  
  # 3. ARCH-LM test for conditional heteroscedasticity
  results$arch_lm <- tryCatch({
    if (n >= 20) {
      arch_test(residuals, lags = min(5, floor(n/10)))
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))
  
  # 4. Breusch-Pagan test for heteroscedasticity vs predictors
  if (!is.null(data) && !is.null(predictors)) {
    available_preds <- intersect(predictors, names(data))
    if (length(available_preds) > 0) {
      bp_data <- cbind(data[available_preds], residuals = residuals)
      
      results$breusch_pagan <- tryCatch({
        lm_fit <- lm(residuals ~ ., data = bp_data)
        lmtest::bptest(lm_fit)
      }, error = function(e) list(statistic = NA, p.value = NA))
    }
  }
  
  # 5. Jarque-Bera test for normality
  results$jarque_bera <- tryCatch({
    if (n >= 10) {
      tseries::jarque.bera.test(residuals)
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))
  
  # 6. Runs test for randomness
  results$runs <- tryCatch({
    if (n >= 10) {
      runs_test(residuals)
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))
  
  # Build tests data frame
  results$tests <- data.frame(
    test = c("Ljung-Box", "ARCH-LM", "Jarque-Bera", "Runs"),
    statistic = c(
      as.numeric(results$ljung_box$statistic),
      as.numeric(results$arch_lm$statistic),
      as.numeric(results$jarque_bera$statistic),
      as.numeric(results$runs$statistic)
    ),
    p_value = c(
      results$ljung_box$p.value,
      results$arch_lm$p.value,
      results$jarque_bera$p.value,
      results$runs$p.value
    ),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(results$breusch_pagan)) {
    results$tests <- rbind(results$tests, data.frame(
      test = "Breusch-Pagan",
      statistic = as.numeric(results$breusch_pagan$statistic),
      p_value = results$breusch_pagan$p.value,
      stringsAsFactors = FALSE
    ))
  }
  
  # Determine if issues detected
  results$issues <- character()
  
  if (!is.na(results$ljung_box$p.value) && results$ljung_box$p.value < 0.05) {
    results$issues <- c(results$issues, "Autocorrelation detected (Ljung-Box)")
  }
  
  if (!is.na(results$arch_lm$p.value) && results$arch_lm$p.value < 0.05) {
    results$issues <- c(results$issues, "ARCH effects detected (conditional heteroscedasticity)")
  }
  
  if (!is.null(results$breusch_pagan) && 
      !is.na(results$breusch_pagan$p.value) && 
      results$breusch_pagan$p.value < 0.05) {
    results$issues <- c(results$issues, "Heteroscedasticity detected (Breusch-Pagan)")
  }
  
  if (!is.na(results$jarque_bera$p.value) && results$jarque_bera$p.value < 0.05) {
    results$issues <- c(results$issues, "Non-normality detected (Jarque-Bera)")
  }
  
  results$residuals <- residuals
  
  class(results) <- "residual_diagnostics"
  
  if (plot && n >= 5) {
    tryCatch({
      plot_residual_diagnostics_panel(results)
    }, error = function(e) {
      warning("Could not produce diagnostic plots: ", e$message)
    })
  }
  
  results
}


#' ARCH-LM Test
#'
#' @keywords internal
arch_test <- function(x, lags = 5) {
  n <- length(x)
  x2 <- x^2
  
  # Regress squared residuals on lagged squared residuals
  y <- x2[(lags + 1):n]
  X <- sapply(1:lags, function(i) x2[(lags + 1 - i):(n - i)])
  
  fit <- lm(y ~ X)
  r2 <- summary(fit)$r.squared
  
  test_stat <- n * r2
  p_value <- 1 - pchisq(test_stat, df = lags)
  
  list(
    statistic = test_stat,
    p.value = p_value,
    df = lags,
    method = "ARCH-LM Test"
  )
}


#' Runs Test for Randomness
#'
#' @keywords internal
runs_test <- function(x) {
  n <- length(x)
  m <- median(x)
  signs <- x > m
  
  # Count runs
  runs <- 1 + sum(diff(as.numeric(signs)) != 0)
  n1 <- sum(signs)
  n2 <- n - n1
  
  # Expected runs and variance
  expected <- 1 + 2 * n1 * n2 / n
  variance <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
  
  if (variance <= 0) {
    return(list(statistic = NA, p.value = NA, runs = runs, expected = expected, method = "Wald-Wolfowitz Runs Test"))
  }
  
  z <- (runs - expected) / sqrt(variance)
  p_value <- 2 * pnorm(-abs(z))
  
  list(
    statistic = z,
    p.value = p_value,
    runs = runs,
    expected = expected,
    method = "Wald-Wolfowitz Runs Test"
  )
}


#' Compute Skewness
#'
#' @keywords internal
compute_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA)
  sum((x - m)^3) / (n * s^3)
}


#' Compute Excess Kurtosis
#'
#' @keywords internal
compute_kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(NA)
  sum((x - m)^4) / (n * s^4) - 3
}


#' Print Residual Diagnostics
#'
#' @param x Object of class residual_diagnostics
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object (called for side effects).
#' @export
print.residual_diagnostics <- function(x, ...) {
  cat("\n=== Residual Diagnostics ===\n\n")
  
  cat("Basic Statistics:\n")
  cat(sprintf("  N: %d\n", x$basic$n))
  cat(sprintf("  Mean: %.6f (should be ~0)\n", x$basic$mean))
  cat(sprintf("  SD: %.6f\n", x$basic$sd))
  cat(sprintf("  Skewness: %.3f\n", x$basic$skewness))
  cat(sprintf("  Excess Kurtosis: %.3f\n", x$basic$kurtosis))
  
  cat("\nStatistical Tests:\n")
  
  cat(sprintf("  Ljung-Box (autocorrelation): stat=%.2f, p=%.4f %s\n",
              x$ljung_box$statistic, x$ljung_box$p.value,
              if (!is.na(x$ljung_box$p.value) && x$ljung_box$p.value < 0.05) "***" else ""))
  
  cat(sprintf("  ARCH-LM (het. conditional): stat=%.2f, p=%.4f %s\n",
              x$arch_lm$statistic, x$arch_lm$p.value,
              if (!is.na(x$arch_lm$p.value) && x$arch_lm$p.value < 0.05) "***" else ""))
  
  if (!is.null(x$breusch_pagan)) {
    cat(sprintf("  Breusch-Pagan (het.): stat=%.2f, p=%.4f %s\n",
                x$breusch_pagan$statistic, x$breusch_pagan$p.value,
                if (!is.na(x$breusch_pagan$p.value) && x$breusch_pagan$p.value < 0.05) "***" else ""))
  }
  
  cat(sprintf("  Jarque-Bera (normality): stat=%.2f, p=%.4f %s\n",
              x$jarque_bera$statistic, x$jarque_bera$p.value,
              if (!is.na(x$jarque_bera$p.value) && x$jarque_bera$p.value < 0.05) "***" else ""))
  
  if (length(x$issues) > 0) {
    cat("\nIssues Detected:\n")
    for (issue in x$issues) {
      cat("  * ", issue, "\n")
    }
  } else {
    cat("\nNo significant issues detected.\n")
  }
  
  invisible(x)
}


#' Plot Residual Diagnostics Panel
#'
#' Creates a multi-panel diagnostic plot for residual analysis.
#'
#' @param x Object of class residual_diagnostics
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisibly returns the input object
#'
#' @export
plot_residual_diagnostics_panel <- function(x, ...) {
  if (!inherits(x, "residual_diagnostics")) {
    stop("x must be of class 'residual_diagnostics'")
  }
  
  residuals <- x$residuals
  n <- length(residuals)
  
  if (n < 5) {
    warning("Too few observations for diagnostic plots")
    return(invisible(x))
  }
  
  # Set up 2x2 panel
  oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(oldpar))
  
  # 1. Residuals vs index
  plot(seq_len(n), residuals, type = "p", pch = 16, col = "steelblue",
       xlab = "Index", ylab = "Residuals", main = "Residuals vs Index")
  abline(h = 0, lty = 2, col = "red")
  
  # 2. Histogram
  hist(residuals, breaks = "Sturges", col = "lightblue", border = "white",
       main = "Histogram of Residuals", xlab = "Residuals", freq = FALSE)
  if (n >= 5) {
    curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), 
          add = TRUE, col = "red", lwd = 2)
  }
  
  # 3. Q-Q plot
  qqnorm(residuals, pch = 16, col = "steelblue", main = "Normal Q-Q Plot")
  qqline(residuals, col = "red", lwd = 2)
  
  # 4. ACF
  if (n >= 10) {
    acf(residuals, main = "Autocorrelation", lag.max = min(20, n - 1))
  } else {
    plot.new()
    text(0.5, 0.5, "ACF requires >= 10 observations", cex = 1.2)
  }
  
  invisible(x)
}


#' Model Conditional Variance
#'
#' Estimates how the residual variance depends on state variables, used
#' for constructing the diffusion term of an SDE.
#'
#' @param residuals Numeric vector of residuals
#' @param predictors Formula or data frame of predictor variables (or vector of names)
#' @param data Data frame (if predictors is a formula or vector of names)
#' @param method Modeling method: "symbolic", "linear", "quadratic", "gam", or "constant"
#' @param transform Transformation of residuals: "squared", "absolute", or "log_squared"
#' @param ... Additional arguments passed to the modeling function
#'
#' @return An object of class "variance_model" containing the fitted model
#'
#' @export
model_conditional_variance <- function(residuals, predictors, data = NULL,
                                       method = c("symbolic", "linear", 
                                                  "quadratic", "gam", "constant"),
                                       transform = c("absolute", "squared", 
                                                     "log_squared"),
                                       ...) {
  
  method <- match.arg(method)
  transform <- match.arg(transform)
  
  # Prepare predictor data
  pred_vars <- NULL
  pred_data <- NULL
  
  if (method != "constant") {
    if (inherits(predictors, "formula")) {
      if (is.null(data)) {
        stop("'data' required when 'predictors' is a formula", call. = FALSE)
      }
      pred_vars <- all.vars(predictors)
      pred_data <- data[, pred_vars, drop = FALSE]
    } else if (is.character(predictors)) {
      if (is.null(data)) {
        stop("'data' required when 'predictors' is a vector of names", call. = FALSE)
      }
      pred_vars <- predictors
      pred_data <- data[, pred_vars, drop = FALSE]
    } else {
      pred_data <- as.data.frame(predictors)
      pred_vars <- names(pred_data)
    }
  }
  
  # Transform residuals
  target <- switch(transform,
                   "squared" = residuals^2,
                   "absolute" = abs(residuals),
                   "log_squared" = log(residuals^2 + 1e-10)
  )
  
  # Fit model
  fit <- switch(method,
                "constant" = {
                  # Fit intercept-only model
                  lm_data <- data.frame(.variance = target)
                  lm(.variance ~ 1, data = lm_data)
                },
                "linear" = {
                  lm_data <- cbind(pred_data, .variance = target)
                  lm(.variance ~ ., data = lm_data)
                },
                "quadratic" = {
                  # Add squared terms
                  for (v in pred_vars) {
                    pred_data[[paste0(v, "_sq")]] <- pred_data[[v]]^2
                  }
                  lm_data <- cbind(pred_data, .variance = target)
                  lm(.variance ~ ., data = lm_data)
                },
                "symbolic" = {
                  # Use symbolic search
                  search_result <- symbolic_search(
                    target = target,
                    predictors = pred_data,
                    n_runs = 3,
                    complexity_penalty = 0.1,
                    verbose = FALSE,
                    ...
                  )
                  select_equation(search_result, criterion = "knee")
                },
                "gam" = {
                  # Simple smooth model using loess (requires no extra packages for basic case)
                  if (ncol(pred_data) == 1) {
                    loess(target ~ pred_data[[1]], span = 0.5)
                  } else if (ncol(pred_data) == 2) {
                    # Simple interaction
                    formula_str <- paste("target ~", pred_vars[1], "*", pred_vars[2])
                    loess(as.formula(formula_str), 
                          data = cbind(pred_data, target = target), span = 0.5)
                  } else {
                    # Fall back to linear for >2 dims in simple gam
                    lm_data <- cbind(pred_data, .variance = target)
                    lm(.variance ~ ., data = lm_data)
                  }
                }
  )
  
  result <- list(
    fit = fit,
    method = method,
    transform = transform,
    predictors = pred_vars,
    target_transform = transform
  )
  
  class(result) <- "variance_model"
  result
}


#' Predict from Variance Model
#'
#' @param object Variance model object
#' @param newdata New data for prediction
#' @param ... Additional arguments
#'
#' @return Numeric vector of predicted standard deviations.
#' @export
predict.variance_model <- function(object, newdata, ...) {
  
  if (object$method == "constant") {
    # Return constant value based on intercept
    # Usually intercept is coef(fit)[1]
    val <- coef(object$fit)[1]
    pred <- rep(val, if(is.null(newdata)) 1 else nrow(newdata))
  } else {
    pred <- if (inherits(object$fit, "symbolic_equation")) {
      predict(object$fit, newdata = newdata)
    } else {
      predict(object$fit, newdata = newdata, ...)
    }
  }
  
  # Transform back to standard deviation scale
  if (object$transform == "squared") {
    sqrt(pmax(pred, 0))
  } else if (object$transform == "log_squared") {
    sqrt(exp(pred))
  } else {
    pmax(pred, 0)
  }
}


#' Construct Stochastic Differential Equation Model
#'
#' Combines a drift equation and diffusion model into a complete SDE:
#' dZ = f(Z, X, Y) dt + g(Z, X, Y) dW
#'
#' @param drift Symbolic equation for the drift term f(.)
#' @param diffusion Variance model for the diffusion term g(.)
#' @param variable Name of the main state variable
#' @param refine_with_gls Use iterative GLS to refine estimates?
#' @param gls_max_iter Maximum iterations for GLS
#' @param gls_tolerance Convergence tolerance for GLS
#' @param data Data frame (required if refine_with_gls = TRUE)
#' @param target Target variable name (required if refine_with_gls = TRUE)
#'
#' @return An object of class "sde_model"
#'
#' @export
construct_sde <- function(drift, diffusion = NULL,
                          variable = NULL,
                          refine_with_gls = FALSE,
                          gls_max_iter = 10,
                          gls_tolerance = 1e-4,
                          data = NULL,
                          target = NULL) {
  
  # If GLS refinement requested, do it first
  if (refine_with_gls) {
    if (is.null(data) || is.null(target)) {
      stop("'data' and 'target' required for GLS refinement", call. = FALSE)
    }
    
    # Get predictor names
    if (inherits(drift, "symbolic_equation")) {
      predictors <- drift$predictors
    } else {
      # Try to infer from data
      d_cols <- grep("^d_|^d[A-Z]", names(data), value = TRUE)
      predictors <- setdiff(names(data)[sapply(data, is.numeric)], 
                            c(d_cols, target))
    }
    
    pred_data <- data[, predictors, drop = FALSE]
    
    refined <- estimate_sde_iterative(
      target = data[[target]],
      predictors = pred_data,
      data = data,
      initial_drift = drift,
      max_iter = gls_max_iter,
      tol = gls_tolerance
    )
    
    return(refined)
  }
  
  # Simple construction without refinement
  sde <- list(
    drift = drift,
    diffusion = diffusion,
    variable = variable,
    estimation_method = "direct",
    refined_gls = FALSE
  )
  
  class(sde) <- "sde_model"
  sde
}


#' Iterative GLS Estimation for SDEs
#'
#' Refines drift and diffusion estimates using iterative Generalized Least
#' Squares, which is more appropriate when heteroscedasticity is substantial.
#'
#' @param target Numeric vector of target values (derivatives)
#' @param predictors Data frame of predictor variables
#' @param data Full data frame
#' @param initial_drift Initial drift equation (optional)
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance (RMSE change in coefficients)
#'
#' @return An sde_model object with refined estimates
#'
#' @export
estimate_sde_iterative <- function(target, predictors, data,
                                   initial_drift = NULL,
                                   max_iter = 10, tol = 1e-4) {
  
  n <- length(target)
  
  # Step 0: Initial drift estimate
  if (is.null(initial_drift)) {
    message("Step 0: Estimating initial drift (OLS)...")
    f_current <- symbolic_search(
      target = target,
      predictors = predictors,
      n_runs = 3,
      verbose = FALSE
    )
    f_current <- select_equation(f_current, criterion = "knee")
  } else {
    f_current <- initial_drift
  }
  
  # Initial residuals
  eps_current <- target - predict(f_current, newdata = data)
  if (any(is.na(eps_current))) {
    eps_current[is.na(eps_current)] <- 0
  }
  
  f_prev <- NULL
  g_current <- NULL
  
  for (i in 1:max_iter) {
    message(sprintf("\n--- Iteration %d ---", i))
    
    # Step i.a: Estimate diffusion
    message("  Estimating diffusion g...")
    g_current <- model_conditional_variance(
      residuals = eps_current,
      predictors = predictors,
      method = "linear",
      transform = "absolute"
    )
    
    # Step i.b: Calculate GLS weights
    variance_pred <- predict(g_current, newdata = data)^2
    variance_pred[variance_pred < 1e-10] <- 1e-10
    weights <- 1 / variance_pred
    weights <- weights / mean(weights)
    
    message(sprintf("  Weight range: [%.2f, %.2f]", min(weights), max(weights)))
    
    # Step i.c: Re-estimate drift with WLS
    message("  Re-estimating drift f (WLS)...")
    f_current <- symbolic_search_weighted(
      target = target,
      predictors = predictors,
      weights = weights,
      n_runs = 2,
      verbose = FALSE
    )
    f_current <- select_equation(f_current, criterion = "knee")
    
    # Step i.d: New residuals
    eps_current <- target - predict(f_current, newdata = data)
    if (any(is.na(eps_current))) {
      eps_current[is.na(eps_current)] <- 0
    }
    
    # Check convergence
    if (!is.null(f_prev)) {
      change <- coefficient_change(f_current, f_prev)
      message(sprintf("  Coefficient change: %.6f", change))
      
      if (change < tol) {
        message(sprintf("\nConvergence reached at iteration %d!", i))
        break
      }
    }
    
    f_prev <- f_current
    
    if (i == max_iter) {
      warning("Maximum iterations reached without convergence.")
    }
  }
  
  # Construct final SDE
  sde <- list(
    drift = f_current,
    diffusion = g_current,
    estimation_method = "iterative_gls",
    n_iterations = i,
    converged = (i < max_iter),
    final_weights = weights,
    final_residuals = eps_current
  )
  
  class(sde) <- "sde_model"
  sde
}


#' Compare OLS and GLS Estimation
#'
#' Produces a comparison of drift estimates from ordinary least squares
#' versus iterative GLS.
#'
#' @param ols_model SDE model estimated with OLS
#' @param gls_model SDE model estimated with iterative GLS
#' @param data Data frame used for estimation
#'
#' @return Invisibly returns comparison statistics
#'
#' @export
compare_estimation_methods <- function(ols_model, gls_model, data) {
  
  cat("\n=== Comparison: OLS vs GLS Iterative ===\n\n")
  
  # Coefficients
  cat("Drift Coefficients:\n")
  ols_coef <- coef(ols_model$drift)
  gls_coef <- coef(gls_model$drift)
  
  if (!is.null(ols_coef) && !is.null(gls_coef)) {
    cat("  OLS:", paste(names(ols_coef), "=", 
                        round(ols_coef, 4), collapse = ", "), "\n")
    cat("  GLS:", paste(names(gls_coef), "=", 
                        round(gls_coef, 4), collapse = ", "), "\n\n")
  } else {
    cat("  OLS equation:", ols_model$drift$string, "\n")
    cat("  GLS equation:", gls_model$drift$string, "\n\n")
  }
  
  # Residuals
  eps_ols <- compute_residuals(ols_model, data)
  eps_gls <- compute_residuals(gls_model, data)
  
  cat("Residual RMSE:\n")
  cat(sprintf("  OLS: %.4f\n", sqrt(mean(eps_ols^2, na.rm = TRUE))))
  cat(sprintf("  GLS: %.4f\n\n", sqrt(mean(eps_gls^2, na.rm = TRUE))))
  
  # Heteroscedasticity test
  cat("Breusch-Pagan Test (H0: homoscedasticity):\n")
  
  # Get a predictor for the test
  if (!is.null(ols_model$drift$predictors) && length(ols_model$drift$predictors) > 0) {
    pred_var <- ols_model$drift$predictors[1]
    if (pred_var %in% names(data)) {
      bp_ols <- tryCatch({
        lmtest::bptest(lm(eps_ols ~ data[[pred_var]]))$p.value
      }, error = function(e) NA)
      
      bp_gls <- tryCatch({
        lmtest::bptest(lm(eps_gls ~ data[[pred_var]]))$p.value
      }, error = function(e) NA)
      
      cat(sprintf("  OLS p-value: %.4f\n", bp_ols))
      cat(sprintf("  GLS p-value: %.4f\n", bp_gls))
    }
  }
  
  invisible(list(
    coef_change = if (!is.null(gls_coef) && !is.null(ols_coef)) {
      gls_coef[names(ols_coef)] - ols_coef
    } else NULL,
    rmse_improvement = sqrt(mean(eps_ols^2, na.rm = TRUE)) - 
      sqrt(mean(eps_gls^2, na.rm = TRUE))
  ))
}


#' Fit Residual Distribution
#'
#' Fits candidate probability distributions to residuals, optionally with
#' parameters that depend on state variables.
#'
#' @param residuals Numeric vector of residuals
#' @param candidates Character vector of distribution families to try
#' @param conditional_on Formula for conditional parameters (optional)
#' @param data Data frame (required if conditional_on specified)
#'
#' @return List with best fitting distribution and parameters
#'
#' @export
fit_residual_distribution <- function(residuals,
                                      candidates = c("normal", "t", 
                                                     "skew-normal"),
                                      conditional_on = NULL,
                                      data = NULL) {
  
  results <- list()
  
  for (dist in candidates) {
    fit <- tryCatch({
      switch(dist,
             "normal" = {
               list(
                 distribution = "normal",
                 params = c(mean = mean(residuals), sd = sd(residuals)),
                 loglik = sum(dnorm(residuals, mean(residuals), sd(residuals), log = TRUE)),
                 n_params = 2
               )
             },
             "t" = {
               # Fit Student's t via maximum likelihood
               fit_t <- fit_t_distribution(residuals)
               fit_t$distribution <- "t"
               fit_t
             },
             "skew-normal" = {
               # Simplified: use normal if skewness is small
               skew <- compute_skewness(residuals)
               if (abs(skew) < 0.5) {
                 list(
                   distribution = "skew-normal",
                   params = c(mean = mean(residuals), sd = sd(residuals), alpha = 0),
                   loglik = sum(dnorm(residuals, mean(residuals), sd(residuals), log = TRUE)),
                   n_params = 3
                 )
               } else {
                 list(
                   distribution = "skew-normal",
                   params = c(mean = mean(residuals), sd = sd(residuals), 
                              alpha = sign(skew) * min(abs(skew), 2)),
                   loglik = NA,
                   n_params = 3
                 )
               }
             }
      )
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      results[[dist]] <- fit
    }
  }
  
  # Compare using AIC
  n <- length(residuals)
  for (dist in names(results)) {
    if (!is.na(results[[dist]]$loglik)) {
      results[[dist]]$aic <- -2 * results[[dist]]$loglik + 2 * results[[dist]]$n_params
    } else {
      results[[dist]]$aic <- Inf
    }
  }
  
  # Find best
  aics <- sapply(results, function(x) x$aic)
  best_dist <- names(which.min(aics))
  
  list(
    best = results[[best_dist]],
    all = results,
    comparison = data.frame(
      distribution = names(results),
      aic = aics,
      stringsAsFactors = FALSE
    )
  )
}


#' Fit Student's t Distribution
#'
#' @keywords internal
fit_t_distribution <- function(x) {
  n <- length(x)
  
  # Initial estimates
  mu <- mean(x)
  sigma <- sd(x)
  nu <- 10  # Initial degrees of freedom
  
  # Simple grid search for nu
  best_ll <- -Inf
  best_nu <- nu
  
  for (test_nu in c(2, 3, 4, 5, 7, 10, 15, 20, 30, 50)) {
    ll <- sum(dt((x - mu) / sigma, df = test_nu, log = TRUE)) - n * log(sigma)
    if (ll > best_ll) {
      best_ll <- ll
      best_nu <- test_nu
    }
  }
  
  list(
    params = c(mu = mu, sigma = sigma, df = best_nu),
    loglik = best_ll,
    n_params = 3
  )
}


# S3 Methods for sde_model

#' @export
print.sde_model <- function(x, ...) {
  cat("\nStochastic Differential Equation Model\n")
  cat("=======================================\n\n")
  
  cat("Drift: dZ = f(Z, X, ...) dt\n")
  if (!is.null(x$drift$string)) {
    cat("  ", x$drift$string, "\n\n")
  } else if (!is.null(x$drift$expression)) {
    cat("  ", x$drift$expression, "\n\n")
  }
  
  if (!is.null(x$diffusion)) {
    cat("Diffusion: + g(Z, X, ...) dW\n")
    cat("  Method:", x$diffusion$method, "\n")
    cat("  Transform:", x$diffusion$transform, "\n\n")
  }
  
  if (!is.null(x$variable)) {
    cat("Variable:", x$variable, "\n")
  }
  
  cat("Estimation:", x$estimation_method, "\n")
  if (x$estimation_method == "iterative_gls") {
    cat("  Iterations:", x$n_iterations, "\n")
    cat("  Converged:", x$converged, "\n")
  }
  
  invisible(x)
}


#' @export
summary.sde_model <- function(object, ...) {
  cat("\n=== SDE Model Summary ===\n\n")
  
  cat("DRIFT COMPONENT\n")
  cat("---------------\n")
  summary(object$drift)
  
  if (!is.null(object$diffusion)) {
    cat("\nDIFFUSION COMPONENT\n")
    cat("-------------------\n")
    cat("Method:", object$diffusion$method, "\n")
    cat("Predictors:", paste(object$diffusion$predictors, collapse = ", "), "\n")
    
    if (inherits(object$diffusion$fit, "lm")) {
      cat("\nCoefficients:\n")
      print(coef(object$diffusion$fit))
    }
  }
  
  invisible(object)
}


#' @export
predict.sde_model <- function(object, newdata = NULL, 
                              component = c("drift", "diffusion", "both"),
                              ...) {
  component <- match.arg(component)
  
  result <- list()
  
  if (component %in% c("drift", "both")) {
    result$drift <- predict(object$drift, newdata = newdata)
  }
  
  if (component %in% c("diffusion", "both") && !is.null(object$diffusion)) {
    result$diffusion <- predict(object$diffusion, newdata = newdata)
  }
  
  if (component == "drift") {
    result$drift
  } else if (component == "diffusion") {
    result$diffusion
  } else {
    result
  }
}