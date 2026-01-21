#' @title Visual Exploration of Dynamical Structure
#' @description Functions for visually exploring the structure of dynamical
#'   systems before formal model fitting. These diagnostics should be used
#'   BEFORE any symbolic search to inform hypotheses about functional forms.
#' @name exploration
NULL


#' Comprehensive Dynamics Exploration
#'
#' Generates a battery of diagnostic plots to explore the dynamical structure
#' of the data and suggests potential functional forms.
#'
#' @param data Data frame containing the time series.
#' @param target Name of the target variable (or its derivative).
#' @param predictors Character vector of predictor variable names.
#' @param time Name of the time column (auto-detected if NULL).
#' @param n_bins Number of bins for conditional analysis.
#' @param include Which plots to include: "all", or subset of 
#'   c("timeseries", "phase", "bivariate", "interactions").
#'
#' @return A list containing:
#' \itemize{
#'   \item suggestions: Character vector of suggested functional forms
#'   \item statistics: Data frame of diagnostic statistics
#'   \item plots: List of ggplot objects (if available)
#' }
#'
#' @examples
#' \donttest{
#' # Toy example
#' data <- data.frame(
#'   time = 1:50,
#'   Z = sin(seq(0, 10, length.out = 50)),
#'   X = cos(seq(0, 10, length.out = 50))
#' )
#' data$dZ <- c(diff(data$Z)/diff(data$time), NA)
#' data <- na.omit(data)
#' 
#' result <- explore_dynamics(data, 
#'   target = "dZ",
#'   predictors = c("Z", "X")
#' )
#' print(result$suggestions)
#' }
#'
#' @export
explore_dynamics <- function(data, target, predictors = NULL, time = NULL,
                             n_bins = 10, include = "all") {
  
  # Validate inputs
  if (!target %in% names(data)) {
    stop("Target variable '", target, "' not found in data", call. = FALSE)
  }
  
  # Auto-detect predictors from variable specification
  if (is.null(predictors)) {
    var_spec <- attr(data, "var_spec")
    if (!is.null(var_spec)) {
      predictors <- var_spec$all_predictors
    } else {
      # Use all numeric columns except target and time
      num_cols <- names(data)[sapply(data, is.numeric)]
      predictors <- setdiff(num_cols, c(target, time))
    }
  }
  
  # Validate predictors exist
  missing <- setdiff(predictors, names(data))
  if (length(missing) > 0) {
    warning("Predictors not found: ", paste(missing, collapse = ", "))
    predictors <- intersect(predictors, names(data))
  }
  
  # Auto-detect time column
  if (is.null(time)) {
    var_spec <- attr(data, "var_spec")
    if (!is.null(var_spec) && !is.null(var_spec$time)) {
      time <- var_spec$time
    } else {
      time_candidates <- c("time", "t", "date", "year", "period")
      time <- intersect(time_candidates, names(data))[1]
    }
  }
  
  # Initialize results
  suggestions <- character()
  statistics <- list()
  plots <- list()
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("  DYNAMICS EXPLORATION: ", target)
  message(paste(rep("=", 60), collapse = ""), "\n")
  
  y <- data[[target]]
  
  # 1. Time series plot
  if (include == "all" || "timeseries" %in% include) {
    message("1. Time Series Analysis")
    message("   ---------------------")
    
    if (!is.null(time)) {
      plots$timeseries <- plot_timeseries(data, target, time)
    }
    
    # Basic statistics
    statistics$target_mean <- mean(y, na.rm = TRUE)
    statistics$target_sd <- stats::sd(y, na.rm = TRUE)
    statistics$target_range <- range(y, na.rm = TRUE)
    
    # Trend detection
    n <- length(y)
    t_idx <- 1:n
    trend_cor <- stats::cor(t_idx, y, use = "complete.obs")
    statistics$trend_correlation <- trend_cor
    
    if (!is.na(trend_cor) && abs(trend_cor) > 0.7) {
      message(sprintf("   Strong trend detected (r = %.3f)", trend_cor))
      suggestions <- c(suggestions, "Consider detrending or including time trend")
    } else if (!is.na(trend_cor) && abs(trend_cor) > 0.3) {
      message(sprintf("   Moderate trend detected (r = %.3f)", trend_cor))
    } else {
      message("   No significant trend detected")
    }
    
    cat("\n")
  }
  
  # 2. Phase diagram (dZ vs Z) for each predictor
  if (include == "all" || "phase" %in% include) {
    message("2. Phase Diagrams")
    message("   ---------------")
    
    for (pred in predictors) {
      if (pred == target) next
      
      x <- data[[pred]]
      
      # Correlation
      r <- stats::cor(x, y, use = "complete.obs")
      message(sprintf("   %s ~ %s: r = %.3f", target, pred, r))
      
      # Test for nonlinearity using RESET-like approach
      lm_lin <- stats::lm(y ~ x)
      x2 <- x^2
      x3 <- x^3
      lm_quad <- stats::lm(y ~ x + x2)
      lm_cubic <- stats::lm(y ~ x + x2 + x3)
      
      # Compare models
      aic_lin <- stats::AIC(lm_lin)
      aic_quad <- stats::AIC(lm_quad)
      aic_cubic <- stats::AIC(lm_cubic)
      
      best_form <- which.min(c(aic_lin, aic_quad, aic_cubic))
      form_names <- c("linear", "quadratic", "cubic")
      
      if (best_form > 1) {
        message(sprintf("      -> %s form suggested (AIC improvement: %.1f)",
                        form_names[best_form], aic_lin - min(aic_quad, aic_cubic)))
        
        if (best_form == 2) {
          suggestions <- c(suggestions, sprintf("%s^2 term for %s", pred, pred))
        } else {
          suggestions <- c(suggestions, sprintf("Polynomial in %s", pred))
        }
      }
      
      statistics[[paste0("cor_", pred)]] <- r
      statistics[[paste0("nonlin_", pred)]] <- form_names[best_form]
      
      # Generate phase plot
      plots[[paste0("phase_", pred)]] <- plot_phase_1d(data, pred, target)
    }
    cat("\n")
  }
  
  # 3. Bivariate relationships
  if (include == "all" || "bivariate" %in% include) {
    message("3. Bivariate Scatter Analysis")
    message("   ---------------------------")
    
    for (pred in predictors) {
      if (pred == target) next
      
      plots[[paste0("bivariate_", pred)]] <- plot_bivariate(data, pred, target)
      
      # Detect saturation/asymptotic behavior
      x <- data[[pred]]
      x_range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
      
      # Check for logistic-type saturation
      # Group into tertiles and check if slope decreases
      breaks <- stats::quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
      # Ensure breaks are unique
      if (length(unique(breaks)) < 4) {
        # Fallback if insufficient unique values
        groups <- as.numeric(cut(x, 3))
      } else {
        groups <- cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      }
      
      slopes <- numeric(3)
      for (g in 1:3) {
        idx <- which(groups == g)
        if (length(idx) > 5) {
          slopes[g] <- stats::coef(stats::lm(y[idx] ~ x[idx]))[2]
        }
      }
      
      if (!any(is.na(slopes))) {
        if (all(diff(slopes) < 0) && slopes[1] > 0) {
          message(sprintf("   %s: Possible saturation/logistic behavior", pred))
          suggestions <- c(suggestions, 
                           sprintf("Logistic saturation in %s", pred))
        } else if (all(diff(slopes) > 0) && slopes[3] > 0) {
          message(sprintf("   %s: Possible exponential growth", pred))
          suggestions <- c(suggestions,
                           sprintf("Exponential term in %s", pred))
        }
      }
    }
    cat("\n")
  }
  
  # 4. Interaction detection
  if (include == "all" || "interactions" %in% include) {
    if (length(predictors) >= 2) {
      message("4. Interaction Detection")
      message("   ----------------------")
      
      for (i in 1:(length(predictors) - 1)) {
        for (j in (i + 1):length(predictors)) {
          pred1 <- predictors[i]
          pred2 <- predictors[j]
          
          if (pred1 == target || pred2 == target) next
          
          x1 <- data[[pred1]]
          x2 <- data[[pred2]]
          
          # Test for interaction
          lm_add <- stats::lm(y ~ x1 + x2)
          lm_int <- stats::lm(y ~ x1 * x2)
          
          # F-test for interaction
          anova_result <- stats::anova(lm_add, lm_int)
          p_interaction <- anova_result$`Pr(>F)`[2]
          
          if (!is.na(p_interaction) && p_interaction < 0.05) {
            message(sprintf("   Significant interaction: %s * %s (p = %.4f)",
                            pred1, pred2, p_interaction))
            suggestions <- c(suggestions,
                             sprintf("Interaction term %s * %s", pred1, pred2))
            
            statistics[[paste0("interaction_", pred1, "_", pred2)]] <- p_interaction
          }
        }
      }
      cat("\n")
    }
  }
  
  # 5. Summary of suggestions
  message("5. Summary of Suggested Functional Forms")
  message("   --------------------------------------")
  
  if (length(suggestions) > 0) {
    suggestions <- unique(suggestions)
    for (s in suggestions) {
      message("   * ", s)
    }
  } else {
    message("   No strong nonlinear patterns detected")
    message("   Consider starting with linear specification")
  }
  
  message("\n", paste(rep("=", 60), collapse = ""))
  
  result <- list(
    suggestions = suggestions,
    statistics = statistics,
    plots = plots,
    target = target,
    predictors = predictors
  )
  
  class(result) <- "dynamics_exploration"
  result
}


#' Time Series Plot
#'
#' Creates a time series plot with optional trend line and change point detection.
#'
#' @param data Data frame.
#' @param var Variable name to plot.
#' @param time Time variable name.
#' @param show_trend Add trend line?
#' @param highlight_changes Highlight potential structural breaks?
#'
#' @return A ggplot object.
#'
#' @export
plot_timeseries <- function(data, var, time = NULL, show_trend = TRUE,
                            highlight_changes = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  if (is.null(time)) {
    data$.time <- seq_len(nrow(data))
    time <- ".time"
  }
  
  y <- data[[var]]
  t <- data[[time]]
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[time]], y = .data[[var]])) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::labs(
      title = paste("Time Series:", var),
      x = "Time", y = var
    ) +
    ed_theme()
  
  if (show_trend) {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, 
                                  color = "red", linetype = "dashed",
                                  linewidth = 0.5)
  }
  
  if (highlight_changes) {
    # Simple change detection using rolling variance
    n <- length(y)
    if (n > 20) {
      window <- max(5, floor(n / 10))
      rolling_var <- sapply(1:n, function(i) {
        start <- max(1, i - window)
        end <- min(n, i + window)
        stats::var(y[start:end], na.rm = TRUE)
      })
      
      # Find peaks in variance (potential change points)
      var_thresh <- mean(rolling_var, na.rm = TRUE) + 2 * stats::sd(rolling_var, na.rm = TRUE)
      change_points <- which(rolling_var > var_thresh & 
                               c(FALSE, diff(rolling_var) > 0) &
                               c(diff(rolling_var) < 0, FALSE))
      
      if (length(change_points) > 0 && length(change_points) <= 5) {
        change_data <- data.frame(x = t[change_points])
        p <- p + ggplot2::geom_vline(data = change_data,
                                     ggplot2::aes(xintercept = x),
                                     color = "orange", linetype = "dotted",
                                     linewidth = 0.8)
      }
    }
  }
  
  p
}


#' 1D Phase Diagram
#'
#' Creates a phase diagram plotting dZ vs Z, useful for visualizing autonomous
#' dynamics and identifying fixed points.
#'
#' @param data Data frame.
#' @param z_var Name of state variable Z.
#' @param dz_var Name of derivative dZ (auto-constructed if starts with "d_").
#' @param show_zero_line Add horizontal line at dZ = 0?
#' @param show_fit Add nonparametric fit?
#' @param fit_method Method for fit: "loess", "gam", or "spline".
#'
#' @return A ggplot object.
#'
#' @details
#' In the phase diagram:
#' \itemize{
#'   \item Points where the curve crosses dZ = 0 are fixed points
#'   \item Negative slope at crossing indicates stability
#'   \item Positive slope indicates instability
#' }
#'
#' @export
plot_phase_1d <- function(data, z_var, dz_var = NULL, show_zero_line = TRUE,
                          show_fit = TRUE, fit_method = "loess") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  # Auto-construct derivative name if needed
  if (is.null(dz_var)) {
    dz_var <- paste0("d_", z_var)
  }
  
  if (!dz_var %in% names(data)) {
    stop("Derivative variable '", dz_var, "' not found. ",
         "Use compute_derivatives() first.", call. = FALSE)
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[z_var]], y = .data[[dz_var]])) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
    ggplot2::labs(
      title = paste("Phase Diagram:", dz_var, "vs", z_var),
      subtitle = "Fixed points where curve crosses zero",
      x = z_var, y = expression(dot(Z))
    ) +
    ed_theme()
  
  if (show_zero_line) {
    p <- p + ggplot2::geom_hline(yintercept = 0, color = "red", 
                                 linetype = "dashed", linewidth = 0.5)
  }
  
  if (show_fit) {
    if (fit_method == "loess") {
      p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE,
                                    color = "darkred", fill = "pink",
                                    alpha = 0.3)
    } else if (fit_method == "gam") {
      # Use gam only if mgcv is available
      if (requireNamespace("mgcv", quietly = TRUE)) {
        p <- p + ggplot2::geom_smooth(method = "gam", 
                                      formula = y ~ s(x, bs = "cs"),
                                      se = TRUE, color = "darkred")
      } else {
        warning("mgcv package not installed. Using default smoothing method.")
        p <- p + ggplot2::geom_smooth(se = TRUE, color = "darkred")
      }
    } else {
      p <- p + ggplot2::geom_smooth(se = TRUE, color = "darkred")
    }
  }
  
  p
}


#' Bivariate Scatter Plot
#'
#' Creates a scatter plot with optional nonparametric fit and marginal
#' distributions.
#'
#' @param data Data frame.
#' @param x_var X variable name.
#' @param y_var Y variable name.
#' @param color_var Optional variable for color mapping.
#' @param show_fit Add smooth fit?
#' @param show_marginals Add marginal histograms?
#'
#' @return A ggplot object.
#'
#' @export
plot_bivariate <- function(data, x_var, y_var, color_var = NULL,
                           show_fit = TRUE, show_marginals = FALSE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  if (!is.null(color_var) && color_var %in% names(data)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], 
                                            y = .data[[y_var]],
                                            color = .data[[color_var]]))
  } else {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], 
                                            y = .data[[y_var]]))
  }
  
  p <- p +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::labs(
      title = paste(y_var, "vs", x_var),
      x = x_var, y = y_var
    ) +
    ed_theme()
  
  if (is.null(color_var)) {
    p <- p + ggplot2::aes(color = NULL) +
      ggplot2::geom_point(color = "steelblue", alpha = 0.6)
  }
  
  if (show_fit) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE,
                                  color = "darkred", fill = "pink",
                                  alpha = 0.3)
  }
  
  p
}


#' 3D Response Surface
#'
#' Creates a 3D surface or contour plot showing how the target variable
#' depends on two predictors.
#'
#' @param data Data frame.
#' @param x_var First predictor variable.
#' @param y_var Second predictor variable.
#' @param z_var Response variable (target).
#' @param type Plot type: "contour", "filled_contour", or "persp".
#' @param n_grid Grid resolution for surface estimation.
#' @param method Surface fitting method: "loess", "gam", or "linear".
#'
#' @return A plot (base graphics for persp, ggplot for contour).
#'
#' @export
plot_surface_3d <- function(data, x_var, y_var, z_var,
                            type = c("contour", "filled_contour", "persp"),
                            n_grid = 30, method = "loess") {
  
  type <- match.arg(type)
  
  x <- data[[x_var]]
  y <- data[[y_var]]
  z <- data[[z_var]]
  
  # Create grid
  x_grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n_grid)
  y_grid <- seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length.out = n_grid)
  grid <- expand.grid(x = x_grid, y = y_grid)
  names(grid) <- c(x_var, y_var)
  
  # Fit surface
  if (method == "loess") {
    fit <- stats::loess(stats::as.formula(paste(z_var, "~", x_var, "*", y_var)), 
                        data = data, span = 0.5)
    grid$z <- stats::predict(fit, newdata = grid)
  } else if (method == "linear") {
    fit <- stats::lm(stats::as.formula(paste(z_var, "~", x_var, "*", y_var)), data = data)
    grid$z <- stats::predict(fit, newdata = grid)
  } else {
    # Simple interpolation / default
    fit <- stats::loess(stats::as.formula(paste(z_var, "~", x_var, "*", y_var)), 
                        data = data)
    grid$z <- stats::predict(fit, newdata = grid)
  }
  
  z_mat <- matrix(grid$z, nrow = n_grid, ncol = n_grid)
  
  if (type == "persp") {
    graphics::persp(x_grid, y_grid, z_mat,
                    xlab = x_var, ylab = y_var, zlab = z_var,
                    main = paste("Response Surface:", z_var),
                    theta = 30, phi = 20,
                    col = "lightblue", shade = 0.5)
    
  } else if (type == "contour") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for contour plots.")
    }
    
    grid$z_clean <- grid$z
    grid$z_clean[is.na(grid$z_clean)] <- mean(grid$z, na.rm = TRUE)
    
    ggplot2::ggplot(grid, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], 
                                       z = z_clean)) +
      ggplot2::geom_contour_filled(bins = 15) +
      ggplot2::geom_point(data = data, 
                          ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]],
                                       z = NULL),
                          alpha = 0.3, size = 0.5) +
      ggplot2::labs(
        title = paste("Contour Plot:", z_var),
        x = x_var, y = y_var, fill = z_var
      ) +
      ed_theme() +
      ggplot2::theme(legend.position = "right")
    
  } else {
    # filled_contour using base graphics
    graphics::filled.contour(x_grid, y_grid, z_mat,
                             xlab = x_var, ylab = y_var,
                             main = paste("Response Surface:", z_var),
                             color.palette = grDevices::colorRampPalette(c("blue", "white", "red")))
  }
}


#' 2D Trajectory Plot
#'
#' Plots the trajectory of a system in the (Z, X) plane, useful for
#' visualizing attractors and limit cycles.
#'
#' @param data Data frame.
#' @param x_var First state variable.
#' @param y_var Second state variable.
#' @param time_var Time variable (for coloring trajectory).
#' @param show_arrows Add direction arrows?
#' @param arrow_spacing Spacing between arrows (every nth point).
#' @param show_start_end Mark start and end points?
#'
#' @return A ggplot object.
#'
#' @export
plot_trajectory_2d <- function(data, x_var, y_var, time_var = NULL,
                               show_arrows = TRUE, arrow_spacing = 10,
                               show_start_end = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  x <- data[[x_var]]
  y <- data[[y_var]]
  n <- length(x)
  
  # Create data for plotting
  plot_data <- data.frame(x = x, y = y, idx = 1:n)
  
  if (!is.null(time_var) && time_var %in% names(data)) {
    plot_data$time <- data[[time_var]]
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = time))
  } else {
    plot_data$time <- 1:n
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = idx))
  }
  
  p <- p +
    ggplot2::geom_path(linewidth = 0.8, alpha = 0.8) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::labs(
      title = paste("Phase Trajectory:", y_var, "vs", x_var),
      x = x_var, y = y_var, color = "Time"
    ) +
    ed_theme()
  
  if (show_arrows) {
    arrow_idx <- seq(1, n - 1, by = arrow_spacing)
    # Ensure arrows don't exceed data length
    arrow_idx <- arrow_idx[arrow_idx < n]
    
    if (length(arrow_idx) > 0) {
      arrow_data <- data.frame(
        x = x[arrow_idx],
        y = y[arrow_idx],
        xend = x[arrow_idx + 1],
        yend = y[arrow_idx + 1]
      )
      
      p <- p +
        ggplot2::geom_segment(
          data = arrow_data,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
          color = "gray40", alpha = 0.6, inherit.aes = FALSE
        )
    }
  }
  
  if (show_start_end) {
    endpoints <- data.frame(
      x = c(x[1], x[n]),
      y = c(y[1], y[n]),
      label = c("Start", "End")
    )
    
    p <- p +
      ggplot2::geom_point(
        data = endpoints,
        ggplot2::aes(x = x, y = y, shape = label),
        size = 4, inherit.aes = FALSE
      ) +
      ggplot2::scale_shape_manual(values = c("Start" = 17, "End" = 15))
  }
  
  p
}


#' Annotate Hypotheses
#'
#' Records researcher hypotheses based on visual exploration, to be used
#' as constraints or guides in subsequent symbolic search.
#'
#' @param data Data frame (with exploration results as attribute).
#' @param hypotheses Character vector of hypotheses.
#'
#' @return Data frame with hypotheses attached as attribute.
#'
#' @examples
#' \donttest{
#' # Toy example
#' data <- data.frame(Z = 1:10)
#' data <- annotate_hypotheses(data, c(
#'   "Z exhibits logistic saturation around Z=100",
#'   "Effect of X appears linear"
#' ))
#' }
#'
#' @export
annotate_hypotheses <- function(data, hypotheses) {
  
  existing <- attr(data, "hypotheses")
  
  if (!is.null(existing)) {
    hypotheses <- c(existing, hypotheses)
  }
  
  attr(data, "hypotheses") <- hypotheses
  
  message("Recorded hypotheses:")
  for (i in seq_along(hypotheses)) {
    message(sprintf("  %d. %s", i, hypotheses[i]))
  }
  
  data
}