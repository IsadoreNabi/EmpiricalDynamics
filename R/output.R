# =============================================================================
# EmpiricalDynamics: output.R
# Step F: Report Generation - LaTeX export, publication-ready outputs
# =============================================================================

#' @title Output and Report Generation
#' @description Functions for generating publication-ready outputs including
#'    LaTeX equations, comprehensive reports, and formatted summaries.
#' @name output
NULL

# =============================================================================
# LATEX EQUATION FORMATTING
# =============================================================================

#' Convert Equation to LaTeX
#'
#' Converts a discovered equation to LaTeX format for publication.
#'
#' @param equation Fitted equation object.
#' @param variable Name of the dependent variable (for dZ/dt notation).
#' @param precision Number of decimal places for coefficients.
#' @param scientific_notation Use scientific notation for large/small coefficients.
#' @param include_uncertainty Include standard errors in parentheses.
#' @param se_values Named vector of standard errors (optional).
#'
#' @return Character string with LaTeX equation.
#' @export
#' @examples
#' # Toy example using a linear model
#' data <- data.frame(Z = 1:10, dZ = 2 * (1:10) + 3)
#' model <- stats::lm(dZ ~ Z, data = data)
#' 
#' # Convert to LaTeX
#' latex_eq <- to_latex(model, variable = "Z")
#' cat(latex_eq)
to_latex <- function(equation, variable = "Z",
                     precision = 3,
                     scientific_notation = TRUE,
                     include_uncertainty = FALSE,
                     se_values = NULL) {
  
  # Extract expression and coefficients
  if (inherits(equation, "symbolic_equation")) {
    expr_str <- equation$expression
    coefs <- stats::coef(equation)
  } else if (inherits(equation, "nls")) {
    expr_str <- deparse(stats::formula(equation)[[3]])
    coefs <- stats::coef(equation)
  } else if (inherits(equation, "lm")) {
    # For lm, try to construct formula representation
    # Handle intercept-only model specially
    if (length(stats::coef(equation)) == 1 && names(stats::coef(equation))[1] == "(Intercept)") {
      expr_str <- as.character(stats::coef(equation)[1])
      coefs <- numeric(0) # No replacements needed beyond the initial value
    } else {
      expr_str <- deparse(stats::formula(equation)[[3]])
      coefs <- stats::coef(equation)
    }
  } else if (inherits(equation, "sde_model")) {
    # Format SDE with drift and diffusion
    drift_latex <- to_latex(equation$drift, variable, precision, 
                            scientific_notation, include_uncertainty, se_values)
    diff_latex <- to_latex(equation$diffusion, variable, precision,
                           scientific_notation, FALSE, NULL)
    
    latex <- sprintf("d%s = \\left(%s\\right) dt + \\left(%s\\right) dW_t",
                     variable, drift_latex, diff_latex)
    return(latex)
    
  } else if (inherits(equation, "variance_model")) {
    # Handle variance models (used in diffusion)
    if (equation$method == "constant") {
      # For constant method, the fit is an intercept-only LM
      # We want the standard deviation (sigma), not variance
      val <- stats::coef(equation$fit)[1]
      
      # Adjust scale based on transform used during fitting
      if (equation$transform == "squared") val <- sqrt(max(0, val))
      else if (equation$transform == "log_squared") val <- sqrt(exp(val))
      else val <- max(0, val) # absolute
      
      return(sprintf("%.*f", precision, val))
    } else {
      # Delegate to internal fit object for linear/quadratic/symbolic models
      return(to_latex(equation$fit, variable, precision, scientific_notation))
    }
    
  } else {
    stop("Unknown equation type: ", class(equation)[1])
  }
  
  # Replace coefficients with formatted values
  latex_expr <- expr_str
  
  for (nm in names(coefs)) {
    val <- coefs[nm]
    
    # Format the coefficient
    if (scientific_notation && (abs(val) > 1000 || (abs(val) < 0.001 && val != 0))) {
      # Scientific notation
      exponent <- floor(log10(abs(val)))
      mantissa <- val / 10^exponent
      formatted <- sprintf("%.*f \\times 10^{%d}", precision, mantissa, exponent)
    } else {
      formatted <- sprintf("%.*f", precision, val)
    }
    
    # Add uncertainty if requested
    if (include_uncertainty && !is.null(se_values) && nm %in% names(se_values)) {
      formatted <- sprintf("%s\\,(%.*f)", formatted, precision, se_values[nm])
    }
    
    # Replace parameter name with value
    # Escape parentheses for regex if present in names (like "(Intercept)")
    safe_nm <- gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", nm))
    
    # Use word boundaries to avoid partial replacements
    latex_expr <- gsub(paste0("\\b", safe_nm, "\\b"), formatted, latex_expr)
    
    # Handle Intercept special case in lm formulas explicitly if regex missed it
    if (nm == "(Intercept)") {
      latex_expr <- sub("1", formatted, latex_expr)
    }
  }
  
  # Convert R syntax to LaTeX
  latex_expr <- r_to_latex_expr(latex_expr, variable)
  
  # Wrap in equation notation ONLY if this isn't a sub-call (heuristic)
  # Ideally, we'd pass a flag, but for now, we assume simple calls want the full equation
  # If called recursively (from SDE block), we just return the expression part.
  # Since sde_model block calls recursively and constructs its own equation, 
  # we only wrap here if it's NOT an SDE model (which we already handled above)
  # or a Variance model (which returns raw expression).
  
  # However, consistent behavior: to_latex usually returns right-hand side 
  # unless it's the top-level SDE which returns full dZ = ...
  
  # If the input was just an equation (lm, nls, symbolic), usually we want "dZ/dt = ..."
  # But if this function is called recursively for the diffusion term, we don't want "dZ/dt ="
  
  # Current simple logic: Always return just the LaTeX expression string.
  # EXCEPT for SDE which formats the full differential.
  # If user wants "dZ/dt =", they can add it or use format_equation().
  
  # Re-reading original code: it wrapped in dZ/dt. Let's keep that behavior 
  # but ONLY if not being called recursively? 
  # Actually, the SDE block calls to_latex() and expects the expression part.
  # So we should return just latex_expr here?
  # Original code: wrapped it.
  # Let's verify: SDE block does: drift_latex <- to_latex(...)
  # If to_latex returns "dZ/dt = ...", the SDE string becomes "dZ = (dZ/dt = ...) dt + ...", which is wrong.
  
  # CORRECTION: We need to know if we are in recursive mode.
  # Since we can't easily change the signature without breaking docs/tests,
  # let's assume if it's an SDE model, we return full equation.
  # If it's a component model, we verify.
  
  # NOTE: To fix the SDE formatting bug, we will just return the expression (RHS) here.
  # The user-facing function might need to add LHS if desired, or we rely on the user understanding.
  # HOWEVER, the original function wrapped it. Let's support both via a heuristic:
  # The SDE block constructs the full equation itself.
  
  return(latex_expr)
}

#' Convert R Expression to LaTeX
#' @keywords internal
r_to_latex_expr <- function(expr, variable = "Z") {
  latex <- expr
  
  # Powers: Z^2 -> Z^{2}
  latex <- gsub("([A-Za-z_][A-Za-z0-9_]*)\\^([0-9]+)", "\\1^{\\2}", latex)
  latex <- gsub("([A-Za-z_][A-Za-z0-9_]*)\\^\\(([^)]+)\\)", "\\1^{\\2}", latex)
  
  # Multiplication: * -> \cdot or space
  latex <- gsub("\\*", " \\\\cdot ", latex)
  
  # Division: / -> \frac{}{}
  # Simple cases: a/b -> \frac{a}{b}
  latex <- gsub("([A-Za-z0-9_.]+)/([A-Za-z0-9_.]+)", "\\\\frac{\\1}{\\2}", latex)
  
  # Functions
  latex <- gsub("\\bsqrt\\(", "\\\\sqrt{", latex)
  latex <- gsub("\\bexp\\(", "\\\\exp\\\\left(", latex)
  latex <- gsub("\\blog\\(", "\\\\ln\\\\left(", latex)
  latex <- gsub("\\bsin\\(", "\\\\sin\\\\left(", latex)
  latex <- gsub("\\bcos\\(", "\\\\cos\\\\left(", latex)
  latex <- gsub("\\btan\\(", "\\\\tan\\\\left(", latex)
  
  # Close function parentheses with \right)
  # This is simplified - complex expressions may need manual adjustment
  latex <- gsub("\\\\left\\(([^)]+)\\)", "\\\\left(\\1\\\\right)", latex)
  
  # Greek letters (common parameter names)
  latex <- gsub("\\balpha\\b", "\\\\alpha", latex)
  latex <- gsub("\\bbeta\\b", "\\\\beta", latex)
  latex <- gsub("\\bgamma\\b", "\\\\gamma", latex)
  latex <- gsub("\\bdelta\\b", "\\\\delta", latex)
  latex <- gsub("\\bsigma\\b", "\\\\sigma", latex)
  latex <- gsub("\\bmu\\b", "\\\\mu", latex)
  latex <- gsub("\\blambda\\b", "\\\\lambda", latex)
  latex <- gsub("\\btheta\\b", "\\\\theta", latex)
  latex <- gsub("\\brho\\b", "\\\\rho", latex)
  latex <- gsub("\\bphi\\b", "\\\\phi", latex)
  
  return(latex)
}

#' Format Equation for Display
#'
#' Creates a nicely formatted string representation of the equation.
#'
#' @param equation Fitted equation object.
#' @param format Output format: "text", "latex", "markdown".
#' @param precision Number of decimal places.
#'
#' @return Formatted string.
#' @export
format_equation <- function(equation, format = c("text", "latex", "markdown"),
                            precision = 4) {
  format <- match.arg(format)
  
  # Handle SDE special case for text/markdown
  if (inherits(equation, "sde_model")) {
    if (format == "latex") return(to_latex(equation, precision = precision))
    
    # Text representation
    drift_str <- format_equation(equation$drift, "text", precision)
    diff_str <- if(!is.null(equation$diffusion)) {
      if(inherits(equation$diffusion, "variance_model")) {
        # Handle simple constant case cleanly
        if(equation$diffusion$method == "constant") {
          val <- stats::coef(equation$diffusion$fit)[1]
          if(equation$diffusion$transform == "squared") val <- sqrt(max(0,val))
          sprintf("%.*f", precision, val)
        } else {
          format_equation(equation$diffusion$fit, "text", precision)
        }
      } else "g(.)"
    } else "sigma"
    
    var_name <- equation$variable %||% "Z"
    return(sprintf("d%s = [%s] dt + [%s] dW", var_name, drift_str, diff_str))
  }
  
  coefs <- stats::coef(equation)
  
  if (inherits(equation, "symbolic_equation")) {
    expr <- equation$expression
  } else if (inherits(equation, "nls") || inherits(equation, "lm")) {
    expr <- deparse(stats::formula(equation)[[3]])
  } else {
    expr <- "Unknown"
  }
  
  # Substitute coefficient values
  for (nm in names(coefs)) {
    val <- sprintf("%.*g", precision, coefs[nm])
    expr <- gsub(paste0("\\b", nm, "\\b"), val, expr)
  }
  
  if (format == "latex") {
    # Note: to_latex returns the RHS expression, add LHS manually here if needed
    return(to_latex(equation, precision = precision))
  } else if (format == "markdown") {
    return(paste0("$", to_latex(equation, precision = precision), "$"))
  } else {
    return(expr)
  }
}


# =============================================================================
# TABLE GENERATION
# =============================================================================

#' Generate Coefficient Table
#'
#' Creates a publication-ready table of estimated coefficients.
#'
#' @param equation Fitted equation object.
#' @param bootstrap_results Results from \code{bootstrap_parameters} (optional).
#' @param format Output format: "data.frame", "latex", "markdown", "html".
#' @param caption Table caption.
#' @param label LaTeX label for referencing.
#'
#' @return Formatted table.
#' @export
coefficient_table <- function(equation, bootstrap_results = NULL,
                              format = c("data.frame", "latex", "markdown", "html"),
                              caption = "Estimated Coefficients",
                              label = "tab:coefficients") {
  format <- match.arg(format)
  
  coefs <- stats::coef(equation)
  
  # Try to get standard errors
  se <- tryCatch({
    if (inherits(equation, "nls")) {
      summary(equation)$parameters[, "Std. Error"]
    } else if (inherits(equation, "lm")) {
      summary(equation)$coefficients[, "Std. Error"]
    } else {
      rep(NA, length(coefs))
    }
  }, error = function(e) rep(NA, length(coefs)))
  
  # Use bootstrap results if provided
  if (!is.null(bootstrap_results)) {
    se <- bootstrap_results$se
    ci_lower <- bootstrap_results$ci_lower
    ci_upper <- bootstrap_results$ci_upper
  } else {
    ci_lower <- coefs - 1.96 * se
    ci_upper <- coefs + 1.96 * se
  }
  
  # Create data frame
  df <- data.frame(
    Parameter = names(coefs),
    Estimate = coefs,
    `Std. Error` = se,
    `CI Lower` = ci_lower,
    `CI Upper` = ci_upper,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (format == "data.frame") {
    return(df)
  }
  
  # Format numbers
  df$Estimate <- sprintf("%.4f", df$Estimate)
  df$`Std. Error` <- ifelse(is.na(df$`Std. Error`), "--", sprintf("%.4f", df$`Std. Error`))
  df$`CI Lower` <- ifelse(is.na(df$`CI Lower`), "--", sprintf("%.4f", df$`CI Lower`))
  df$`CI Upper` <- ifelse(is.na(df$`CI Upper`), "--", sprintf("%.4f", df$`CI Upper`))
  
  if (format == "latex") {
    return(df_to_latex(df, caption, label))
  } else if (format == "markdown") {
    return(df_to_markdown(df))
  } else if (format == "html") {
    return(df_to_html(df, caption))
  }
}

#' Generate Model Comparison Table
#'
#' Creates a table comparing multiple candidate equations.
#'
#' @param equations Named list of fitted equation objects.
#' @param data Data for computing fit statistics.
#' @param derivative_col Response variable column.
#' @param format Output format.
#' @param caption Table caption.
#'
#' @return Comparison table.
#' @export
model_comparison_table <- function(equations, data, derivative_col,
                                   format = c("data.frame", "latex", "markdown"),
                                   caption = "Model Comparison") {
  format <- match.arg(format)
  
  results <- lapply(names(equations), function(nm) {
    eq <- equations[[nm]]
    
    # Compute metrics
    pred <- tryCatch(stats::predict(eq, newdata = data), error = function(e) rep(NA, nrow(data)))
    actual <- data[[derivative_col]]
    resid <- actual - pred
    
    n <- sum(!is.na(resid))
    k <- length(stats::coef(eq))
    
    ss_res <- sum(resid^2, na.rm = TRUE)
    ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - ss_res / ss_tot
    adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - k - 1)
    rmse <- sqrt(mean(resid^2, na.rm = TRUE))
    
    aic <- tryCatch(stats::AIC(eq), error = function(e) NA)
    bic <- tryCatch(stats::BIC(eq), error = function(e) NA)
    
    data.frame(
      Model = nm,
      Parameters = k,
      `R-squared` = r2,
      `Adj. R-squared` = adj_r2,
      RMSE = rmse,
      AIC = aic,
      BIC = bic,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  
  # Sort by BIC
  df <- df[order(df$BIC), ]
  
  if (format == "data.frame") {
    return(df)
  } else if (format == "latex") {
    return(df_to_latex(df, caption, "tab:comparison"))
  } else {
    return(df_to_markdown(df))
  }
}

#' Data Frame to LaTeX Table
#' @keywords internal
df_to_latex <- function(df, caption = NULL, label = NULL) {
  n_col <- ncol(df)
  
  # Header
  latex <- "\\begin{table}[htbp]\n\\centering\n"
  if (!is.null(caption)) {
    latex <- paste0(latex, "\\caption{", caption, "}\n")
  }
  if (!is.null(label)) {
    latex <- paste0(latex, "\\label{", label, "}\n")
  }
  
  # Column specification
  col_spec <- paste(rep("c", n_col), collapse = "")
  latex <- paste0(latex, "\\begin{tabular}{", col_spec, "}\n\\hline\n")
  
  # Header row
  latex <- paste0(latex, paste(names(df), collapse = " & "), " \\\\\n\\hline\n")
  
  # Data rows
  for (i in 1:nrow(df)) {
    row_str <- paste(df[i, ], collapse = " & ")
    latex <- paste0(latex, row_str, " \\\\\n")
  }
  
  latex <- paste0(latex, "\\hline\n\\end{tabular}\n\\end{table}")
  
  return(latex)
}

#' Data Frame to Markdown Table
#' @keywords internal
df_to_markdown <- function(df) {
  # Header
  md <- paste0("| ", paste(names(df), collapse = " | "), " |\n")
  md <- paste0(md, "| ", paste(rep("---", ncol(df)), collapse = " | "), " |\n")
  
  # Data rows
  for (i in 1:nrow(df)) {
    md <- paste0(md, "| ", paste(df[i, ], collapse = " | "), " |\n")
  }
  
  return(md)
}

#' Data Frame to HTML Table
#' @keywords internal
df_to_html <- function(df, caption = NULL) {
  html <- "<table class='table'>\n"
  
  if (!is.null(caption)) {
    html <- paste0(html, "<caption>", caption, "</caption>\n")
  }
  
  # Header
  html <- paste0(html, "<thead><tr>\n")
  for (nm in names(df)) {
    html <- paste0(html, "<th>", nm, "</th>\n")
  }
  html <- paste0(html, "</tr></thead>\n<tbody>\n")
  
  # Data rows
  for (i in 1:nrow(df)) {
    html <- paste0(html, "<tr>\n")
    for (j in 1:ncol(df)) {
      html <- paste0(html, "<td>", df[i, j], "</td>\n")
    }
    html <- paste0(html, "</tr>\n")
  }
  
  html <- paste0(html, "</tbody>\n</table>")
  
  return(html)
}


# =============================================================================
# COMPREHENSIVE REPORTS
# =============================================================================

#' Generate Analysis Report
#'
#' Creates a comprehensive report of the entire analysis workflow.
#'
#' @param results List containing analysis results with elements:
#'    \itemize{
#'      \item data: Original data frame
#'      \item derivatives: Computed derivatives
#'      \item exploration: Results from \code{explore_dynamics}
#'      \item equation: Best fitted equation
#'      \item sde: SDE model (optional)
#'      \item validation: Results from \code{validate_model}
#'    }
#' @param output_file Path for output file (required, no default to comply with CRAN policy).
#' @param format Report format: "markdown", "html", "latex".
#' @param title Report title.
#' @param author Author name.
#' @param include_plots Include diagnostic plots.
#'
#' @return Path to generated report.
#' @export
#' @examples
#' \donttest{
#' # Toy example to demonstrate report generation
#' # Using a temporary file to avoid writing to user's working directory
#' tmp_file <- tempfile("report_example")
#' 
#' # Mock results object
#' mock_results <- list(
#'   data = data.frame(time = 1:10, Z = runif(10)),
#'   equation = stats::lm(Z ~ time, data = data.frame(time = 1:10, Z = runif(10)))
#' )
#' 
#' # Generate report
#' report_path <- generate_report(mock_results, output_file = tmp_file, format = "markdown")
#' if(file.exists(report_path)) unlink(report_path)
#' }
generate_report <- function(results, output_file,
                            format = c("markdown", "html", "latex"),
                            title = "Empirical Dynamics Analysis Report",
                            author = "EmpiricalDynamics",
                            include_plots = TRUE) {
  format <- match.arg(format)
  
  # Build report content
  report <- character(0)
  
  # Header
  if (format == "latex") {
    report <- c(report, 
                "\\documentclass{article}",
                "\\usepackage{amsmath,graphicx,booktabs}",
                sprintf("\\title{%s}", title),
                sprintf("\\author{%s}", author),
                "\\date{\\today}",
                "\\begin{document}",
                "\\maketitle"
    )
  } else if (format == "markdown") {
    report <- c(report,
                sprintf("# %s", title),
                sprintf("**Author:** %s", author),
                sprintf("**Date:** %s", Sys.Date()),
                ""
    )
  } else if (format == "html") {
    report <- c(report,
                "<!DOCTYPE html>",
                "<html><head>",
                sprintf("<title>%s</title>", title),
                "<style>",
                "body { font-family: Arial, sans-serif; max-width: 800px; margin: auto; padding: 20px; }",
                "table { border-collapse: collapse; width: 100%; }",
                "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
                "th { background-color: #f2f2f2; }",
                "</style>",
                "</head><body>",
                sprintf("<h1>%s</h1>", title),
                sprintf("<p><strong>Author:</strong> %s<br>", author),
                sprintf("<strong>Date:</strong> %s</p>", Sys.Date())
    )
  }
  
  # Section 1: Data Summary
  report <- c(report, section_header("Data Summary", format))
  
  if (!is.null(results$data)) {
    n_obs <- nrow(results$data)
    n_vars <- ncol(results$data)
    
    if (format == "latex") {
      report <- c(report,
                  sprintf("Number of observations: %d\\\\", n_obs),
                  sprintf("Number of variables: %d\\\\", n_vars)
      )
    } else if (format == "markdown") {
      report <- c(report,
                  sprintf("- **Observations:** %d", n_obs),
                  sprintf("- **Variables:** %d", n_vars),
                  ""
      )
    } else {
      report <- c(report,
                  sprintf("<ul><li><strong>Observations:</strong> %d</li>", n_obs),
                  sprintf("<li><strong>Variables:</strong> %d</li></ul>", n_vars)
      )
    }
  }
  
  # Section 2: Discovered Equation
  report <- c(report, section_header("Discovered Equation", format))
  
  if (!is.null(results$equation)) {
    # Use to_latex here to get just the RHS expression
    eq_latex <- to_latex(results$equation)
    # Add LHS manually for clarity in report
    target_var <- results$equation$target %||% "dZ"
    full_eq_latex <- sprintf("\\frac{d%s}{dt} = %s", "Z", eq_latex) # Assuming Z, user can refine
    
    if (format == "latex") {
      report <- c(report,
                  "\\begin{equation}",
                  full_eq_latex,
                  "\\end{equation}"
      )
    } else if (format == "markdown") {
      report <- c(report,
                  "$$",
                  full_eq_latex,
                  "$$",
                  ""
      )
    } else {
      report <- c(report,
                  sprintf("<p>$$%s$$</p>", full_eq_latex)
      )
    }
    
    # Coefficient table
    coef_table <- coefficient_table(results$equation, format = format)
    report <- c(report, coef_table)
  }
  
  # Section 3: SDE Model
  if (!is.null(results$sde)) {
    report <- c(report, section_header("Stochastic Differential Equation", format))
    
    sde_latex <- to_latex(results$sde)
    
    if (format == "latex") {
      report <- c(report,
                  "\\begin{equation}",
                  sde_latex,
                  "\\end{equation}"
      )
    } else if (format == "markdown") {
      report <- c(report,
                  "$$",
                  sde_latex,
                  "$$",
                  ""
      )
    } else {
      report <- c(report,
                  sprintf("<p>$$%s$$</p>", sde_latex)
      )
    }
  }
  
  # Section 4: Validation Results
  if (!is.null(results$validation)) {
    report <- c(report, section_header("Model Validation", format))
    
    val <- results$validation
    
    metrics <- data.frame(
      Metric = c("Cross-validation RMSE", "In-sample R-squared", 
                 "AIC", "BIC"),
      Value = c(
        sprintf("%.4f", val$summary$cv_rmse),
        sprintf("%.4f", val$summary$in_sample_r2),
        sprintf("%.2f", val$fit_stats$aic),
        sprintf("%.2f", val$fit_stats$bic)
      ),
      stringsAsFactors = FALSE
    )
    
    if (format == "latex") {
      report <- c(report, df_to_latex(metrics, "Validation Metrics", "tab:validation"))
    } else if (format == "markdown") {
      report <- c(report, df_to_markdown(metrics))
    } else {
      report <- c(report, df_to_html(metrics, "Validation Metrics"))
    }
  }
  
  # Section 5: Qualitative Analysis
  if (!is.null(results$validation$qualitative)) {
    report <- c(report, section_header("Qualitative Analysis", format))
    
    qual <- results$validation$qualitative
    
    if (!is.null(qual$checks$fixed_points) && nrow(qual$checks$fixed_points) > 0) {
      fps <- qual$checks$fixed_points
      fps$fixed_point <- sprintf("%.4f", fps$fixed_point)
      fps$eigenvalue <- sprintf("%.4f", fps$eigenvalue)
      
      if (format == "latex") {
        report <- c(report, df_to_latex(fps, "Fixed Points", "tab:fixedpoints"))
      } else if (format == "markdown") {
        report <- c(report, "### Fixed Points", "", df_to_markdown(fps))
      } else {
        report <- c(report, "<h3>Fixed Points</h3>", df_to_html(fps))
      }
    }
  }
  
  # Close document
  if (format == "latex") {
    report <- c(report, "\\end{document}")
    ext <- ".tex"
  } else if (format == "html") {
    report <- c(report, "</body></html>")
    ext <- ".html"
  } else {
    ext <- ".md"
  }
  
  # Write to file
  output_path <- paste0(output_file, ext)
  writeLines(report, output_path)
  
  message("Report written to: ", output_path)
  
  return(output_path)
}

#' Create Section Header
#' @keywords internal
section_header <- function(title, format) {
  if (format == "latex") {
    return(sprintf("\\section{%s}", title))
  } else if (format == "markdown") {
    return(c("", sprintf("## %s", title), ""))
  } else {
    return(sprintf("<h2>%s</h2>", title))
  }
}


# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

#' Export Results to Multiple Formats
#'
#' Exports analysis results to various file formats.
#'
#' @param results Analysis results list.
#' @param output_dir Output directory (required, no default to comply with CRAN policy).
#' @param prefix File name prefix.
#' @param formats Vector of formats: "rds", "csv", "json", "latex".
#'
#' @return List of paths to created files.
#' @export
#' @examples
#' \donttest{
#' # Toy example
#' tmp_dir <- tempdir()
#' mock_results <- list(
#'   equation = stats::lm(mpg ~ wt, data = mtcars)
#' )
#' 
#' # Export
#' paths <- export_results(mock_results, output_dir = tmp_dir, formats = c("csv", "rds"))
#' }
export_results <- function(results, output_dir,
                           prefix = "empirical_dynamics",
                           formats = c("rds", "csv")) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  paths <- list()
  
  # RDS (full R object)
  if ("rds" %in% formats) {
    path <- file.path(output_dir, paste0(prefix, "_results.rds"))
    saveRDS(results, path)
    paths$rds <- path
  }
  
  # CSV (coefficient table)
  if ("csv" %in% formats && !is.null(results$equation)) {
    path <- file.path(output_dir, paste0(prefix, "_coefficients.csv"))
    coef_df <- coefficient_table(results$equation, format = "data.frame")
    utils::write.csv(coef_df, path, row.names = FALSE)
    paths$csv_coef <- path
    
    # Also export validation metrics if available
    if (!is.null(results$validation)) {
      path <- file.path(output_dir, paste0(prefix, "_validation.csv"))
      val_df <- data.frame(
        metric = names(results$validation$summary),
        value = unlist(results$validation$summary)
      )
      utils::write.csv(val_df, path, row.names = FALSE)
      paths$csv_val <- path
    }
  }
  
  # JSON
  if ("json" %in% formats) {
    path <- file.path(output_dir, paste0(prefix, "_results.json"))
    
    # Create JSON-friendly structure
    json_data <- list(
      equation = if (!is.null(results$equation)) {
        list(
          expression = if (inherits(results$equation, "symbolic_equation")) 
            results$equation$expression else deparse(stats::formula(results$equation)[[3]]),
          coefficients = as.list(stats::coef(results$equation))
        )
      } else NULL,
      validation = if (!is.null(results$validation)) {
        results$validation$summary
      } else NULL
    )
    
    json_str <- to_json(json_data)
    writeLines(json_str, path)
    paths$json <- path
  }
  
  # LaTeX
  if ("latex" %in% formats && !is.null(results$equation)) {
    path <- file.path(output_dir, paste0(prefix, "_equation.tex"))
    latex <- to_latex(results$equation)
    writeLines(latex, path)
    paths$latex <- path
  }
  
  message("Exported results to: ", output_dir)
  
  return(paths)
}

#' Simple JSON Conversion
#' @keywords internal
to_json <- function(x, indent = 0) {
  spaces <- paste(rep("  ", indent), collapse = "")
  
  if (is.null(x)) {
    return("null")
  } else if (is.list(x) && !is.null(names(x))) {
    # Named list -> object
    items <- sapply(names(x), function(nm) {
      sprintf('"%s": %s', nm, to_json(x[[nm]], indent + 1))
    })
    return(paste0("{\n", spaces, "  ", 
                  paste(items, collapse = paste0(",\n", spaces, "  ")),
                  "\n", spaces, "}"))
  } else if (is.list(x)) {
    # Unnamed list -> array
    items <- sapply(x, to_json, indent + 1)
    return(paste0("[", paste(items, collapse = ", "), "]"))
  } else if (is.character(x)) {
    return(sprintf('"%s"', gsub('"', '\\"', x)))
  } else if (is.numeric(x)) {
    if (length(x) == 1) {
      return(sprintf("%.6g", x))
    } else {
      return(paste0("[", paste(sprintf("%.6g", x), collapse = ", "), "]"))
    }
  } else if (is.logical(x)) {
    return(tolower(as.character(x)))
  } else {
    return(sprintf('"%s"', as.character(x)))
  }
}


# =============================================================================
# PLOT EXPORT
# =============================================================================

#' Save Diagnostic Plots
#'
#' Saves all diagnostic plots to files.
#'
#' @param results Analysis results list.
#' @param output_dir Output directory (required, no default to comply with CRAN policy).
#' @param prefix File name prefix.
#' @param format Image format: "png", "pdf", "svg".
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Resolution for raster formats.
#'
#' @return List of paths to created files.
#' @export
save_plots <- function(results, output_dir,
                       prefix = "empirical_dynamics",
                       format = c("png", "pdf"),
                       width = 8, height = 6, dpi = 300) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot export.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  paths <- list()
  
  # Helper function to save a plot
  save_one <- function(plot_expr, name) {
    for (fmt in format) {
      path <- file.path(output_dir, paste0(prefix, "_", name, ".", fmt))
      
      if (fmt == "pdf") {
        grDevices::pdf(path, width = width, height = height)
      } else if (fmt == "png") {
        grDevices::png(path, width = width * dpi, height = height * dpi, res = dpi)
      } else if (fmt == "svg") {
        grDevices::svg(path, width = width, height = height)
      }
      
      tryCatch({
        print(eval(plot_expr))
      }, error = function(e) {
        message("Failed to create plot '", name, "': ", e$message)
      })
      
      grDevices::dev.off()
      paths[[paste0(name, "_", fmt)]] <- path
    }
  }
  
  # Exploration plots
  if (!is.null(results$exploration)) {
    # Phase diagrams
    if (!is.null(results$exploration$phase_plots)) {
      for (i in seq_along(results$exploration$phase_plots)) {
        p <- results$exploration$phase_plots[[i]]
        save_one(quote(p), paste0("phase_", i))
      }
    }
  }
  
  # Validation plots
  if (!is.null(results$validation)) {
    # CV results
    if (!is.null(results$validation$cv)) {
      save_one(quote(plot(results$validation$cv)), "cv_results")
    }
    
    # Residual diagnostics
    if (!is.null(results$validation$residuals)) {
      save_one(quote(plot(results$validation$residuals)), "residual_diagnostics")
    }
    
    # Trajectory comparison
    if (!is.null(results$validation$trajectory)) {
      save_one(quote(plot(results$validation$trajectory$simulation)), "trajectory")
    }
  }
  
  message("Saved plots to: ", output_dir)
  
  return(paths)
}


# =============================================================================
# SUMMARY AND PRINTING
# =============================================================================

#' Create Analysis Summary
#'
#' Creates a concise text summary of the analysis.
#'
#' @param results Analysis results list.
#' @param verbose Include additional details.
#'
#' @return Character string with summary.
#' @export
analysis_summary <- function(results, verbose = TRUE) {
  
  lines <- character(0)
  
  lines <- c(lines,
             "=================================================================",
             "          EMPIRICAL DYNAMICS ANALYSIS SUMMARY",
             "=================================================================",
             ""
  )
  
  # Data info
  if (!is.null(results$data)) {
    lines <- c(lines,
               sprintf("Data: %d observations, %d variables", 
                       nrow(results$data), ncol(results$data)),
               ""
    )
  }
  
  # Discovered equation
  if (!is.null(results$equation)) {
    lines <- c(lines,
               "DISCOVERED EQUATION:",
               "--------------------"
    )
    
    eq_str <- format_equation(results$equation, format = "text")
    lines <- c(lines, sprintf("  dZ/dt = %s", eq_str), "")
    
    if (verbose) {
      coefs <- stats::coef(results$equation)
      lines <- c(lines, "  Coefficients:")
      for (nm in names(coefs)) {
        lines <- c(lines, sprintf("    %s = %.6g", nm, coefs[nm]))
      }
      lines <- c(lines, "")
    }
  }
  
  # SDE
  if (!is.null(results$sde)) {
    lines <- c(lines,
               "STOCHASTIC MODEL:",
               "-----------------"
    )
    
    drift_str <- format_equation(results$sde$drift, format = "text")
    diff_str <- if (!is.null(results$sde$diffusion)) {
      if(inherits(results$sde$diffusion, "variance_model") && results$sde$diffusion$method == "constant") {
        # Simple formatting for constant
        val <- stats::coef(results$sde$diffusion$fit)[1]
        if(results$sde$diffusion$transform == "squared") val <- sqrt(max(0, val))
        sprintf("%.4f", val)
      } else {
        # Generic case
        "g(.)"
      }
    } else {
      "sigma"
    }
    
    lines <- c(lines,
               sprintf("  Drift:     %s", drift_str),
               sprintf("  Diffusion: %s", diff_str),
               ""
    )
  }
  
  # Validation
  if (!is.null(results$validation)) {
    lines <- c(lines,
               "VALIDATION:",
               "-----------"
    )
    
    val <- results$validation
    
    lines <- c(lines,
               sprintf("  CV RMSE:        %.4f", val$summary$cv_rmse),
               sprintf("  In-sample R^2:   %.4f", val$summary$in_sample_r2)
    )
    
    if (!is.null(val$fit_stats)) {
      lines <- c(lines,
                 sprintf("  AIC:            %.2f", val$fit_stats$aic),
                 sprintf("  BIC:            %.2f", val$fit_stats$bic)
      )
    }
    
    # Residual tests
    if (!is.null(val$residuals) && verbose) {
      n_pass <- sum(val$residuals$tests$p_value > 0.05)
      n_total <- nrow(val$residuals$tests)
      lines <- c(lines,
                 sprintf("  Residual tests: %d/%d passed", n_pass, n_total)
      )
    }
    
    # Qualitative
    if (!is.null(val$qualitative)) {
      lines <- c(lines,
                 sprintf("  Qualitative:    %d/%d checks passed",
                         val$qualitative$n_passed, val$qualitative$n_total)
      )
    }
    
    lines <- c(lines, "")
  }
  
  # Fixed points
  if (!is.null(results$validation$qualitative$checks$fixed_points)) {
    fps <- results$validation$qualitative$checks$fixed_points
    if (nrow(fps) > 0) {
      lines <- c(lines,
                 "FIXED POINTS:",
                 "-------------"
      )
      for (i in 1:nrow(fps)) {
        lines <- c(lines,
                   sprintf("  Z* = %.4f (%s, lambda = %.4f)",
                           fps$fixed_point[i], fps$stability[i], fps$eigenvalue[i])
        )
      }
      lines <- c(lines, "")
    }
  }
  
  lines <- c(lines,
             "=================================================================",
             sprintf("Generated: %s", Sys.time()),
             "================================================================="
  )
  
  summary_text <- paste(lines, collapse = "\n")
  
  return(summary_text)
}

#' Print Analysis Summary
#' @param results Analysis results list.
#' @return Invisibly returns the input results object (called for side effects).
#' @export
print_summary <- function(results) {
  cat(analysis_summary(results), "\n")
  invisible(results)
}


# =============================================================================
# WORKFLOW TEMPLATE
# =============================================================================

#' Get Analysis Template
#'
#' Returns a template script for running a complete analysis.
#'
#' @param output_file Path to save template.
#'
#' @return Template code as character string (invisibly).
#' @export
#' @examples
#' \donttest{
#' # Save template to a temporary file
#' tmp_file <- tempfile("analysis_template", fileext = ".R")
#' get_analysis_template(tmp_file)
#' 
#' # Clean up
#' if (file.exists(tmp_file)) unlink(tmp_file)
#' }
get_analysis_template <- function(output_file = NULL) {
  
  template <- '
# =============================================================================
# EmpiricalDynamics Analysis Template
# =============================================================================

library(EmpiricalDynamics)

# -----------------------------------------------------------------------------
# Step A: Load and Preprocess Data
# -----------------------------------------------------------------------------

# Load your data
data <- read_empirical_data("your_data.csv")

# Specify variables
vars <- specify_variables(
  data = data,
  endogenous = "Z",            # Main variable of interest
  exogenous = c("X", "Y"),     # External drivers
  time_col = "time"
)

# Compute derivatives (TVR recommended for economic data)
deriv <- compute_derivative(
  data[[vars$endogenous]],
  time = data[[vars$time_col]],
  method = "tvr",
  lambda = "auto"
)

# Add derivative to data
data$dZ <- deriv$derivative

# -----------------------------------------------------------------------------
# Step B: Visual Exploration
# -----------------------------------------------------------------------------

# Explore relationships
exploration <- explore_dynamics(
  data = data,
  response = "dZ",
  predictors = c("Z", "X", "Y")
)

# Review plots and suggestions
print(exploration$suggestions)

# -----------------------------------------------------------------------------
# Step C: Equation Discovery
# -----------------------------------------------------------------------------

# Option 1: Symbolic search (automated)
search_result <- symbolic_search(
  data = data,
  response = "dZ",
  predictors = c("Z", "X", "Y"),
  max_complexity = 20,
  n_runs = 3
)

# Select best equation
best_eq <- select_equation(search_result, criterion = "bic")

# Option 2: Fit specified equation (theory-guided)
# eq <- fit_specified_equation(
#   "a + b * Z + c * Z^2 + d * X",
#   data = data,
#   derivative_col = "dZ"
# )

# -----------------------------------------------------------------------------
# Step D: Residual Analysis and SDE Construction
# -----------------------------------------------------------------------------

# Compute residuals
resid <- compute_residuals(best_eq, data, "dZ")

# Run diagnostics
diag <- residual_diagnostics(resid, data)
print(diag)

# Model conditional variance (diffusion)
var_model <- model_conditional_variance(
  resid,
  data = data,
  predictors = c("Z"),
  method = "linear"
)

# Construct SDE
sde <- construct_sde(
  drift = best_eq,
  diffusion = var_model,
  variable = "Z"
)

# -----------------------------------------------------------------------------
# Step E: Validation
# -----------------------------------------------------------------------------

# Comprehensive validation
validation <- validate_model(
  equation = best_eq,
  sde = sde,
  data = data,
  derivative_col = "dZ",
  variable = "Z",
  cv_folds = 5
)

print(validation)

# -----------------------------------------------------------------------------
# Step F: Output and Reporting
# -----------------------------------------------------------------------------
  
# Collect all results
results <- list(
  data = data,
  exploration = exploration,
  equation = best_eq,
  sde = sde,
  validation = validation
)

# Generate report
# NOTE: Replace with your desired output path
generate_report(
  results,
  output_file = file.path(tempdir(), "analysis_report"),  # Change path as needed
  format = "markdown",
  title = "My Empirical Dynamics Analysis"
)

# Export results
# NOTE: Replace with your desired output directory
export_results(results, output_dir = tempdir(), formats = c("rds", "csv"))

# Print summary
print_summary(results)

# Get LaTeX equation for publication
cat("\\n\\nLaTeX equation:\\n")
cat(to_latex(best_eq), "\\n")
'

if (!is.null(output_file)) {
  writeLines(template, output_file)
  message("Template saved to: ", output_file)
}

cat(template)
invisible(template)
}