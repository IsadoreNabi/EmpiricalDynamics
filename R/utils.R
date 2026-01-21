#' @title Utility Functions for EmpiricalDynamics
#' @description Internal utility functions and helpers used across the package.
#' @name utils
#' @keywords internal
#' 
#' @importFrom methods as
#' @importFrom stats acf coef confint cor fft fitted lm mad median mvfft na.omit nls nls.control optim predict qnorm quantile resid rnorm runif sd smooth.spline spline var AIC BIC logLik residuals simulate ts Box.test anova approx as.formula dnorm dt formula loess pchisq pnorm setNames spectrum time uniroot qqline qqnorm
#' @importFrom graphics abline grid legend lines mtext par plot points polygon rect text title contour persp image filled.contour curve hist plot.new
#' @importFrom grDevices colorRampPalette dev.off pdf png rainbow adjustcolor svg
#' @importFrom utils head tail write.csv modifyList
NULL

utils::globalVariables(c(
  ".data", "actual", "predicted", "fold", "rmse", "time", "value",
  "q05", "q25", "q50", "q75", "q95", "parameter_value", "fixed_point",
  "stability", "complexity", "equation", "idx", "xend", "yend", 
  "label", "x", "z_clean"
))

#' Read Empirical Data from File
#'
#' Reads time series or panel data from various file formats and prepares
#' it for use with EmpiricalDynamics functions.
#'
#' @param file Path to the data file (CSV, RDS, or RData).
#' @param time_col Name of the time/date column (auto-detected if NULL).
#' @param date_format Date format string if time column is character.
#' @param ... Additional arguments passed to read.csv.
#'
#' @return A data.frame with time column converted to numeric if needed.
#'
#' @export
read_empirical_data <- function(file, time_col = NULL, date_format = NULL, ...) {
  
  # Determine file type and read
  ext <- tolower(tools::file_ext(file))
  
  data <- switch(ext,
                 "csv" = utils::read.csv(file, stringsAsFactors = FALSE, ...),
                 "rds" = readRDS(file),
                 "rdata" = {
                   env <- new.env()
                   load(file, envir = env)
                   obj_names <- ls(env)
                   if (length(obj_names) == 1) {
                     get(obj_names[1], envir = env)
                   } else {
                     stop("RData file contains multiple objects. Please use RDS format.")
                   }
                 },
                 stop("Unsupported file format: ", ext)
  )
  
  # Auto-detect time column if not specified
  if (is.null(time_col)) {
    time_candidates <- c("time", "t", "date", "year", "period", "Time", "Date", "Year")
    time_col <- intersect(time_candidates, names(data))[1]
    if (!is.na(time_col)) {
      message("Auto-detected time column: ", time_col)
    }
  }
  
  # Convert date column to numeric if needed
  if (!is.null(time_col) && time_col %in% names(data)) {
    if (inherits(data[[time_col]], "Date")) {
      data[[time_col]] <- as.numeric(data[[time_col]])
    } else if (is.character(data[[time_col]])) {
      if (!is.null(date_format)) {
        data[[time_col]] <- as.numeric(as.Date(data[[time_col]], format = date_format))
      } else {
        # Try common formats
        for (fmt in c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%m/%d/%Y", "%Y")) {
          parsed <- tryCatch(
            as.Date(data[[time_col]], format = fmt),
            error = function(e) NULL
          )
          if (!is.null(parsed) && !all(is.na(parsed))) {
            data[[time_col]] <- as.numeric(parsed)
            break
          }
        }
      }
    }
  }
  
  # Add class for method dispatch
  class(data) <- c("empirical_data", class(data))
  attr(data, "time_col") <- time_col
  
  data
}


#' Calculate Coefficient Change (Internal)
#' @noRd
coefficient_change <- function(eq_new, eq_old) {
  coef_new <- stats::coef(eq_new)
  coef_old <- stats::coef(eq_old)
  
  if (is.null(coef_new) || is.null(coef_old)) {
    return(Inf)
  }
  
  common <- intersect(names(coef_new), names(coef_old))
  
  if (length(common) == 0) {
    return(Inf)
  }
  
  sqrt(mean((coef_new[common] - coef_old[common])^2))
}


#' Default ggplot2 Theme for EmpiricalDynamics
#'
#' A clean, publication-ready theme for all diagnostic plots.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme object.
#'
#' @export
ed_theme <- function(base_size = 11, base_family = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to use ed_theme.")
  }
  
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey70"),
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey90", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}


#' Safely Check Package Availability (Internal)
#' @noRd
check_package <- function(pkg, purpose = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(purpose)) {
      msg <- paste0(msg, " for ", purpose)
    }
    msg <- paste0(msg, ".")
    stop(msg, call. = FALSE)
  }
  TRUE
}


#' Validate Time Series Data (Internal)
#' @noRd
validate_timeseries <- function(x, name = "x", allow_na = FALSE, min_length = 10) {
  if (!is.numeric(x)) {
    stop(paste0("'", name, "' must be numeric"), call. = FALSE)
  }
  
  if (length(x) < min_length) {
    stop(paste0("'", name, "' must have at least ", min_length, " observations"), 
         call. = FALSE)
  }
  
  if (!allow_na && any(is.na(x))) {
    n_na <- sum(is.na(x))
    stop(paste0("'", name, "' contains ", n_na, " NA values. ",
                "Use na.omit() or set allow_na = TRUE"), call. = FALSE)
  }
  
  invisible(x)
}


#' Construct Integration Matrix for TVR (Internal)
#' @noRd
build_integration_matrix <- function(n, dt) {
  A <- matrix(0, n, n)
  for (i in 2:n) {
    A[i, 1:(i-1)] <- dt[1:(i-1)]
  }
  A
}


#' Construct Difference Matrix (Internal)
#' @noRd
build_difference_matrix <- function(n) {
  D <- matrix(0, n - 1, n)
  for (i in 1:(n - 1)) {
    D[i, i] <- -1
    D[i, i + 1] <- 1
  }
  D
}


#' Safe Division (Internal)
#' @noRd
safe_divide <- function(x, y, floor = 1e-10) {
  x / pmax(abs(y), floor) * sign(y + floor)
}


#' Extract Formula Variables (Internal)
#' @noRd
extract_formula_vars <- function(formula) {
  all.vars(formula)
}


#' Create Safe Evaluation Environment (Internal)
#' @noRd
create_eval_env <- function(data, params = list()) {
  env <- new.env(parent = baseenv())
  
  # Add mathematical functions
  env$exp <- base::exp
  env$log <- base::log
  env$sqrt <- base::sqrt
  env$sin <- base::sin
  env$cos <- base::cos
  env$tan <- base::tan
  env$abs <- base::abs
  env$sign <- base::sign
  env$inv <- function(x) 1/x
  env$square <- function(x) x^2
  env$logistic <- function(x, k = 1, x0 = 0) 1 / (1 + base::exp(-k * (x - x0)))
  env$threshold <- function(x, c) ifelse(x > c, 1, 0)
  env$softplus <- function(x) base::log(1 + base::exp(x))
  env$relu <- function(x) base::pmax(0, x)
  
  # Add data
  for (nm in names(data)) {
    env[[nm]] <- data[[nm]]
  }
  
  # Add parameters
  for (nm in names(params)) {
    env[[nm]] <- params[[nm]]
  }
  
  env
}


#' Format Number (Internal)
#' @noRd
format_number <- function(x, digits = 4, scientific = TRUE) {
  if (is.na(x)) return("NA")
  if (scientific && (abs(x) > 1e4 || (abs(x) < 1e-3 && x != 0))) {
    sprintf(paste0("%.", digits, "e"), x)
  } else {
    sprintf(paste0("%.", digits, "f"), x)
  }
}


#' Generate Block Bootstrap Indices (Internal)
#' @noRd
block_bootstrap_indices <- function(n, block_size = NULL, n_samples = 200) {
  if (is.null(block_size)) {
    block_size <- max(5, floor(n^(1/3)))
  }
  n_blocks <- ceiling(n / block_size)
  max_start <- n - block_size + 1
  
  lapply(1:n_samples, function(i) {
    starts <- sample(1:max_start, n_blocks, replace = TRUE)
    indices <- unlist(lapply(starts, function(s) s:(s + block_size - 1)))
    indices[1:n]
  })
}


#' Print Rule (Internal)
#' @noRd
print_rule <- function(char = "=", width = 60) {
  cat(paste(rep(char, width), collapse = ""), "\n")
}


#' Safe Namespace Call (Internal)
#' @noRd
ns_call <- function(pkg, fn, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg, "' required but not installed."), call. = FALSE)
  }
  getExportedValue(pkg, fn)(...)
}


#' Load Example Dataset
#'
#' Load one of the example datasets included with the package.
#'
#' @param name Name of the dataset to load. Available datasets:
#'   \itemize{
#'     \item "logistic_growth" - Logistic population growth
#'     \item "predator_prey" - Lotka-Volterra predator-prey dynamics
#'     \item "interest_rate" - Vasicek mean-reverting interest rate
#'     \item "epidemic_data" - SIR epidemic model
#'     \item "oscillator_data" - Van der Pol oscillator
#'     \item "business_cycle" - Kaldor-type business cycle
#'   }
#'
#' @return A data.frame containing the time series data.
#'
#' @examples
#' \donttest{
#' # Load logistic growth data if available
#' if(requireNamespace("utils", quietly = TRUE)) {
#'   try({
#'     data <- load_example_data("logistic_growth")
#'     head(data)
#'   })
#' }
#' }
#'
#' @export
load_example_data <- function(name) {
  available <- c("logistic_growth", "predator_prey", "interest_rate",
                 "epidemic_data", "oscillator_data", "business_cycle")
  
  if (!name %in% available) {
    stop("Unknown dataset: '", name, "'. Available datasets: ",
         paste(available, collapse = ", "))
  }
  
  file <- system.file("extdata", paste0(name, ".csv"), 
                      package = "EmpiricalDynamics")
  
  if (file == "") {
    stop("Dataset file not found. Package may not be properly installed or data is missing.")
  }
  
  utils::read.csv(file, stringsAsFactors = FALSE)
}


#' List Available Example Datasets
#'
#' Returns information about the example datasets included with the package.
#'
#' @return A data.frame with dataset names and descriptions.
#'
#' @examples
#' list_example_data()
#'
#' @export
list_example_data <- function() {
  data.frame(
    name = c("logistic_growth", "predator_prey", "interest_rate",
             "epidemic_data", "oscillator_data", "business_cycle"),
    description = c(
      "Logistic population growth: dZ/dt = r*Z*(1-Z/K)",
      "Lotka-Volterra predator-prey: dX/dt = aX - bXY, dY/dt = cXY - dY",
      "Vasicek interest rate: dr = k(theta-r)dt + sigma*dW",
      "SIR epidemic model: dS/dt = -beta*S*I/N, dI/dt = beta*S*I/N - gamma*I",
      "Van der Pol oscillator: dy/dt = mu*(1-x^2)*y - x",
      "Kaldor business cycle: nonlinear investment-savings dynamics"
    ),
    variables = c("time, Z", "time, X, Y", "time, rate", 
                  "time, S, I, R", "time, x, y", "time, Y, K"),
    stringsAsFactors = FALSE
  )
}