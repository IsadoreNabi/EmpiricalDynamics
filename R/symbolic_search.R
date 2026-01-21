#' @title Symbolic Regression and Equation Discovery
#' @description Functions for discovering functional forms through symbolic
#'   regression using genetic algorithms. Interfaces with Julia's 
#'   SymbolicRegression.jl for advanced search, with fallback to R-native
#'   methods for simpler cases.
#' @name symbolic_search
NULL


#' Symbolic Equation Discovery
#'
#' Discovers the functional form of a differential equation from data using
#' genetic/evolutionary algorithms. Returns a Pareto front of equations
#' trading off complexity against fit.
#'
#' @param target Numeric vector of target values (typically derivatives).
#' @param predictors Data frame of predictor variables.
#' @param operators List specifying allowed operators:
#'   \itemize{
#'     \item binary: c("+", "-", "*", "/")
#'     \item unary: c("exp", "log", "sqrt", "inv", "square")
#'     \item custom: Custom function names (must be defined)
#'   }
#' @param constraints List of constraints:
#'   \itemize{
#'     \item forced: Formula of terms that must appear
#'     \item forbidden: Formula of terms that must not appear
#'     \item max_complexity: Maximum expression complexity
#'   }
#' @param n_runs Number of independent runs for robustness.
#' @param complexity_penalty Penalty per unit complexity.
#' @param parsimony_pressure Type of parsimony: "constant", "adaptive", or "none".
#' @param backend Computation backend: "julia", "r_genetic", or "r_exhaustive".
#' @param julia_options List of options passed to SymbolicRegression.jl.
#' @param weights Optional weight vector for weighted regression.
#' @param verbose Print progress messages?
#'
#' @return An object of class "symbolic_search_result" containing:
#'   \itemize{
#'     \item pareto_front: Data frame of Pareto-optimal equations
#'     \item all_equations: All discovered equations
#'     \item best_by_complexity: Best equation at each complexity level
#'     \item run_diagnostics: Information about each run
#'   }
#'
#' @examples
#' \donttest{
#' # Toy example using R-native exhaustive search (fastest for demo)
#' data <- data.frame(
#'   x = seq(1, 10, length.out = 20),
#'   y = seq(1, 10, length.out = 20)^2 + rnorm(20, sd = 0.1)
#' )
#' 
#' # Discover y ~ x^2
#' results <- symbolic_search(
#'   target = data$y,
#'   predictors = data["x"],
#'   backend = "r_exhaustive"
#' )
#' 
#' print(head(results$pareto_front))
#' }
#'
#' @export
symbolic_search <- function(target, predictors,
                            operators = NULL,
                            constraints = NULL,
                            n_runs = 5,
                            complexity_penalty = 0.05,
                            parsimony_pressure = c("adaptive", "constant", "none"),
                            backend = c("r_genetic", "julia", "r_exhaustive"),
                            julia_options = NULL,
                            weights = NULL,
                            verbose = TRUE) {
  
  backend <- match.arg(backend)
  parsimony_pressure <- match.arg(parsimony_pressure)
  
  # Validate inputs
  if (any(is.na(target))) stop("Target contains NA values.")
  
  if (!is.data.frame(predictors)) {
    predictors <- as.data.frame(predictors)
  }
  
  if (nrow(predictors) != length(target)) {
    stop("Number of rows in 'predictors' must equal length of 'target'", 
         call. = FALSE)
  }
  
  # Default operators
  if (is.null(operators)) {
    operators <- list(
      binary = c("+", "-", "*", "/"),
      unary = c("inv", "square", "sqrt", "log", "exp"),
      custom = c()
    )
  }
  
  # Default constraints
  if (is.null(constraints)) {
    constraints <- list(
      forced = NULL,
      forbidden = NULL,
      max_complexity = 30
    )
  }
  
  # Set up weights
  if (is.null(weights)) {
    weights <- rep(1, length(target))
  }
  
  # Dispatch to backend
  result <- switch(backend,
                   "julia" = symbolic_search_julia(target, predictors, operators, 
                                                   constraints, n_runs, complexity_penalty,
                                                   parsimony_pressure, julia_options, 
                                                   weights, verbose),
                   "r_genetic" = symbolic_search_r_genetic(target, predictors, operators,
                                                           constraints, n_runs, 
                                                           complexity_penalty,
                                                           parsimony_pressure, weights,
                                                           verbose),
                   "r_exhaustive" = symbolic_search_r_exhaustive(target, predictors, 
                                                                 operators, constraints,
                                                                 weights, verbose)
  )
  
  class(result) <- "symbolic_search_result"
  result
}


#' R-Native Genetic Algorithm for Symbolic Regression
#'
#' A pure R implementation of symbolic regression using a genetic algorithm.
#'
#' @keywords internal
symbolic_search_r_genetic <- function(target, predictors, operators, 
                                      constraints, n_runs, complexity_penalty,
                                      parsimony_pressure, weights, verbose) {
  
  n <- length(target)
  var_names <- names(predictors)
  n_vars <- length(var_names)
  
  # Expression node types
  NODE_CONST <- 0
  NODE_VAR <- 1
  NODE_UNARY <- 2
  NODE_BINARY <- 3
  
  # --- Helper Functions ---
  
  # Create initial population
  create_random_expr <- function(max_depth = 4, depth = 0) {
    if (depth >= max_depth || (depth > 0 && stats::runif(1) < 0.3)) {
      # Terminal node
      if (stats::runif(1) < 0.5) {
        # Constant
        list(type = NODE_CONST, value = stats::runif(1, -5, 5))
      } else {
        # Variable
        list(type = NODE_VAR, var = sample(var_names, 1))
      }
    } else if (stats::runif(1) < 0.3 && length(operators$unary) > 0) {
      # Unary operator
      list(
        type = NODE_UNARY,
        op = sample(operators$unary, 1),
        child = create_random_expr(max_depth, depth + 1)
      )
    } else {
      # Binary operator
      list(
        type = NODE_BINARY,
        op = sample(operators$binary, 1),
        left = create_random_expr(max_depth, depth + 1),
        right = create_random_expr(max_depth, depth + 1)
      )
    }
  }
  
  # Evaluate expression
  eval_expr <- function(expr, data) {
    tryCatch({
      if (expr$type == NODE_CONST) {
        rep(expr$value, nrow(data))
      } else if (expr$type == NODE_VAR) {
        data[[expr$var]]
      } else if (expr$type == NODE_UNARY) {
        child_val <- eval_expr(expr$child, data)
        switch(expr$op,
               "inv" = 1 / child_val,
               "square" = child_val^2,
               "sqrt" = sqrt(abs(child_val)),
               "log" = log(abs(child_val) + 1e-10),
               "exp" = exp(pmin(child_val, 20)),
               "abs" = abs(child_val),
               "neg" = -child_val,
               child_val
        )
      } else {
        left_val <- eval_expr(expr$left, data)
        right_val <- eval_expr(expr$right, data)
        switch(expr$op,
               "+" = left_val + right_val,
               "-" = left_val - right_val,
               "*" = left_val * right_val,
               "/" = left_val / (right_val + sign(right_val) * 1e-10),
               left_val + right_val
        )
      }
    }, error = function(e) rep(NA, nrow(data)))
  }
  
  # Convert expression to string
  expr_to_string <- function(expr) {
    if (expr$type == NODE_CONST) {
      sprintf("%.4f", expr$value)
    } else if (expr$type == NODE_VAR) {
      expr$var
    } else if (expr$type == NODE_UNARY) {
      sprintf("%s(%s)", expr$op, expr_to_string(expr$child))
    } else {
      sprintf("(%s %s %s)", expr_to_string(expr$left), 
              expr$op, expr_to_string(expr$right))
    }
  }
  
  # Calculate complexity
  expr_complexity <- function(expr) {
    if (expr$type %in% c(NODE_CONST, NODE_VAR)) {
      1
    } else if (expr$type == NODE_UNARY) {
      1 + expr_complexity(expr$child)
    } else {
      1 + expr_complexity(expr$left) + expr_complexity(expr$right)
    }
  }
  
  # Fitness function (minimize)
  fitness <- function(expr, data, target, weights, complexity_pen) {
    pred <- eval_expr(expr, data)
    if (any(is.na(pred)) || any(!is.finite(pred))) return(Inf)
    
    mse <- sum(weights * (target - pred)^2) / sum(weights)
    complexity <- expr_complexity(expr)
    
    # Check max complexity constraint
    if (!is.null(constraints$max_complexity) && complexity > constraints$max_complexity) {
      return(Inf)
    }
    
    mse + complexity_pen * complexity
  }
  
  # Mutation
  mutate <- function(expr, mutation_rate = 0.2) {
    if (stats::runif(1) < mutation_rate) {
      # Replace with new random subtree
      create_random_expr(max_depth = 3)
    } else if (expr$type == NODE_CONST) {
      expr$value <- expr$value + stats::rnorm(1, sd = 0.5)
      expr
    } else if (expr$type == NODE_VAR) {
      if (stats::runif(1) < 0.3) {
        expr$var <- sample(var_names, 1)
      }
      expr
    } else if (expr$type == NODE_UNARY) {
      expr$child <- mutate(expr$child, mutation_rate * 0.8)
      expr
    } else {
      expr$left <- mutate(expr$left, mutation_rate * 0.8)
      expr$right <- mutate(expr$right, mutation_rate * 0.8)
      expr
    }
  }
  
  # Crossover
  crossover <- function(expr1, expr2) {
    # Simplified: replace random subtree if binary, else return mutated parent
    if (expr1$type == NODE_BINARY && stats::runif(1) < 0.7) {
      if (stats::runif(1) < 0.5) {
        expr1$left <- expr2
      } else {
        expr1$right <- expr2
      }
      return(expr1)
    }
    # Fallback to mutation if crossover not applicable
    mutate(expr1, 0.1)
  }
  
  # --- GA Main Loop ---
  
  run_ga <- function(run_id, pop_size = 100, n_gen = 50) {
    if (verbose) message(sprintf("  Run %d/%d", run_id, n_runs))
    
    # Initialize population
    population <- lapply(1:pop_size, function(i) create_random_expr(max_depth = 5))
    
    # Current complexity penalty (for adaptive parsimony)
    current_penalty <- complexity_penalty
    
    best_by_gen <- list()
    
    for (gen in 1:n_gen) {
      # Evaluate fitness
      fits <- sapply(population, fitness, data = predictors, 
                     target = target, weights = weights, 
                     complexity_pen = current_penalty)
      
      # Handle infinities
      fits[!is.finite(fits)] <- max(fits[is.finite(fits)], 1e10)
      
      # Adaptive parsimony update
      if (parsimony_pressure == "adaptive" && gen %% 10 == 0) {
        complexities <- sapply(population, expr_complexity)
        mean_complexity <- mean(complexities)
        if (mean_complexity > 15) {
          current_penalty <- current_penalty * 1.5
        } else if (mean_complexity < 8) {
          current_penalty <- current_penalty * 0.8
        }
      }
      
      # Selection (tournament)
      selected <- numeric(pop_size)
      for (i in 1:pop_size) {
        candidates <- sample(pop_size, 3, replace = TRUE)
        selected[i] <- candidates[which.min(fits[candidates])]
      }
      
      # Create new population
      new_pop <- list()
      for (i in seq(1, pop_size, by = 2)) {
        parent1 <- population[[selected[i]]]
        # Ensure we have a second parent
        idx2 <- if (i + 1 <= pop_size) selected[i + 1] else selected[1]
        parent2 <- population[[idx2]]
        
        # Crossover
        if (stats::runif(1) < 0.7) {
          child1 <- crossover(parent1, parent2)
          child2 <- crossover(parent2, parent1)
        } else {
          child1 <- parent1
          child2 <- parent2
        }
        
        # Mutation
        child1 <- mutate(child1, 0.15)
        child2 <- mutate(child2, 0.15)
        
        new_pop[[i]] <- child1
        if (i + 1 <= pop_size) new_pop[[i + 1]] <- child2
      }
      
      # Elitism: keep best solution found so far
      best_idx <- which.min(fits)
      new_pop[[1]] <- population[[best_idx]]
      
      population <- new_pop
      
      # Track best of generation
      best_expr <- population[[1]]
      # Recalculate MSE without penalty for reporting
      best_mse <- fitness(best_expr, predictors, target, weights, 0)
      
      best_by_gen[[gen]] <- list(
        expr = best_expr,
        mse = best_mse,
        complexity = expr_complexity(best_expr)
      )
    }
    
    # Return results from final population
    all_results <- lapply(population, function(expr) {
      pred <- eval_expr(expr, predictors)
      mse <- if (any(is.na(pred))) Inf else mean((target - pred)^2)
      list(
        expr = expr,
        string = expr_to_string(expr),
        mse = mse,
        rmse = sqrt(mse),
        complexity = expr_complexity(expr)
      )
    })
    
    # Filter valid results
    all_results <- all_results[sapply(all_results, function(x) is.finite(x$mse))]
    
    all_results
  }
  
  # Run multiple times
  if (verbose) message("Running genetic algorithm...")
  
  all_runs <- lapply(1:n_runs, run_ga)
  all_results <- do.call(c, all_runs)
  
  # Build Pareto front
  pareto_front <- build_pareto_front(all_results)
  
  # Find best at each complexity
  best_by_complexity <- list()
  for (res in all_results) {
    c <- res$complexity
    key <- as.character(c)
    if (is.null(best_by_complexity[[key]]) || 
        res$mse < best_by_complexity[[key]]$mse) {
      best_by_complexity[[key]] <- res
    }
  }
  
  list(
    pareto_front = pareto_front,
    all_equations = all_results,
    best_by_complexity = best_by_complexity,
    backend = "r_genetic",
    n_runs = n_runs,
    predictors = var_names
  )
}


#' Build Pareto Front from Results
#'
#' @keywords internal
build_pareto_front <- function(results) {
  if (length(results) == 0) return(data.frame())
  
  # Extract mse and complexity
  df <- data.frame(
    idx = 1:length(results),
    mse = sapply(results, function(x) x$mse),
    complexity = sapply(results, function(x) x$complexity),
    string = sapply(results, function(x) x$string)
  )
  
  # Order by complexity then MSE
  df <- df[order(df$complexity, df$mse), ]
  
  # Find Pareto-optimal points
  pareto_idx <- c()
  current_best_mse <- Inf
  
  for (i in 1:nrow(df)) {
    if (df$mse[i] < current_best_mse) {
      pareto_idx <- c(pareto_idx, df$idx[i])
      current_best_mse <- df$mse[i]
    }
  }
  
  pareto_results <- results[pareto_idx]
  
  data.frame(
    complexity = sapply(pareto_results, function(x) x$complexity),
    mse = sapply(pareto_results, function(x) x$mse),
    rmse = sapply(pareto_results, function(x) x$rmse),
    equation = sapply(pareto_results, function(x) x$string),
    stringsAsFactors = FALSE
  )
}


#' Exhaustive Search for Simple Equations
#'
#' @keywords internal
symbolic_search_r_exhaustive <- function(target, predictors, operators,
                                         constraints, weights, verbose) {
  
  var_names <- names(predictors)
  n <- length(target)
  
  candidates <- list()
  
  # Level 1: Single variables
  for (v in var_names) {
    candidates[[length(candidates) + 1]] <- list(
      formula = stats::as.formula(paste("target ~", v)),
      string = v,
      complexity = 1
    )
  }
  
  # Level 2: Simple transformations
  for (v in var_names) {
    if ("inv" %in% operators$unary) {
      candidates[[length(candidates) + 1]] <- list(
        formula = stats::as.formula(paste("target ~ I(1/", v, ")")),
        string = paste0("1/", v),
        complexity = 2
      )
    }
    if ("square" %in% operators$unary) {
      candidates[[length(candidates) + 1]] <- list(
        formula = stats::as.formula(paste("target ~ I(", v, "^2)")),
        string = paste0(v, "^2"),
        complexity = 2
      )
    }
    if ("log" %in% operators$unary) {
      candidates[[length(candidates) + 1]] <- list(
        formula = stats::as.formula(paste("target ~ log(abs(", v, ") + 1e-10)")),
        string = paste0("log(", v, ")"),
        complexity = 2
      )
    }
  }
  
  # Level 3: Pairwise combinations
  if (length(var_names) >= 2) {
    for (i in 1:(length(var_names) - 1)) {
      for (j in (i + 1):length(var_names)) {
        v1 <- var_names[i]
        v2 <- var_names[j]
        
        # Sum
        candidates[[length(candidates) + 1]] <- list(
          formula = stats::as.formula(paste("target ~", v1, "+", v2)),
          string = paste(v1, "+", v2),
          complexity = 3
        )
        
        # Product
        candidates[[length(candidates) + 1]] <- list(
          formula = stats::as.formula(paste("target ~", v1, "*", v2)),
          string = paste(v1, "*", v2),
          complexity = 3
        )
      }
    }
  }
  
  # Fit all candidates
  data_fit <- cbind(predictors, target = target)
  
  results <- lapply(candidates, function(cand) {
    tryCatch({
      fit <- stats::lm(cand$formula, data = data_fit, weights = weights)
      pred <- stats::fitted(fit)
      mse <- sum(weights * (target - pred)^2) / sum(weights)
      
      list(
        expr = NULL,
        string = cand$string,
        formula = cand$formula,
        mse = mse,
        rmse = sqrt(mse),
        complexity = cand$complexity,
        fit = fit
      )
    }, error = function(e) NULL)
  })
  
  results <- results[!sapply(results, is.null)]
  
  pareto_front <- build_pareto_front(results)
  
  list(
    pareto_front = pareto_front,
    all_equations = results,
    best_by_complexity = list(),
    backend = "r_exhaustive",
    predictors = var_names
  )
}


#' Julia Backend for Symbolic Search
#'
#' Uses SymbolicRegression.jl for advanced symbolic regression.
#'
#' @keywords internal
symbolic_search_julia <- function(target, predictors, operators, constraints,
                                  n_runs, complexity_penalty, parsimony_pressure,
                                  julia_options, weights, verbose) {
  
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    warning("JuliaCall not available. Falling back to R backend.")
    return(symbolic_search_r_genetic(target, predictors, operators, constraints,
                                     n_runs, complexity_penalty, parsimony_pressure,
                                     weights, verbose))
  }
  
  # Initialize Julia
  tryCatch({
    JuliaCall::julia_setup()
    JuliaCall::julia_library("SymbolicRegression")
  }, error = function(e) {
    warning("Julia setup failed: ", e$message, "\nFalling back to R backend.")
    return(symbolic_search_r_genetic(target, predictors, operators, constraints,
                                     n_runs, complexity_penalty, parsimony_pressure,
                                     weights, verbose))
  })
  
  # Transfer data to Julia
  X <- as.matrix(predictors)
  y <- as.numeric(target)
  
  JuliaCall::julia_assign("X", t(X))  # Julia uses column-major
  JuliaCall::julia_assign("y", y)
  
  # Build options
  options_code <- sprintf("
    options = SymbolicRegression.Options(
        binary_operators = [+, -, *, /],
        unary_operators = [inv, square, sqrt, log, exp],
        complexity_of_operators = Dict(inv => 1, square => 1),
        parsimony = %f,
        maxsize = %d,
        niterations = %d
    )
  ", complexity_penalty, constraints$max_complexity %||% 30, n_runs * 100)
  
  JuliaCall::julia_command(options_code)
  
  # Run symbolic regression
  JuliaCall::julia_command("
    hall_of_fame = EquationSearch(X, y; options=options, niterations=100)
  ")
  
  # Extract results
  results <- JuliaCall::julia_eval("
    [(string(member.tree), member.loss, compute_complexity(member.tree)) 
     for member in hall_of_fame]
  ")
  
  # Convert to R format
  all_results <- lapply(results, function(r) {
    list(
      expr = NULL,
      string = r[[1]],
      mse = r[[2]],
      rmse = sqrt(r[[2]]),
      complexity = r[[3]]
    )
  })
  
  pareto_front <- build_pareto_front(all_results)
  
  list(
    pareto_front = pareto_front,
    all_equations = all_results,
    best_by_complexity = list(),
    backend = "julia",
    n_runs = n_runs,
    predictors = names(predictors)
  )
}


#' Setup Julia Backend
#'
#' Checks if Julia and the required SymbolicRegression.jl package are installed.
#'
#' @return Logical indicating if the backend is ready.
#'
#' @export
setup_julia_backend <- function() {
  
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    message("Package 'JuliaCall' is not installed. Please install it to use the Julia backend.")
    return(FALSE)
  }
  
  message("Checking Julia backend configuration...")
  
  tryCatch({
    JuliaCall::julia_setup()
    
    # Check if SymbolicRegression.jl is installed
    check_code <- '
      try
        import Pkg
        haskey(Pkg.dependencies(), Base.UUID("8254be44-1295-4e6a-a16d-e31fe2c4a48b"))
      catch
        false
      end
    '
    installed <- JuliaCall::julia_eval(check_code)
    
    if (!installed) {
      message("\n'SymbolicRegression.jl' is not installed in Julia.")
      message("To install it, please run the following in R:")
      message("  JuliaCall::julia_library('Pkg')")
      message("  JuliaCall::julia_command('Pkg.add(\"SymbolicRegression\")')")
      return(FALSE)
    }
    
    JuliaCall::julia_library("SymbolicRegression")
    message("Julia backend is ready!")
    TRUE
    
  }, error = function(e) {
    warning("Julia setup failed. Error: ", e$message)
    message("Please ensure Julia is installed and accessible.")
    FALSE
  })
}


#' Weighted Symbolic Search
#'
#' Performs symbolic search with weighted least squares.
#'
#' @inheritParams symbolic_search
#'
#' @return A symbolic_search_result object
#'
#' @export
symbolic_search_weighted <- function(target, predictors, weights,
                                     operators = NULL, constraints = NULL,
                                     n_runs = 3, complexity_penalty = 0.05,
                                     verbose = TRUE) {
  
  symbolic_search(
    target = target,
    predictors = predictors,
    operators = operators,
    constraints = constraints,
    n_runs = n_runs,
    complexity_penalty = complexity_penalty,
    backend = "r_genetic",
    weights = weights,
    verbose = verbose
  )
}


#' Fit Specified Equation
#'
#' Fits a researcher-specified functional form, estimating only the
#' parameters. Uses Levenberg-Marquardt algorithm for robustness.
#'
#' @param expression Character string specifying the equation (e.g., "a + b * Z").
#' @param data Data frame with predictor variables.
#' @param response Name of the response/target column.
#' @param derivative_col Alias for response (for compatibility).
#' @param start List of starting values for parameters (auto-estimated if NULL).
#' @param method Optimization method: "LM" (Levenberg-Marquardt, recommended),
#'   "nls" (standard), or "optim" (general optimization).
#' @param weights Optional weight vector.
#' @param lower Lower bounds for parameters (for "optim" method).
#' @param upper Upper bounds for parameters (for "optim" method).
#'
#' @return An object of class "symbolic_equation" containing the fitted model.
#'
#' @examples
#' \donttest{
#' # Toy example
#' data <- data.frame(Z = seq(1, 10, length.out = 20))
#' data$dZ <- 0.5 * data$Z * (1 - data$Z / 20) + rnorm(20, sd = 0.01)
#'
#' # Fit logistic equation
#' eq <- fit_specified_equation(
#'   expression = "r * Z * (1 - Z/K)",
#'   data = data,
#'   response = "dZ",
#'   start = list(r = 0.5, K = 20)
#' )
#' print(eq)
#' }
#'
#' @export
fit_specified_equation <- function(expression, data, response = NULL,
                                   derivative_col = NULL, start = NULL,
                                   method = c("LM", "nls", "optim"),
                                   weights = NULL, lower = -Inf, 
                                   upper = Inf) {
  
  # Compatibility: derivative_col is alias for response
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  if (is.null(response)) {
    stop("Must specify 'response' or 'derivative_col'", call. = FALSE)
  }
  
  # Validate that response exists in data
  if (!response %in% names(data)) {
    stop("Column '", response, "' not found in data", call. = FALSE)
  }
  
  method <- match.arg(method)
  
  # Auto-estimate starting values if not provided
  if (is.null(start)) {
    start <- estimate_initial_values(expression, data, response)
  }
  
  # Build nls formula from expression string
  nls_formula <- stats::as.formula(paste(response, "~", expression))
  
  # Define fitting functions
  try_lm <- function() {
    if (!requireNamespace("minpack.lm", quietly = TRUE)) {
      stop("Package 'minpack.lm' is required for LM optimization.")
    }
    # nls.lm.control DOES NOT ACCEPT warnOnly
    minpack.lm::nlsLM(
      formula = nls_formula,
      data = data,
      start = start,
      weights = weights,
      control = minpack.lm::nls.lm.control(maxiter = 500)
    )
  }
  
  try_nls <- function() {
    stats::nls(
      formula = nls_formula,
      data = data,
      start = start,
      weights = weights,
      control = stats::nls.control(maxiter = 200, warnOnly = TRUE)
    )
  }
  
  try_optim <- function() {
    fit_with_optim(response, expression, data, start, weights, lower, upper)
  }
  
  # Try fitting with fallback
  result <- tryCatch({
    if (method == "LM") {
      try_lm()
    } else if (method == "nls") {
      try_nls()
    } else {
      try_optim()
    }
  }, error = function(e) {
    # If primary method fails, try standard NLS
    message("Primary optimization method failed: ", e$message)
    message("Attempting fallback to standard NLS...")
    tryCatch({
      try_nls()
    }, error = function(e2) {
      # If standard NLS fails, try general optimization
      message("Standard NLS failed. Attempting general optimization (optim)...")
      try_optim()
    })
  })
  
  # Wrap in symbolic_equation class
  eq <- list(
    fit = result,
    expression = expression,
    formula = nls_formula,
    nls_formula = nls_formula,
    coefficients = stats::coef(result),
    string = format_equation_string(expression, stats::coef(result)),
    method = method,
    target = response,
    predictors = setdiff(all.vars(nls_formula), c(names(start), response))
  )
  
  class(eq) <- "symbolic_equation"
  eq
}


#' Fit Using General Optimization
#'
#' @keywords internal
fit_with_optim <- function(target_name, expression, data, start, 
                           weights, lower, upper) {
  
  if (is.null(weights)) weights <- rep(1, nrow(data))
  
  target_vals <- data[[target_name]]
  
  objective <- function(params) {
    names(params) <- names(start)
    env <- create_eval_env(data, as.list(params))
    pred <- eval(parse(text = expression), envir = env)
    if(any(!is.finite(pred))) return(Inf)
    sum(weights * (target_vals - pred)^2)
  }
  
  result <- stats::optim(
    par = unlist(start),
    fn = objective,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )
  
  names(result$par) <- names(start)
  
  structure(
    list(
      m = list(
        getPars = function() result$par,
        getAllPars = function() result$par
      ),
      convergence = result$convergence,
      message = result$message,
      coefficients = result$par
    ),
    class = "optim_fit"
  )
}


#' Automatic Initial Value Estimation
#'
#' Estimates reasonable starting values for nonlinear least squares.
#'
#' @param expression Character string with the equation expression.
#' @param data Data frame with variables.
#' @param response Name of the response variable.
#' @param derivative_col Alias for response (for compatibility).
#' @param method Estimation method: "grid_search", "random", or "heuristic".
#' @param n_tries Number of attempts for random method.
#'
#' @return Named list of initial values.
#'
#' @export
estimate_initial_values <- function(expression, data, response = NULL,
                                    derivative_col = NULL,
                                    method = c("grid_search", "random", "heuristic"),
                                    n_tries = 100) {
  
  # Compatibility
  if (is.null(response) && !is.null(derivative_col)) {
    response <- derivative_col
  }
  
  method <- match.arg(method)
  
  # Parse expression to find parameters
  expr_parsed <- parse(text = expression)
  formula_vars <- all.vars(expr_parsed)
  data_vars <- names(data)
  params <- setdiff(formula_vars, data_vars)
  
  if (length(params) == 0) {
    stop("No parameters to estimate in expression", call. = FALSE)
  }
  
  target <- data[[response]]
  
  eval_at <- function(param_values) {
    tryCatch({
      env <- create_eval_env(data, param_values)
      pred <- eval(parse(text = expression), envir = env)
      if (any(!is.finite(pred))) return(Inf)
      mean((target - pred)^2)
    }, error = function(e) Inf)
  }
  
  if (method == "grid_search") {
    # Coarse grid search
    grid_vals <- c(-10, -1, -0.1, 0.1, 1, 10)
    grid <- expand.grid(setNames(
      replicate(length(params), grid_vals, simplify = FALSE),
      params
    ))
    
    errors <- apply(grid, 1, function(row) {
      eval_at(as.list(row))
    })
    
    best_idx <- which.min(errors)
    if (length(best_idx) == 0 || !is.finite(errors[best_idx])) {
      return(setNames(as.list(rep(0.1, length(params))), params))
    }
    as.list(grid[best_idx, ])
    
  } else if (method == "random") {
    best_error <- Inf
    best_start <- NULL
    
    for (i in 1:n_tries) {
      start <- setNames(
        lapply(params, function(p) stats::runif(1, -10, 10)),
        params
      )
      error <- eval_at(start)
      
      if (error < best_error) {
        best_error <- error
        best_start <- start
      }
    }
    
    if (is.null(best_start)) {
      return(setNames(as.list(rep(0.1, length(params))), params))
    }
    best_start
    
  } else {
    # Simple heuristic: all 0.1
    setNames(as.list(rep(0.1, length(params))), params)
  }
}


#' Format Equation String
#'
#' @keywords internal
format_equation_string <- function(expression, coefficients) {
  expr <- expression
  
  for (name in names(coefficients)) {
    val <- sprintf("%.4f", coefficients[name])
    expr <- gsub(paste0("\\b", name, "\\b"), val, expr)
  }
  
  expr
}


#' Plot Pareto Front
#'
#' Visualizes the trade-off between equation complexity and fit quality.
#'
#' @param results A symbolic_search_result object.
#' @param highlight_selection Which equation to highlight ("knee", "BIC", "AIC", or index).
#' @param show_all Show all equations or just Pareto front?
#'
#' @return A ggplot object.
#'
#' @export
plot_pareto_front <- function(results, highlight_selection = "knee",
                              show_all = FALSE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  if (!inherits(results, "symbolic_search_result")) {
    stop("'results' must be a symbolic_search_result object", call. = FALSE)
  }
  
  pareto <- results$pareto_front
  
  p <- ggplot2::ggplot(pareto, 
                       ggplot2::aes(x = complexity, y = rmse, label = equation)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Pareto Front: Complexity vs. Fit",
      subtitle = "Lower RMSE = better fit, Lower complexity = simpler model",
      x = "Expression Complexity",
      y = "RMSE (log scale)"
    ) +
    ed_theme()
  
  if (show_all && length(results$all_equations) > 0) {
    all_df <- data.frame(
      complexity = sapply(results$all_equations, function(x) x$complexity),
      rmse = sapply(results$all_equations, function(x) x$rmse)
    )
    all_df <- all_df[is.finite(all_df$rmse), ]
    
    p <- p + ggplot2::geom_point(data = all_df, alpha = 0.2, size = 1)
  }
  
  # Highlight selected equation
  if (is.character(highlight_selection)) {
    sel_idx <- switch(highlight_selection,
                      "knee" = find_knee_point(pareto$complexity, pareto$rmse),
                      1  # Default to first
    )
  } else {
    sel_idx <- highlight_selection
  }
  
  if (sel_idx >= 1 && sel_idx <= nrow(pareto)) {
    p <- p + ggplot2::geom_point(
      data = pareto[sel_idx, ],
      color = "red", size = 5, shape = 1, stroke = 2
    )
  }
  
  p
}


#' Find Knee Point in Pareto Front
#'
#' Uses the maximum curvature method to find the "elbow" of the Pareto front.
#'
#' @keywords internal
find_knee_point <- function(x, y) {
  # Normalize to [0, 1]
  x_norm <- (x - min(x)) / (max(x) - min(x) + 1e-10)
  y_norm <- (y - min(y)) / (max(y) - min(y) + 1e-10)
  
  # Distance from line connecting first and last points
  n <- length(x)
  if (n < 2) return(1)
  
  # Line from (0, y_norm[1]) to (1, y_norm[n])
  a <- y_norm[n] - y_norm[1]
  b <- x_norm[1] - x_norm[n]
  c <- x_norm[n] * y_norm[1] - x_norm[1] * y_norm[n]
  
  distances <- abs(a * x_norm + b * y_norm + c) / sqrt(a^2 + b^2)
  
  which.max(distances)
}


#' Select Equation from Pareto Front
#'
#' Selects the best equation from the Pareto front using a specified criterion.
#'
#' @param results A symbolic_search_result object.
#' @param criterion Selection criterion: "knee", "BIC", "AIC", "min_complexity",
#'   or "min_error".
#' @param n Sample size (required for BIC/AIC if not stored).
#'
#' @return A symbolic_equation object.
#'
#' @export
select_equation <- function(results, criterion = c("knee", "BIC", "AIC",
                                                   "min_complexity", "min_error"),
                            n = NULL) {
  
  criterion <- match.arg(criterion)
  
  pareto <- results$pareto_front
  if (nrow(pareto) == 0) return(NULL)
  
  idx <- switch(criterion,
                "knee" = find_knee_point(pareto$complexity, pareto$rmse),
                "min_complexity" = 1,
                "min_error" = which.min(pareto$mse),
                "BIC" = {
                  if (is.null(n)) n <- 100  # Default
                  k <- pareto$complexity
                  bic <- n * log(pareto$mse) + k * log(n)
                  which.min(bic)
                },
                "AIC" = {
                  if (is.null(n)) n <- 100
                  k <- pareto$complexity
                  aic <- n * log(pareto$mse) + 2 * k
                  which.min(aic)
                }
  )
  
  selected <- if (length(results$all_equations) > 0) {
    # Find matching equation in all results
    matching <- which(sapply(results$all_equations, function(x) {
      x$complexity == pareto$complexity[idx] && 
        abs(x$mse - pareto$mse[idx]) < 1e-10
    }))
    if (length(matching) > 0) {
      results$all_equations[[matching[1]]]
    } else {
      list(
        string = pareto$equation[idx],
        mse = pareto$mse[idx],
        rmse = pareto$rmse[idx],
        complexity = pareto$complexity[idx]
      )
    }
  } else {
    list(
      string = pareto$equation[idx],
      mse = pareto$mse[idx],
      rmse = pareto$rmse[idx],
      complexity = pareto$complexity[idx]
    )
  }
  
  eq <- list(
    string = selected$string,
    expression = selected$string,
    expr = selected$expr,
    mse = selected$mse,
    rmse = selected$rmse,
    complexity = selected$complexity,
    selection_criterion = criterion,
    predictors = results$predictors
  )
  
  class(eq) <- "symbolic_equation"
  eq
}


#' Get Full Pareto Set
#'
#' Returns all equations on the Pareto front as a list.
#'
#' @param results A symbolic_search_result object.
#'
#' @return List of symbolic_equation objects.
#'
#' @export
get_pareto_set <- function(results) {
  
  pareto <- results$pareto_front
  
  lapply(1:nrow(pareto), function(i) {
    eq <- list(
      string = pareto$equation[i],
      expression = pareto$equation[i],
      mse = pareto$mse[i],
      rmse = pareto$rmse[i],
      complexity = pareto$complexity[i]
    )
    class(eq) <- "symbolic_equation"
    eq
  })
}


#' Define Custom Operators
#'
#' Defines custom mathematical operators for use in symbolic search.
#'
#' @param ... Named functions to add as operators.
#'
#' @return List of operator definitions suitable for symbolic_search.
#'
#' @examples
#' \donttest{
#' ops <- define_custom_operators(
#'   logistic = function(x, k = 1, x0 = 0) 1 / (1 + exp(-k * (x - x0))),
#'   threshold = function(x, c) ifelse(x > c, 1, 0)
#' )
#' }
#'
#' @export
define_custom_operators <- function(...) {
  funcs <- list(...)
  
  list(
    binary = c("+", "-", "*", "/"),
    unary = c("inv", "square", "sqrt", "log", "exp"),
    custom = names(funcs),
    custom_definitions = funcs
  )
}


# S3 Methods for symbolic_equation

#' @export
print.symbolic_equation <- function(x, ...) {
  cat("Symbolic Equation\n")
  cat("=================\n")
  cat("Expression:", x$string, "\n")
  if (!is.null(x$rmse)) cat(sprintf("RMSE: %.6f\n", x$rmse))
  if (!is.null(x$complexity)) cat(sprintf("Complexity: %d\n", x$complexity))
  if (!is.null(x$coefficients)) {
    cat("Coefficients:\n")
    print(x$coefficients)
  }
  invisible(x)
}


#' @export
coef.symbolic_equation <- function(object, ...) {
  if (!is.null(object$coefficients)) {
    object$coefficients
  } else if (!is.null(object$fit)) {
    stats::coef(object$fit)
  } else {
    NULL
  }
}


#' @export
predict.symbolic_equation <- function(object, newdata = NULL, ...) {
  if (!is.null(object$fit)) {
    stats::predict(object$fit, newdata = newdata, ...)
  } else if (!is.null(object$expression) && !is.null(object$coefficients)) {
    # Evaluate expression with coefficients
    env <- create_eval_env(newdata, as.list(object$coefficients))
    eval(parse(text = object$expression), envir = env)
  } else if (!is.null(object$expr)) {
    warning("Direct expression evaluation not yet implemented for prediction")
    NULL
  } else {
    warning("Cannot predict without fitted model or expression")
    NULL
  }
}


#' @export
summary.symbolic_equation <- function(object, ...) {
  cat("\nSymbolic Equation Summary\n")
  cat("=========================\n\n")
  
  cat("Formula:", object$string, "\n\n")
  
  if (!is.null(object$fit)) {
    cat("Model Summary:\n")
    print(summary(object$fit))
  }
  
  if (!is.null(object$rmse)) {
    cat(sprintf("\nGoodness of fit:\n  RMSE: %.6f\n  MSE: %.6f\n", 
                object$rmse, object$mse))
  }
  
  if (!is.null(object$complexity)) {
    cat(sprintf("  Complexity: %d\n", object$complexity))
  }
  
  invisible(object)
}