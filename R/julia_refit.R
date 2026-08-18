# =============================================================================
# EmpiricalDynamics: julia_refit.R
# Re-estimating the constants of a discovered equation through the same engine
# that discovered it.
# =============================================================================

#' @title Julia Refitting Engine
#' @description Support for re-estimating the constants of a discovered equation
#'   with SymbolicRegression.jl's own constant optimisation, which is the engine
#'   and the loss the search itself used.
#' @name julia_refit
#' @keywords internal
NULL

# Session-level state: setting Julia up and sourcing the backend costs seconds,
# and a cross-validation asks for a refit once per fold per equation.
.ed_julia_state <- new.env(parent = emptyenv())


#' Is the Julia refitting engine usable in this session?
#'
#' @return \code{TRUE} if Julia, JuliaCall and the shipped backend are all in
#'   place. Sets up the session the first time it is asked, and remembers.
#' @keywords internal
#' @noRd
ed_julia_refit_ready <- function() {
  if (!is.null(.ed_julia_state$ready)) return(.ed_julia_state$ready)

  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    .ed_julia_state$ready <- FALSE
    .ed_julia_state$reason <- "package 'JuliaCall' is not installed"
    return(FALSE)
  }

  backend_path <- system.file("julia", "symbolic_backend.jl",
                              package = "EmpiricalDynamics")
  if (!nzchar(backend_path)) {
    .ed_julia_state$ready <- FALSE
    .ed_julia_state$reason <- "the Julia backend shipped with the package could not be located"
    return(FALSE)
  }

  ok <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    JuliaCall::julia_library("SymbolicRegression")
    JuliaCall::julia_library("DynamicExpressions")
    JuliaCall::julia_source(backend_path)
    TRUE
  }, error = function(e) {
    .ed_julia_state$reason <- conditionMessage(e)
    FALSE
  })

  .ed_julia_state$ready <- isTRUE(ok)
  .ed_julia_state$ready
}


#' Why the Julia refitting engine is not usable
#' @keywords internal
#' @noRd
ed_julia_refit_reason <- function() {
  reason <- .ed_julia_state$reason
  if (is.null(reason)) "Julia could not be started" else reason
}


#' Refit an equation's constants in Julia and predict on held-out rows
#'
#' @param expression Character string with the discovered equation, constants
#'   included as literals: this engine optimises the constant nodes of the
#'   expression tree directly and needs no parameterisation.
#' @param train_data,test_data Data frames with the training and the held-out
#'   rows.
#' @param response Name of the response column.
#' @param weights Optional observation weights for the training rows.
#'
#' @return A list with \code{ok}, \code{message}, \code{predictions} (on
#'   \code{test_data}), \code{constants}, \code{n_constants} and
#'   \code{expression} (the refitted equation).
#' @keywords internal
#' @noRd
ed_julia_refit <- function(expression, train_data, test_data, response,
                           weights = NULL) {

  fail <- function(msg) {
    list(ok = FALSE, message = msg, predictions = numeric(0),
         constants = numeric(0), constants_initial = numeric(0),
         n_constants = 0L, expression = expression, moved = FALSE,
         loss_before = NA_real_, loss_after = NA_real_)
  }

  if (!ed_julia_refit_ready()) {
    return(fail(paste0("the Julia refitting engine is unavailable: ",
                       ed_julia_refit_reason())))
  }

  parsed_vars <- tryCatch(all.vars(parse(text = expression)),
                          error = function(e) NULL)
  if (is.null(parsed_vars)) {
    return(fail("the equation could not be parsed"))
  }

  unknown <- setdiff(parsed_vars, names(train_data))
  if (length(unknown) > 0L) {
    return(fail(paste0("the equation mentions variables that are not in the data: ",
                       paste(unknown, collapse = ", "))))
  }

  vars <- intersect(parsed_vars, names(train_data))
  if (length(vars) == 0L) {
    # An equation with no variables is still refittable -- its constant moves --
    # but the search engine wants a feature matrix, so one unused column is
    # carried along to give it a shape.
    numeric_cols <- names(train_data)[vapply(train_data, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, response)
    if (length(numeric_cols) == 0L) {
      return(fail("no numeric predictor is available to shape the data for Julia"))
    }
    vars <- numeric_cols[1]
  }

  x_train <- unname(t(as.matrix(train_data[, vars, drop = FALSE])))
  x_test <- unname(t(as.matrix(test_data[, vars, drop = FALSE])))
  y_train <- as.numeric(train_data[[response]])

  if (!all(is.finite(x_train)) || !all(is.finite(y_train))) {
    return(fail("the training rows carry non-finite values"))
  }

  transferred <- tryCatch({
    JuliaCall::julia_assign("_ed_rf_expr", expression)
    JuliaCall::julia_assign("_ed_rf_xtr", x_train)
    JuliaCall::julia_assign("_ed_rf_xte", x_test)
    JuliaCall::julia_assign("_ed_rf_ytr", y_train)
    JuliaCall::julia_assign("_ed_rf_vars", vars)
    # A single variable is transferred as a scalar String, and the search
    # engine requires a Vector{String}.
    JuliaCall::julia_command(
      "_ed_rf_vars = _ed_rf_vars isa AbstractVector ? String.(_ed_rf_vars) : String[String(_ed_rf_vars)]; nothing")
    # A single row or a single variable arrives without its matrix shape.
    JuliaCall::julia_command(
      "_ed_rf_xtr = _ed_rf_xtr isa AbstractMatrix ? Matrix{Float64}(_ed_rf_xtr) : reshape(Float64.(_ed_rf_xtr), length(_ed_rf_vars), :); nothing")
    JuliaCall::julia_command(
      "_ed_rf_xte = _ed_rf_xte isa AbstractMatrix ? Matrix{Float64}(_ed_rf_xte) : reshape(Float64.(_ed_rf_xte), length(_ed_rf_vars), :); nothing")
    if (is.null(weights)) {
      JuliaCall::julia_command("_ed_rf_w = nothing; nothing")
    } else {
      JuliaCall::julia_assign("_ed_rf_w", as.numeric(weights))
    }
    TRUE
  }, error = function(e) {
    .ed_julia_state$last_error <- conditionMessage(e)
    FALSE
  })

  if (!isTRUE(transferred)) {
    return(fail(paste0("the data could not be transferred to Julia: ",
                       .ed_julia_state$last_error)))
  }

  called <- tryCatch({
    JuliaCall::julia_command(paste0(
      "_ed_rf_out = ed_refit_constants(_ed_rf_expr, _ed_rf_xtr, _ed_rf_ytr, ",
      "_ed_rf_xte, _ed_rf_vars; weights = _ed_rf_w); nothing"))
    TRUE
  }, error = function(e) {
    .ed_julia_state$last_error <- conditionMessage(e)
    FALSE
  })

  if (!isTRUE(called)) {
    return(fail(paste0("the Julia refit could not be run: ",
                       .ed_julia_state$last_error)))
  }

  out <- tryCatch({
    list(
      ok = as.logical(JuliaCall::julia_eval("_ed_rf_out.ok")),
      message = as.character(JuliaCall::julia_eval("_ed_rf_out.message")),
      predictions = as.numeric(JuliaCall::julia_eval("_ed_rf_out.predictions")),
      constants = as.numeric(JuliaCall::julia_eval("_ed_rf_out.constants")),
      constants_initial = as.numeric(JuliaCall::julia_eval("_ed_rf_out.constants_initial")),
      n_constants = as.integer(JuliaCall::julia_eval("_ed_rf_out.n_constants")),
      expression = as.character(JuliaCall::julia_eval("_ed_rf_out.expression")),
      moved = as.logical(JuliaCall::julia_eval("_ed_rf_out.moved")),
      loss_before = as.numeric(JuliaCall::julia_eval("_ed_rf_out.loss_before")),
      loss_after = as.numeric(JuliaCall::julia_eval("_ed_rf_out.loss_after"))
    )
  }, error = function(e) NULL)

  if (is.null(out)) {
    return(fail("the result of the Julia refit could not be read back"))
  }
  if (!isTRUE(out$ok)) {
    return(fail(if (nzchar(out$message)) out$message else "the Julia refit failed"))
  }
  if (length(out$predictions) != nrow(test_data)) {
    return(fail(sprintf(
      "the Julia refit returned %d predictions for %d held-out rows",
      length(out$predictions), nrow(test_data))))
  }

  out
}
