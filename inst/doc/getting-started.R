## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
# # Install from local source
# install.packages("EmpiricalDynamics", repos = NULL, type = "source")
# 
# # Or using devtools
# devtools::install_github("your-repo/EmpiricalDynamics")

## ----quickstart-prep, eval=FALSE----------------------------------------------
# library(EmpiricalDynamics)
# 
# # Load example data
# logistic_growth <- load_example_data("logistic_growth")
# 
# # View the data
# head(logistic_growth)
# 
# # Compute derivative using Total Variation Regularization (recommended)
# deriv <- compute_derivative(
#   logistic_growth$Z,
#   time = logistic_growth$time,
#   method = "tvr",
#   lambda = "auto"
# )
# 
# # Add derivative to data
# logistic_growth$dZ <- deriv$derivative
# 
# # Plot diagnostic
# plot(deriv)

## ----quickstart-explore, eval=FALSE-------------------------------------------
# # Explore relationships between dZ/dt and Z
# exploration <- explore_dynamics(
#   data = logistic_growth,
#   response = "dZ",
#   predictors = "Z"
# )
# 
# # View suggestions for functional forms
# print(exploration$suggestions)

## ----quickstart-fit, eval=FALSE-----------------------------------------------
# # Fit the theoretical equation
# equation <- fit_specified_equation(
#   "r * Z * (1 - Z / K)",
#   data = logistic_growth,
#   derivative_col = "dZ",
#   start = list(r = 0.5, K = 100)
# )
# 
# # View results
# summary(equation)
# print(coef(equation))

## ----quickstart-validate, eval=FALSE------------------------------------------
# # Run comprehensive validation
# validation <- validate_model(
#   equation = equation,
#   data = logistic_growth,
#   derivative_col = "dZ",
#   variable = "Z",
#   cv_folds = 5
# )
# 
# print(validation)

## ----quickstart-output, eval=FALSE--------------------------------------------
# # Get LaTeX equation
# latex_eq <- to_latex(equation, variable = "Z")
# cat(latex_eq)
# 
# # Generate full report (using tempdir() per CRAN policy)
# generate_report(
#   results = list(
#     data = logistic_growth,
#     equation = equation,
#     validation = validation
#   ),
#   output_file = file.path(tempdir(), "logistic_analysis"),
#   format = "markdown"
# )

## ----tvr-example, eval=FALSE--------------------------------------------------
# deriv_tvr <- compute_derivative(
#   y = data$Z,
#   time = data$time,
#   method = "tvr",
#   lambda = "auto"  # Automatic regularization selection
# )
# 
# # Plot diagnostic: shows original, derivative, reconstruction, residuals
# plot_tvr_diagnostic(deriv_tvr)

## ----other-methods, eval=FALSE------------------------------------------------
# # Savitzky-Golay filter (preserves peaks, good for smooth data)
# deriv_sg <- compute_derivative(y, time, method = "savgol",
#                                 window = 11, order = 3)
# 
# # Smoothing spline (good for high noise)
# deriv_spline <- compute_derivative(y, time, method = "spline",
#                                     spar = 0.7)
# 
# # Finite differences (simple, fast, sensitive to noise)
# deriv_fd <- compute_derivative(y, time, method = "fd")
# 
# # Spectral (for periodic data ONLY - warns about Gibbs phenomenon otherwise)
# deriv_spec <- compute_derivative(y, time, method = "spectral")

## ----suggest-method, eval=FALSE-----------------------------------------------
# # Get recommendation based on data characteristics
# suggestion <- suggest_differentiation_method(data$Z, data$time)
# print(suggestion)
# 
# # Compare all methods visually
# compare_differentiation_methods(data$Z, data$time)

## ----exploration, eval=FALSE--------------------------------------------------
# # Comprehensive exploration
# exploration <- explore_dynamics(
#   data = data,
#   response = "dZ",
#   predictors = c("Z", "X", "Y")
# )
# 
# # Individual plots
# plot_phase_1d(data, x = "Z", y = "dZ")
# plot_bivariate(data, x = "Z", y = "dZ")
# plot_surface_3d(data, x1 = "Z", x2 = "X", y = "dZ")

## ----symbolic-search, eval=FALSE----------------------------------------------
# # Genetic algorithm search (pure R)
# search_result <- symbolic_search(
#   data = data,
#   response = "dZ",
#   predictors = c("Z", "X"),
#   backend = "r_genetic",
#   max_complexity = 15,
#   n_generations = 50,
#   population_size = 100,
#   n_runs = 3
# )
# 
# # View Pareto front
# plot_pareto_front(search_result)
# 
# # Select best equation
# best_eq <- select_equation(search_result, criterion = "bic")

## ----theory-guided, eval=FALSE------------------------------------------------
# # Fit a specific functional form
# equation <- fit_specified_equation(
#   "alpha + beta * Z + gamma * Z^2 + delta * X",
#   data = data,
#   derivative_col = "dZ",
#   method = "levenberg-marquardt"  # Robust optimization
# )
# 
# summary(equation)

## ----julia-backend, eval=FALSE------------------------------------------------
# # Setup Julia (one-time)
# setup_julia_backend()
# 
# # Run with Julia's SymbolicRegression.jl
# search_result <- symbolic_search(
#   data = data,
#   response = "dZ",
#   predictors = c("Z", "X"),
#   backend = "julia",
#   max_complexity = 25
# )

## ----residual-diag, eval=FALSE------------------------------------------------
# # Compute residuals
# resid <- compute_residuals(equation, data, "dZ")
# 
# # Run diagnostic tests
# diag <- residual_diagnostics(resid, data)
# print(diag)
# plot(diag)

## ----diffusion, eval=FALSE----------------------------------------------------
# # Model conditional variance
# var_model <- model_conditional_variance(
#   resid,
#   data = data,
#   predictors = c("Z"),
#   method = "linear"  # or "quadratic", "symbolic"
# )
# 
# # Construct SDE: dZ = f(Z,X) dt + g(Z) dW
# sde <- construct_sde(
#   drift = equation,
#   diffusion = var_model,
#   variable = "Z"
# )
# 
# print(sde)

## ----iterative-gls, eval=FALSE------------------------------------------------
# # Iterative estimation for unbiased drift coefficients
# sde <- estimate_sde_iterative(
#   drift_formula = "a + b * Z + c * Z^2",
#   data = data,
#   derivative_col = "dZ",
#   diffusion_formula = "sigma * (1 + tau * abs(Z))",
#   max_iter = 20,
#   tolerance = 1e-4
# )
# 
# # Check convergence
# print(sde$converged)
# print(sde$iterations)

## ----cv, eval=FALSE-----------------------------------------------------------
# # Block cross-validation (recommended for time series)
# cv <- cross_validate(
#   equation,
#   data = data,
#   derivative_col = "dZ",
#   k = 5,
#   method = "block"
# )
# 
# print(cv)
# plot(cv)

## ----trajectory, eval=FALSE---------------------------------------------------
# # Simulate from the SDE
# sim <- simulate_trajectory(
#   sde,
#   initial_conditions = c(Z = data$Z[1]),
#   times = data$time,
#   n_sims = 100,
#   method = "euler"
# )
# 
# # Plot with observed data
# plot(sim, observed_data = data)
# 
# # Compare quantitatively
# metrics <- compare_trajectories(sim, data, "time", "Z")
# print(metrics)

## ----qualitative, eval=FALSE--------------------------------------------------
# # Analyze fixed points
# fps <- analyze_fixed_points(equation, "Z", range = c(-10, 150))
# print(fps)
# 
# # Bifurcation analysis (vary a parameter)
# bifurc <- analyze_bifurcations(
#   equation,
#   variable = "Z",
#   parameter = "K",
#   param_range = c(50, 200)
# )
# plot(bifurc)
# 
# # Check against expected behavior
# qual_check <- check_qualitative_behavior(
#   equation, data, "Z",
#   expected_features = list(
#     n_fixed_points = 2,
#     stability_pattern = c("unstable", "stable"),
#     bounded = TRUE
#   )
# )
# print(qual_check)

## ----full-validation, eval=FALSE----------------------------------------------
# validation <- validate_model(
#   equation = equation,
#   sde = sde,
#   data = data,
#   derivative_col = "dZ",
#   variable = "Z",
#   time_col = "time",
#   cv_folds = 5,
#   n_sims = 100,
#   expected_features = list(n_fixed_points = 2)
# )
# 
# print(validation)
# plot(validation)

## ----latex, eval=FALSE--------------------------------------------------------
# # Convert to LaTeX
# latex <- to_latex(equation, variable = "Z", precision = 4)
# cat(latex)
# 
# # Full SDE in LaTeX
# sde_latex <- to_latex(sde, variable = "Z")
# cat(sde_latex)

## ----tables, eval=FALSE-------------------------------------------------------
# # Coefficient table
# coef_table <- coefficient_table(equation, format = "latex")
# cat(coef_table)
# 
# # Model comparison
# models <- list(
#   "Linear" = fit_specified_equation("a + b*Z", data, "dZ"),
#   "Quadratic" = fit_specified_equation("a + b*Z + c*Z^2", data, "dZ"),
#   "Logistic" = fit_specified_equation("r*Z*(1-Z/K)", data, "dZ")
# )
# comparison <- model_comparison_table(models, data, "dZ", format = "markdown")
# cat(comparison)

## ----report, eval=FALSE-------------------------------------------------------
# # Collect all results
# results <- list(
#   data = data,
#   exploration = exploration,
#   equation = equation,
#   sde = sde,
#   validation = validation
# )
# 
# # Generate markdown report (using tempdir() per CRAN policy)
# generate_report(
#   results,
#   output_file = file.path(tempdir(), "analysis_report"),
#   format = "markdown",
#   title = "My Dynamics Analysis"
# )
# 
# # Export to multiple formats (using tempdir() per CRAN policy)
# export_results(results, output_dir = tempdir(), formats = c("rds", "csv", "json"))
# 
# # Save plots (using tempdir() per CRAN policy)
# save_plots(results, output_dir = tempdir(), format = c("png", "pdf"))

## ----summary, eval=FALSE------------------------------------------------------
# # Print summary to console
# print_summary(results)
# 
# # Get template for new analyses (using tempdir() per CRAN policy)
# get_analysis_template(file.path(tempdir(), "my_analysis.R"))

## ----session-info-------------------------------------------------------------
sessionInfo()

