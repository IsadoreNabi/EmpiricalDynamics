Based on the detailed documentation provided for **EmpiricalDynamics**, the addition of strategic diagrams would significantly enhance user understanding, particularly for the complex workflows involved in SDE discovery and the hybrid architecture.

Here is the revised documentation with suggested image tags inserted at high-value locations.

# EmpiricalDynamics

**High-Performance & Robust Empirical Discovery of Differential Equations from Time Series Data**

EmpiricalDynamics is a comprehensive toolkit for discovering differential and difference equations from empirical time series data. It combines the statistical power of R with a high-performance Julia backend (via SymbolicRegression.jl) to offer a robust engine capable of recovering physical laws, economic models, and stochastic differential equations from noisy data.

## Performance Benchmarks:

* **Deterministic ODEs:** R² > 0.93 on chaotic systems (Lorenz attractor)
* **Stochastic SDEs:** Drift R² > 0.86, Diffusion R² > 0.59 with GLS refinement
* **Physics Constants:** Precision of 10⁻⁸ recovering π and e from noisy data

## 🚀 Key Features

| Feature | Description |
| --- | --- |
| **Hybrid R/Julia Architecture** | Seamlessly leverages Julia's speed for symbolic regression while keeping R's statistical ecosystem for analysis |
| **Multi-threaded Search** | Automatic parallelization across CPU cores via Julia threads |
| **SDE Discovery** | Full pipeline for Stochastic Differential Equations: drift + diffusion recovery |
| **GLS Iterative Refinement** | Iterative Generalized Least Squares for heteroscedastic data |
| **TVR Differentiation** | Total Variation Regularized derivatives for robust estimation in noisy data |
| **Physics-Informed** | Automatic detection of fundamental constants (π, e, φ) |
| **Production Ready** | Handles edge cases (NaNs, outliers, single-variable inputs) gracefully |

---

## 📦 Installation

### Prerequisites

* **R** (>= 4.0.0)
* **Julia** (>= 1.9): Must be installed and available in your system path

### Install from Source

```r
# 1. Install R dependencies
install.packages(c("CVXR", "minpack.lm", "signal", "lmtest", "tseries", 
                   "ggplot2", "gridExtra", "JuliaCall", "devtools"))

# 2. Install EmpiricalDynamics
devtools::install()

# 3. (Optional) Enable Julia parallelization
# Add to ~/.Renviron:
# JULIA_NUM_THREADS=8

```

### First-Time Setup

```r
library(EmpiricalDynamics)

# Initialize Julia backend (installs required Julia packages)
# This may take a few minutes on first run
julia <- setup_julia_backend()

```

---

## ⚡ Quick Start

### Example 1: Discovering a Physical Law

Recovering  from noisy observations:

```r
library(EmpiricalDynamics)

# Generate noisy data
set.seed(42)
t <- seq(-3, 3, length.out = 200)
y <- pi * sin(t) + exp(1) * t + rnorm(200, 0, 0.05)
data <- data.frame(y = y, t = t)

# Symbolic search
results <- symbolic_search(
  target = data$y,
  predictors = data[, "t", drop = FALSE],
  complexity_limit = 15,
  n_iterations = 500,
  parsimony = 0.001
)

# Best equation
best <- select_equation(results, criterion = "knee")
print(best)
# Expression: (sin(t) * 3.14159) + (t * 2.71828)
# R²: 0.9999

```

### Example 2: Chaotic System (Lorenz Attractor)

Recovering the Lorenz equations from trajectory data:

```r
# Assuming you have trajectory data with X, Y, Z and their TVR derivatives
data_lorenz <- data.frame(X = X, Y = Y, Z = Z, dX = dX_tvr, dY = dY_tvr, dZ = dZ_tvr)

# Add auxiliary transformations
data_lorenz$XY <- X * Y
data_lorenz$XZ <- X * Z
data_lorenz$YmX <- Y - X

# Search for dx/dt = σ(y - x)
results_dX <- symbolic_search(
  target = data_lorenz$dX,
  predictors = data_lorenz[, c("X", "Y", "Z", "YmX")],
  complexity_limit = 20,
  n_iterations = 1000,
  parsimony = 0.0001
)
# Achieves R² > 0.93

```

---

## 🔬 SDE Discovery: Complete Workflow

Recovering Stochastic Differential Equations of the form:

Where  is the drift and  is the diffusion coefficient.

### The Challenge

In SDEs, the observed derivative contains both signal (drift) and noise (diffusion):


Direct symbolic regression on  recovers drift poorly when Signal-to-Noise Ratio (SNR) < 1.

### Solution: Hybrid GLS + Symbolic Search

```r
library(EmpiricalDynamics)

# ═══════════════════════════════════════════════════════════════════
# STEP 1: Compute derivatives with TVR
# ═══════════════════════════════════════════════════════════════════

data_spec <- specify_variables(data, endogenous = "Z", exogenous = "X", time_col = "time")

tvr_result <- compute_derivative(
  Z = data_spec$Z,
  t = data_spec$time,
  method = "tvr",
  lambda = "auto",
  solver = "osqp"
)

dZ_tvr <- tvr_result$derivative

# ═══════════════════════════════════════════════════════════════════
# STEP 2: Initial drift search (symbolic)
# ═══════════════════════════════════════════════════════════════════

# Prepare predictors with transformations
data_full <- create_transformations(data, variables = c("X", "Z"))
data_full$dZ <- dZ_tvr
data_full$sinX <- sin(data_full$X)
data_full$absX <- abs(data_full$X)
data_full$Z3 <- data_full$Z^3

drift_initial <- symbolic_search(
  target = data_full$dZ,
  predictors = data_full[, c("X", "Z", "sinX", "Z3")],
  complexity_limit = 20,
  n_iterations = 600,
  parsimony = 0.0002
)

# ═══════════════════════════════════════════════════════════════════
# STEP 3: GLS iterative refinement
# ═══════════════════════════════════════════════════════════════════

sde_refined <- estimate_sde_iterative(
  target = data_full$dZ,
  predictors = data_full[, c("X", "Z", "sinX", "Z3", "absX")],
  data = data_full,
  initial_drift = select_equation(drift_initial, "best"),
  max_iter = 8,
  tol = 0.001
)

# GLS iteratively:
# 1. Estimates diffusion from residuals
# 2. Calculates weights = 1/σ²
# 3. Re-estimates drift with weighted least squares
# 4. Repeats until convergence

# ═══════════════════════════════════════════════════════════════════
# STEP 4: Final diffusion search (symbolic)
# ═══════════════════════════════════════════════════════════════════

# Use refined residuals
residuals_refined <- sde_refined$final_residuals
resid_magnitude <- abs(residuals_refined) / sqrt(dt)

# Smooth with median filter (more robust than mean)
resid_smoothed <- stats::runmed(resid_magnitude, k = 21)

diffusion_final <- symbolic_search(
  target = resid_smoothed,
  predictors = data_full[, c("X", "absX")],
  complexity_limit = 12,
  n_iterations = 400,
  parsimony = 0.0005
)

# ═══════════════════════════════════════════════════════════════════
# STEP 5: Construct final SDE model
# ═══════════════════════════════════════════════════════════════════

sde_final <- construct_sde(
  drift = sde_refined$drift,
  diffusion = select_equation(diffusion_final, "best"),
  variable = "Z"
)

print(sde_final)

```

### Benchmark Results

**Test system:** 

| Method | Drift R² | Diffusion R² |
| --- | --- | --- |
| **Direct symbolic search** | 0.82 | 0.07 |
| **GLS + Symbolic (recommended)** | **0.86** | **0.59** |

### When to Use Each Method

| Scenario | Recommended Method |
| --- | --- |
| Deterministic ODE (no noise) | Direct `symbolic_search()` |
| SDE with SNR > 2 | Direct `symbolic_search()` on dZ, then residuals |
| SDE with SNR < 1 | GLS iterative + symbolic (as above) |
| Unknown noise structure | `model_conditional_variance(method = "symbolic")` |

---

## 🧮 Differentiation Methods

Computing derivatives from noisy time series is critical. We offer multiple methods:

```r
# Total Variation Regularization (recommended for noisy data)
deriv_tvr <- compute_derivative(Z, t, method = "tvr", lambda = "auto", solver = "osqp")

# Savitzky-Golay filter (fast fallback)
deriv_sg <- compute_derivative(Z, t, method = "savgol", p = 3, n = 15)

# Cubic spline (smooth data)
deriv_spline <- compute_derivative(Z, t, method = "spline")

# Spectral (periodic data)
deriv_spectral <- compute_derivative(Z, t, method = "spectral")

```

### Comparison

| Method | Noise Robustness | Discontinuities | Speed | Best For |
| --- | --- | --- | --- | --- |
| **TVR** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | Noisy economic/financial data |
| **Savgol** | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | General purpose, fallback |
| **Spline** | ⭐⭐ | ⭐ | ⭐⭐⭐⭐ | Smooth physical data |
| **Spectral** | ⭐ | ⭐ | ⭐⭐⭐⭐⭐ | Periodic signals |

---

## 🛠️ Advanced Configuration

### Julia Backend Parameters

```r
# Direct Julia call with full control
results <- julia_call(
  "run_symbolic_search_robust",
  data_path,
  response_col,
  predictor_cols,
  max_complexity = 20L,      # Tree size limit
  n_iterations = 500L,       # Genetic algorithm generations
  parsimony = 0.001,         # Complexity penalty
  custom_operators = ""      # e.g., "safe_log,safe_sqrt"
)

```

### Parsimony Guidelines

| Value | Use Case |
| --- | --- |
| **0.1** | Quick exploration, very simple models |
| **0.01** | Standard regression |
| **0.001** | Physics discovery (recommended) |
| **0.0001** | Complex systems (Lorenz, etc.) |

### Parallelization

```r
# Check Julia threads
julia_eval("Threads.nthreads()")

# Enable more threads (before julia_setup)
Sys.setenv(JULIA_NUM_THREADS = "8")

```

---

## 📊 Diagnostic Tools

### Residual Analysis

```r
# After fitting drift
residuals <- data$dZ - predict(drift_model, data)

# Comprehensive diagnostics
diag <- residual_diagnostics(
  residuals = residuals,
  data = data,
  predictors = c("X", "Z"),
  max_lag = 20,
  plot = TRUE
)

print(diag)
# Ljung-Box (autocorrelation): p=0.XXX
# ARCH-LM (conditional het.): p=0.XXX  
# Breusch-Pagan (het.): p=0.XXX
# Jarque-Bera (normality): p=0.XXX

```

### Model Validation

```r
# Cross-validation
cv_result <- cross_validate(
  equation = best_drift,
  data = data,
  target_col = "dZ",
  k = 5
)

# Trajectory simulation
sim <- simulate_trajectory(
  sde_model = sde_final,
  initial_state = c(Z = 0),
  t_span = c(0, 30),
  dt = 0.01,
  n_paths = 100
)

plot(sim)

```

---

## 📈 Complete Example: Economic Growth Model

Discovering a Solow-Swan type growth equation from GDP data:

```r
library(EmpiricalDynamics)

# Load data
data <- read_empirical_data("gdp_quarterly.csv")
data <- specify_variables(data, endogenous = "GDP", exogenous = c("Capital", "Labor"), 
                          time_col = "Quarter")

# Compute growth rate
deriv <- compute_derivative(data$GDP, data$Quarter, method = "tvr", lambda = "auto")
data$dGDP <- deriv$derivative

# Create transformations
data <- create_transformations(data, variables = c("GDP", "Capital", "Labor"))

# Search for growth equation
results <- symbolic_search(
  target = data$dGDP,
  predictors = data[, c("GDP", "Capital", "Labor", "log_Capital", "log_Labor")],
  complexity_limit = 20,
  n_iterations = 800,
  parsimony = 0.001
)

# Select and validate
best <- select_equation(results, criterion = "knee")
cv <- cross_validate(best, data, "dGDP", k = 5)

# Export results
export_results(best, cv, format = "latex", file = "growth_model.tex")

```

---

## 📚 Function Reference

### Core Functions

| Function | Description |
| --- | --- |
| `symbolic_search()` | Main symbolic regression interface |
| `select_equation()` | Select from Pareto front (best, simplest, knee) |
| `compute_derivative()` | Numerical differentiation (TVR, Savgol, etc.) |
| `specify_variables()` | Annotate data with variable roles |
| `create_transformations()` | Auto-generate log, sqrt, square variants |

### SDE Functions

| Function | Description |
| --- | --- |
| `estimate_sde_iterative()` | GLS iterative refinement for SDEs |
| `model_conditional_variance()` | Estimate diffusion from residuals |
| `construct_sde()` | Combine drift + diffusion into SDE object |
| `simulate_trajectory()` | Monte Carlo SDE simulation |

### Diagnostics

| Function | Description |
| --- | --- |
| `residual_diagnostics()` | Ljung-Box, ARCH-LM, Breusch-Pagan, Jarque-Bera |
| `cross_validate()` | K-fold cross-validation |
| `validate_model()` | Out-of-sample prediction metrics |

---

## 🔧 Troubleshooting

### Julia Threads Not Working

```r
# Check current threads
julia_eval("Threads.nthreads()")  # Shows 1?

# Solution: Set BEFORE starting R
# Windows: System Environment Variables → JULIA_NUM_THREADS = 8
# Linux/Mac: Add to ~/.Renviron: JULIA_NUM_THREADS=8
# Then restart RStudio completely

```

### TVR Solver Fails

```r
# If OSQP fails, try SCS (slower but more robust)
deriv <- compute_derivative(Z, t, method = "tvr", solver = "scs")

# If all solvers fail, use Savitzky-Golay
deriv <- compute_derivative(Z, t, method = "savgol", p = 3, n = 15)

```

### Low R² on SDE Diffusion

This is expected when SNR < 1. The diffusion term is inherently harder to recover because:

1. It's estimated from residuals (which contain drift estimation error)
2. The stochastic noise masks the systematic pattern

**Solutions:**

* Use GLS iterative refinement (improves diffusion R² by ~8x)
* Increase data quantity (more observations)
* Reduce measurement noise if possible
* Accept that diffusion R² > 0.5 is often the practical ceiling

---

## 📖 Citation

If you use EmpiricalDynamics in your research, please cite:

```bibtex
@software{EmpiricalDynamics,
  title = {EmpiricalDynamics: Robust Discovery of Differential Equations with R and Julia},
  author = {Gómez Julián, José Mauricio},
  year = {2025},
  url = {https://github.com/IsadoreNabi/EmpiricalDynamics},
  version = {0.1.0}
}

```

## 📜 License

MIT License. See LICENSE for details.

---

## 🗺️ Roadmap

* [ ] GPU acceleration for symbolic search
* [ ] Automatic hyperparameter tuning
* [ ] Neural-guided symbolic regression
* [ ] Time-varying parameter detection
* [ ] Integration with Stan for Bayesian SDEs

*Built with ❤️ for Science and Economics.*