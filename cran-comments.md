## Submission of EmpiricalDynamics 0.1.6 (bug-fix release)

This is a patch release. It repairs the cross-validation of equations produced by
the package's own symbolic search, which failed for every candidate, and the
class of fitted object that the fallback optimiser returns, which had no methods.
No user-facing function was removed or renamed and no argument changed meaning;
the new arguments are additive and default to the previous behaviour. Details in
NEWS.md:

* A search returns its equations with the constants already fitted, as numeric
  literals, so `cross_validate()` was asking `fit_specified_equation()` to
  re-estimate an expression with nothing left to estimate. The new exported
  function `parameterize_equation()` turns those literals into free parameters
  starting at the discovered values, and `cross_validate()` now applies it before
  the folds begin. `bootstrap_parameters()` and `sensitivity_analysis()` failed
  on the same equations for the same reason and are fixed the same way.

* `fit_with_optim()` built objects of class `optim_fit` for which the package
  defined no methods at all, so prediction from them failed with "no applicable
  method for 'predict'". `predict()`, `coef()`, `fitted()`, `residuals()` and
  `print()` methods are now registered for the class.

* `nlsLM()` and `nls()` decide whether they were given weights with `missing()`,
  and `fit_specified_equation()` always passed `weights = weights`; with the
  default `NULL` this produced a zero-length weight vector and every
  Levenberg-Marquardt fit failed into the fallback. The argument is now supplied
  only when it has a value.

* `predict.symbolic_equation()` refused to evaluate an equation that carries its
  constants as literals, on the grounds that it had no coefficients; it now
  evaluates it.

* `cross_validate()` and `bootstrap_parameters()` gained an optional `weights`
  argument, and `cross_validate()` gained `refit_engine`, which can hand the
  per-fold re-estimation to SymbolicRegression.jl's own constant optimisation.
  The default is the R engine, so nothing changes for a user who does not ask,
  and no example, test or vignette requires Julia.

## Test environments

* Local: Fedora Linux 44, R 4.6.1, with Julia 1.12.1 + SymbolicRegression.jl
  1.13.2 present (the Julia paths are optional and are skipped when absent).

## R CMD check results

`R CMD check --as-cran --run-donttest` on the built tarball:

    Status: OK

0 errors | 0 warnings | 0 notes

## Reverse dependencies

None; the package has no reverse dependencies on CRAN.
