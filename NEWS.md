# EmpiricalDynamics 0.1.6

## Cross-validating a discovered equation

A symbolic search returns its equations with the constants already fitted, written
in as numeric literals. `cross_validate()` then handed such an equation to
`fit_specified_equation()` to be re-estimated fold by fold -- and there was
nothing in it left to estimate. Every candidate of every front failed, including
the candidate that carried the true structure, and the failure was silent enough
to look like a property of the data. Three pieces were missing, and all three are
additive.

* **`parameterize_equation()` is new and exported.** It turns the fitted constants
  of a discovered equation into free parameters, starting at the values the search
  found: `(x * x) * 3.9305` becomes `(x * x) * c1` with `start = list(c1 = 3.9305)`.
  The expression is taken apart with R's parser and rebuilt from its syntax tree,
  so a digit inside a variable name (`Z3`) is never mistaken for a coefficient. A
  leading minus is folded into the starting value, and an integer exponent is left
  alone as structure rather than turned into a parameter.

* **`optim_fit` objects have methods.** `fit_with_optim()` -- the last fallback of
  `fit_specified_equation()` -- built an object of class `optim_fit` for which the
  package defined nothing at all, so the moment that fallback was taken prediction
  died with *no applicable method for 'predict' applied to an object of class
  "optim_fit"*. There are now `predict()`, `coef()`, `fitted()`, `residuals()` and
  `print()` methods, and the object carries the expression, the fitted values and
  the residuals it needs to serve them.

* **`cross_validate()` parameterises discovered equations before the folds begin.**
  An equation the user specified, which already has named free parameters, keeps
  taking its starting values from its own coefficients; an equation with no
  constants at all is not refitted -- there is nothing to refit -- and is scored by
  evaluating it on the held-out rows. The returned object reports which of these
  happened in its new `refitted` and `parameterization` entries.

`bootstrap_parameters()` and `sensitivity_analysis()` failed on discovered
equations for the same reason and are fixed the same way. For the bootstrap, the
reported estimate is now the constant as re-estimated on the full sample rather
than the literal the search wrote down, so that the interval is drawn around the
point it belongs to.

## Refitting through the engine that did the discovering

* **`cross_validate(refit_engine = "julia")` is new.** It hands the equation to
  `SymbolicRegression.jl`'s own constant optimisation, which is the optimiser and
  the loss the search itself used, rather than to a reimplementation of it. The
  default remains `"r"`, so that the figure a given equation gets does not depend
  on what happens to be installed on the machine, and so that the whole path
  keeps working where Julia is not available.

  The optimiser and the loss are the search's; the *budget* is not. Inside a
  search the constants are re-optimised on every candidate of every generation,
  so the shipped default is deliberately cheap -- eight BFGS iterations and two
  restarts -- and stops short of the optimum. A refit happens once and is asked
  for the constants that actually minimise the loss, so it is given a budget to
  reach them. Measured over 100 equations from two real runs, the two engines
  then agree on the out-of-sample error to a median relative difference of about
  9e-09.

  Where the Julia engine cannot help, it says so instead of pretending. The
  search's restarts perturb the constants multiplicatively, so a constant the
  search left at zero cannot move, and `abs(c)` with `c` at zero is a kink the
  optimiser does not leave at any budget -- measured, at up to 32 restarts and
  2000 iterations. When that happens the fold is listed in the new
  `stalled_folds` entry of the returned object and a warning is raised, because
  the error reported for such a fold is the error of the equation as the search
  left it and not of the equation refitted.

## Other fixes

* **The Levenberg-Marquardt method was never actually used.** `nlsLM()` and
  `nls()` both decide whether they were given weights with `missing()`, and
  `fit_specified_equation()` always passed `weights = weights`; with the default
  `NULL` that yielded a zero-length weight vector and the fit died with
  *evaluation of fn function returns non-sensible value!*, dropping every call to
  the fallback. The argument is now supplied only when there is something to
  supply, so `method = "LM"` runs as documented.

* **`predict()` on a discovered equation returned `NULL` with a warning.**
  `predict.symbolic_equation()` required coefficients before it would evaluate an
  expression, which an equation carrying its constants as literals does not have,
  although it is perfectly evaluable. It now evaluates such equations.

* **`cross_validate()` and `bootstrap_parameters()` accept `weights`**, so a
  search run under weights can be validated under them too.

* **`square()` and `inv()`** -- which the genetic search emits and which are not R
  functions -- are defined in the package namespace, so equations using them can be
  evaluated and refitted wherever they appear.

* The Julia backend no longer warns about ambiguous soft scope when it is loaded.

# EmpiricalDynamics 0.1.5

## Bug fixes (Julia / SymbolicRegression.jl backend)

This release repairs the JuliaCall path of the Julia backend, which had drifted out
of sync with current `SymbolicRegression.jl` (tested against 1.13.2) and failed for
all users on that path. The underlying engine was working; the R-side glue was not.

* **`setup_julia_backend()` reported the backend as unavailable even when it was
  installed.** The installed-package check used an incorrect hard-coded UUID for
  `SymbolicRegression.jl` (`8254be44-1295-4e6a-a16d-e31fe2c4a48b`); the correct UUID
  is `8254be44-1295-4e6a-a16d-46603ac705cb`. The check now succeeds, so
  `setup_julia_backend()` returns `TRUE` when the package is present.

* **`symbolic_search(backend = "julia")` failed with a single predictor.** A
  one-column `predictors` produced a length-1 character vector that JuliaCall
  transfers to Julia as a scalar `String`, while `SymbolicRegression`'s
  `variable_names` argument requires an `AbstractVector{String}`
  (`TypeError: in keyword argument variable_names ...`). `variable_names` is now
  coerced to a `Vector{String}` on the Julia side, so single-predictor searches work.

* **`symbolic_search(backend = "julia")` failed when extracting results.** The
  Hall-of-Fame extraction sent a multi-statement Julia block through a `julia_command`
  call that parses a single expression, raising
  `ParseError("extra token after end of expression")` after a successful search. The
  extraction block is now wrapped in `begin ... end` (a single expression), so results
  are returned correctly as the Pareto frontier.

These fixes restore end-to-end operation of `symbolic_search(backend = "julia")`
(both single- and multi-predictor) and of `setup_julia_backend()` on the JuliaCall
path. No user-facing API changed.
