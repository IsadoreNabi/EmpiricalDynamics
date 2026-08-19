# EmpiricalDynamics 0.1.9

## The lm path of `analyze_fixed_points()` ignored the fitted coefficients

For an `lm` (or `glm`) equation, the function deparsed the bare right-hand
side of the formula and injected the coefficients into the evaluation
environment under their names. But `lm` coefficient names are *term labels*
-- `I(Z)`, `I(Z^2)`, `(Intercept)` -- which are not symbols of the deparsed
text, so the bindings were never consulted and the bare terms were evaluated
with their coefficients effectively set to 1. Fitted to `dZ = 2*Z - Z^2`,
with coefficients recovered exactly as `2` and `-1`, it returned the fixed
points of `Z + Z^2` -- without a word. Written without `I()` the failure
changed shape instead of going away: a coefficient named `Z` collided with
the search variable itself. The example on the function's own manual page
took this path.

`lm` and `glm` equations are now evaluated through their own `predict()`
machinery (on the response scale for `glm`), so the coefficients enter
exactly as fitted -- `I()` terms, factors and link functions included, and
no name collision is possible. The `nls` and `symbolic_equation` paths,
whose coefficient names are genuine symbols of their expressions, were
correct and are unchanged.

## A name that binds nothing is now reported instead of ignored

Every name the equation consults must now be accounted for -- the search
variable, a fitted coefficient, or a declared exogenous value. Previously,
each of the following was silent: an unnamed `exogenous_values` entry (bound
nothing -- `for (nm in names(x))` over a `NULL` runs zero times); declaring
the search variable (overwritten at every grid point); declaring a name the
equation never uses (no effect); declaring a name that collides with a
fitted coefficient (silently *overwrote the fit*); and a free symbol with no
declared value (evaluation failed point by point into `NA`s, or worse,
picked up whatever the name meant in the enclosing environment). The first
two and the last two are errors now; the never-consulted name is a warning.
A new `coefficient_values` argument takes the deliberate case the collision
error forbids: overriding a fitted coefficient by name, validated against
the coefficients that exist.

## A root the grid landed on exactly was reported twice

Fixed points were detected as `diff(sign(f)) != 0` over the grid. A grid
point falling exactly on a root has sign `0`, which that test counts as two
changes -- one into the zero and one out of it -- so `uniroot()` refined the
same root from both adjacent brackets and the point appeared twice in the
output. An exact zero is now taken as a root in its own right, a crossing
requires a strict sign flip, and two detections closer than half a grid step
-- the instrument's own resolution -- are merged.

## `analyze_bifurcations()` swept parameters that were not there

The parameter was injected as a binding under its name and handed down as an
exogenous value; if that name appeared nowhere in the equation, the binding
was never consulted and the sweep returned `n_param` identical copies of the
same fixed points, each labelled with the absent parameter's name. The
parameter is now validated before the sweep: a fitted coefficient is varied
as a coefficient (which also makes coefficients of `lm` equations, like
`I(Z)`, valid bifurcation parameters -- previously impossible), a variable
of the equation is varied as an exogenous value, anything else is an error,
as is `parameter == variable`.

Two more silences in the same loop: a parameter value whose analysis found
no fixed points was dropped from the output entirely, so the requested grid
could not be reconstructed from the result; and any error was caught and
converted into the same empty frame, so a failure was indistinguishable from
an absence. The returned object now carries `$status`, one row per requested
parameter value -- `ok`, `no_fixed_points`, or `error` with the error's
message -- plus `$mode`, and a `print()` method states the accounting.

# EmpiricalDynamics 0.1.8

## The iterative GLS loop could not converge, and did not say so

`estimate_sde_iterative()` measured convergence with `coefficient_change()`,
which compares the *named coefficients* of two equations. An equation that came
out of a symbolic search carries its constants as literals -- `(Z * Z) * 3.93`,
not `(Z * Z) * a` -- so it has no coefficients, `coef()` returns `NULL` and the
comparison returns `Inf`. Handed the *same equation twice*, where the true
change is exactly zero, it still returned `Inf`. No positive tolerance exceeds
`Inf`, so the loop never stopped early: it always exhausted `max_iter` and left
through the warning.

Convergence is now declared only when **two** conditions hold together: the
relative change of the **weighted deviance** -- the quantity the GLS loop is
actually minimising, measured under the weights of the iteration rather than
with an unweighted RMSE -- and the relative distance between successive fitted
functions. Either alone has a failure mode the other covers: measured on a
production run, the objective moved 0.065% between two iterations whose
predictions differed by 9.2%, while the structure changed and moved *away* from
the known truth. A criterion that cannot be evaluated, because an equation
predicts non-finite values, is now reported as such instead of standing in for
a convergence that never arrives.

## The loop returned its last iteration rather than its best

It returned `f_current`, whatever the final round happened to produce -- and
with the criterion above, the final round was always `max_iter`. Measured on a
production run, the last iteration was further from the known truth than the one
before it.

Keeping "the best so far" as the loop runs would have been wrong too, and less
visibly: the GLS weights change every iteration, so each iteration's deviance is
computed on a different scale and the numbers are not comparable with one
another. Measured on clean data, iteration 2 scored 1.47587 against iteration
1's 1.74385 under their own weights, and 1.47587 against 1.47587 -- identical --
once both were scored under the same ones.

So the choice is made once, at the end, over the whole trajectory. Which rule
does the ranking was settled by measurement rather than by argument: three
candidate rules -- weighted deviance, a description length charging for
constants and structure, and blocked cross-validation that **refits** the
constants on the training blocks -- were compared over sixteen real
trajectories against a known truth, under a decision rule written before the
numbers were seen. Cross-validation with refitting won on both the mean and the
median excess error (0.386 and 0.228, against 0.407 and 0.337 for the plain
deviance), and is now the default, as `selection = "blocked_cv"`. The blocks are
contiguous and not random, which is the part of the name that carries
information: these series have serial dependence, and a randomly held-out point
sits between its own neighbours. `selection = "deviance"` remains available and
is what the default falls back to when no block can be scored.

Refitting is what makes the held-out blocks informative. An equation from a
search carries its constants as literals, so a split that does not re-estimate
them leaves the equation unchanged and the held-out error is the in-sample error
cut into pieces. The blocks are contiguous, because these series carry serial
dependence, and deterministic, so the choice needs no seed; they are aggregated
by their mass of weights rather than by a flat mean.

Candidates that predict non-finite values are excluded from the choice, and
listed in `selection_excluded`: whatever else they are, they cannot be the drift
the function returns.

## The Julia backends discarded the observation weights

`symbolic_search(backend = "julia")` accepted `weights` and never used them.
Both Julia arms built their search from `X` and `y` alone, and the argument
survived only in the fallback to `r_genetic`. Because `symbolic_search_weighted()`
is what `estimate_sde_iterative()` calls, asking for a weighted fit could
silently cancel the GLS step it exists to apply.

The weights now travel through a `Dataset`, which is how this package already
hands them to SymbolicRegression.jl elsewhere. Verified with two halves obeying
different laws: with weights of 1000 and 1e-6, the search recovers the law of
the heavy half instead of a compromise between them.

## Two functions were defined twice

`coefficient_change()` and `block_bootstrap_indices()` each existed in both
`utils.R` and `validation.R`. With no `Collate` field, collation is alphabetical
and `validation.R` won both, silently. The pair of `coefficient_change()` was
harmless -- the two bodies computed the same thing -- but the other was not: the
dead `block_bootstrap_indices()` returned a *list* of 200 resamples and made
`block_size` optional, where the live one returns a single vector and requires
it. Reading the dead copy and calling it as written there raised "argument
\"block_size\" is missing".

The dead copies are gone, and a test now censuses every top-level definition in
`R/` with the parser -- not with a regular expression, which counts assignments
inside strings -- and fails naming the pair if any name is defined twice.

# EmpiricalDynamics 0.1.7

## The records the search produces are accepted by the search's own validator

`symbolic_search()` returns its candidates in `all_equations` as plain lists --
a string, a complexity, an error -- while `get_pareto_set()` and
`select_equation()` wrap the same content in a `symbolic_equation`.
`cross_validate()` insisted on the second and refused the first with *Unknown
equation type. Expected symbolic_equation, nls, or lm object*, so the records
the package itself hands back could not be handed to the package's own
cross-validation. Anything carrying an expression is now accepted and classed
at the door, in `cross_validate()`, `bootstrap_parameters()` and
`sensitivity_analysis()` alike; a list that carries no equation is still
refused rather than guessed at.

This is what stood between a search and its own out-of-sample score for any
caller reading `all_equations` rather than the accessors.

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
