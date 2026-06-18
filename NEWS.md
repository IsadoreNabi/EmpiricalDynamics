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
