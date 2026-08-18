## Submission of EmpiricalDynamics 0.1.8 (bug-fix release)

This is a patch release. It repairs the iterative GLS estimator for stochastic
differential equations, whose convergence test could never be satisfied and
which returned its last iteration rather than its best; the observation weights
that the Julia symbolic-search backends accepted and silently discarded; and two
internal functions that were each defined twice in the package sources. No
exported function was removed or renamed, no argument changed meaning, and the
one new argument is additive with a default. Details in NEWS.md:

* `estimate_sde_iterative()` measured convergence by comparing the named
  coefficients of successive drift equations. An equation produced by a symbolic
  search carries its constants as numeric literals, so it has no coefficients
  and the comparison returned `Inf` -- even for an equation compared against
  itself, where the true change is zero. No positive tolerance exceeds `Inf`, so
  the loop never stopped early and always exhausted `max_iter`. Convergence now
  requires two conditions together: the relative change of the weighted deviance
  (the quantity the GLS loop minimises) and the relative distance between
  successive fitted functions. A criterion that cannot be evaluated, because a
  candidate predicts non-finite values, is reported as such.

* The same function returned the drift of its final iteration, which with the
  above was always the iteration limit. It now selects over the whole
  trajectory, scored under one common set of weights, since the GLS weights
  change every iteration and deviances computed under different weights are not
  comparable. The default rule, `selection = "blocked_cv"`, is blocked
  cross-validation that refits the constants on the training blocks; blocks are
  contiguous because these are serially dependent series, and deterministic, so
  no seed is involved. The returned object now records `converged`,
  `stop_reason`, the full `history` and which candidate was selected.

* `symbolic_search(backend = "julia")` accepted a `weights` argument and did not
  use it: both Julia arms built the search from the design matrix and response
  alone. Because `symbolic_search_weighted()` is what the iterative estimator
  calls, requesting a weighted fit could silently cancel the weighting step. The
  weights are now passed through a `Dataset`.

* `coefficient_change()` and `block_bootstrap_indices()` were each defined twice,
  in `utils.R` and in `validation.R`, so collation order decided which one was
  used. The dead copy of the second returned a list where the live one returns a
  vector. The dead copies were removed, and a test now censuses every top-level
  definition in `R/` with the R parser and fails if any name is defined twice.

## Test environments

* local: Fedora Linux 44, R 4.6.0

## R CMD check results

0 errors | 0 warnings | 0 notes

The package test suite (268 assertions) passes with no failures, warnings or
skips, both standalone and under `R CMD check --as-cran --run-donttest`.
