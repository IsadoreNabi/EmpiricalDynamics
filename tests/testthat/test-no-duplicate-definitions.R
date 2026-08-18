## A name defined twice in R/ is decided by collation order, in silence.
##
## This package had two of them, found in August 2026 by censusing every
## top-level definition rather than by reading: `coefficient_change()` and
## `block_bootstrap_indices()`, both in `utils.R` and `validation.R`, both won
## by `validation.R` because collation is alphabetical and the later definition
## overwrites the earlier one. Neither produced an error. The first pair was
## harmless -- the two bodies computed the same thing -- but the second pair was
## not: the `utils.R` version returned a LIST of 200 resamples and made
## `block_size` optional, while the one that actually runs returns a single
## vector and requires it. Anyone reading `utils.R` and calling the function as
## written there gets "argument \"block_size\" is missing".
##
## The failure mode is that the dead copy looks alive. Someone fixing a defect
## edits it, the tests stay green, and nothing changes -- because the edited
## function is the one that never runs.
##
## The census is done with the parser and never with grep: a regular expression
## over the sources counts assignments inside strings and comments.

test_that("no name is defined twice at top level in R/", {
  r_dir <- normalizePath(file.path(test_path(), "..", ".."), mustWork = FALSE)
  r_dir <- file.path(r_dir, "R")
  skip_if_not(dir.exists(r_dir), "R/ not reachable from the test directory")

  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  expect_gt(length(files), 0L)

  seen <- list()
  for (f in files) {
    exprs <- parse(f, keep.source = FALSE)
    for (i in seq_along(exprs)) {
      x <- exprs[[i]]
      if (!is.call(x)) next
      op <- as.character(x[[1L]])[1L]
      if (!op %in% c("<-", "=", "<<-")) next
      lhs <- x[[2L]]
      if (!is.name(lhs)) next
      nm <- as.character(lhs)
      seen[[nm]] <- c(seen[[nm]], basename(f))
    }
  }

  duplicated_names <- names(seen)[vapply(seen, length, integer(1)) > 1L]
  report <- vapply(duplicated_names, function(nm) {
    paste0(nm, " (", paste(seen[[nm]], collapse = ", "), ")")
  }, character(1))

  expect_identical(sort(duplicated_names), character(0),
                   info = paste("defined more than once:",
                                paste(report, collapse = "; ")))
})
