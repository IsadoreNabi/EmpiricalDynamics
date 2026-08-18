# =============================================================================
# Builds tests/testthat/fixtures/ed_fronts.rds
#
# The regression fixture is cut from two real runs of this package, so that the
# tests are exercised against the equations a search actually writes rather than
# against equations written by hand to be easy:
#
#   INFIERNO_SDE.RData  <- HARD_RECOVERY_SDE_GLS_SYMBOLIC.R  (simulated SDE)
#   LORENZ.RData        <- HARD_RECOVERY_LORENZ.R            (Lorenz system)
#
# Both sources are simulated data, so nothing empirical travels in the package.
# The rows are thinned to keep the fixture small while leaving the time
# structure that block cross-validation depends on intact.
#
# This script is not part of the package (see .Rbuildignore) and is kept for the
# record: it says exactly where the fixture came from.
# =============================================================================

SOURCE_DIR <- file.path(
  "/mnt/kingston/Carpeta de Estudio/[Teoría Marxista]",
  "6. [Mis Investigaciones]",
  "ANÁLISIS DINÁMICO Y PROBABILÍSTICO DE LOS PRECIOS DE PRODUCCIÓN",
  "PREDICTIONS OF PRICES OF PRODUCTION/ECUACIONES DIFERENCIALES")

THIN <- 10L   # keep every tenth row

read_env <- function(file) {
  env <- new.env(parent = emptyenv())
  load(file.path(SOURCE_DIR, file), envir = env)
  env
}

# Every variable an equation of the front mentions, plus the response.
columns_needed <- function(fronts, data, response) {
  mentioned <- unique(unlist(lapply(fronts, function(front) {
    unlist(lapply(front$expression, function(e) all.vars(parse(text = e))))
  })))
  intersect(unique(c(response, mentioned)), names(data))
}

build_case <- function(name, fronts, data, response) {
  keep_cols <- columns_needed(fronts, data, response)
  rows <- seq(1L, nrow(data), by = THIN)
  list(
    name = name,
    response = response,
    data = data[rows, keep_cols, drop = FALSE],
    equations = do.call(rbind, lapply(seq_along(fronts), function(i) {
      front <- fronts[[i]]
      data.frame(
        front = names(fronts)[i],
        expression = front$expression,
        complexity = front$complexity,
        loss = front$loss,
        rmse_in_sample = sqrt(pmax(front$mse, 0)),
        stringsAsFactors = FALSE
      )
    }))
  )
}

sde <- read_env("INFIERNO_SDE.RData")
lorenz <- read_env("LORENZ.RData")

cases <- list(
  build_case(
    name = "infierno_sde",
    fronts = list(drift = sde$drift_results, drift_iter = sde$drift_results_iter),
    data = sde$data_full,
    response = "dZ"),
  build_case(
    name = "lorenz",
    fronts = list(dX = lorenz$results_dX, dY = lorenz$results_dY,
                  dZ = lorenz$results_dZ),
    data = lorenz$data_lorenz,
    response = "dX")
)

# The Lorenz fronts are one per state variable, and each is validated against
# its own response; they are split so that every equation is scored against the
# derivative it was discovered for.
cases <- list(
  cases[[1]],
  build_case("lorenz_dX", list(dX = lorenz$results_dX), lorenz$data_lorenz, "dX"),
  build_case("lorenz_dY", list(dY = lorenz$results_dY), lorenz$data_lorenz, "dY"),
  build_case("lorenz_dZ", list(dZ = lorenz$results_dZ), lorenz$data_lorenz, "dZ")
)

target <- file.path("tests", "testthat", "fixtures", "ed_fronts.rds")
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
saveRDS(cases, target, version = 3, compress = "xz")

cat("wrote", target, "\n")
for (case in cases) {
  cat(sprintf("  %-14s %3d equations  %4d rows  %d columns  response %s\n",
              case$name, nrow(case$equations), nrow(case$data),
              ncol(case$data), case$response))
}
cat("size:", format(file.size(target), big.mark = ","), "bytes\n")
