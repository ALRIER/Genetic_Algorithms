# =============================================================================
# ESTIMATOR REGISTRY AND ROBUST LOCATION FUNCTIONALS
# =============================================================================
# This module defines the estimator “building blocks” used by the GA to construct a weighted composite
# location estimator. Each component is implemented defensively so the pipeline remains stable under
# heavy contamination: estimators should not crash on pathological samples and should return a sensible
# fallback (typically the median) if an algorithm fails or yields non-finite outputs.
# =============================================================================



# ===== MODULE RUN TRACKING =====
# Tracks whether this module has already been sourced in the current R session. The shared .MOD_STATUS
# environment is used across modules so long orchestration runs can confirm load order and detect missing
# definitions (e.g., estimator registry not initialized) without altering any stochastic computation.
.MOD_STATUS <- if (exists(".MOD_STATUS", inherits = TRUE) && is.environment(.MOD_STATUS)) {
  .MOD_STATUS
} else {
  new.env(parent = emptyenv())
}

# Records a timestamped “module loaded” flag and prints a standardized console line. This is useful for
# reproducibility and debugging in multi-file projects: the log makes it obvious which modules executed,
# in what order, and whether any module silently failed to source.
mark_module_done <- function(module_id, extra = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .MOD_STATUS[[module_id]] <- list(done = TRUE, time = ts, extra = extra)
  
  cat(sprintf("[MODULE DONE] %s | %s%s\n",
              ts, module_id,
              if (!is.null(extra)) paste0(" | ", extra) else ""))
  flush.console()
  
  invisible(TRUE)
}

# Returns TRUE if the specified module_id has been marked as done. This supports sanity checks before
# running expensive GA stages: if the estimator registry is missing, downstream code would fail or (worse)
# run with incomplete components, so explicit checks can protect experiment integrity.
is_module_done <- function(module_id) {
  x <- try(.MOD_STATUS[[module_id]], silent = TRUE)
  is.list(x) && isTRUE(x$done)
}



# =============================================================================
# HELPER UTILITIES
# =============================================================================

# Computes the 25th, 50th, and 75th percentiles with stable defaults (type=8) and NA handling. These
# quartiles are used to construct robust summary statistics such as the trimean, which combines median
# and quartiles to reduce sensitivity to outliers while retaining reasonable efficiency under clean data.
.qs3 <- function(x) {
  stats::quantile(x, probs = c(0.25, 0.50, 0.75), names = FALSE, type = 8, na.rm = TRUE)
}



# =============================================================================
# ROBUST MEAN VARIANTS
# =============================================================================

# Harmonic mean is highly sensitive to zeros/negatives and non-finite values. This “safe” version restricts
# to positive finite observations; if none exist, it falls back to the sample median so the estimator can
# still be evaluated in contaminated regimes without returning Inf/NaN that would break fitness scoring.
harmonic_mean_safe <- function(x) {
  x <- as.numeric(x); if (!length(x)) return(NA_real_)
  xp <- x[is.finite(x) & x > 0]
  if (!length(xp)) return(stats::median(x, na.rm = TRUE))
  length(xp) / sum(1 / xp)
}

# Geometric mean requires positive values because it uses log(). This safe variant filters to positive finite
# entries and computes exp(mean(log(x))). If positivity is violated (common under some contaminations), it
# returns the sample median as a robust fallback rather than erroring or producing NaNs.
geometric_mean_safe <- function(x) {
  x <- as.numeric(x); if (!length(x)) return(NA_real_)
  xp <- x[is.finite(x) & x > 0]
  if (!length(xp)) return(stats::median(x, na.rm = TRUE))
  exp(base::mean(log(xp)))
}



# =============================================================================
# ROBUST MODE ESTIMATORS
# =============================================================================

# Half-sample mode (HSM) is a robust mode estimator designed to resist outliers. We wrap it in tryCatch
# because mode estimation can fail on degenerate inputs (e.g., all NA, length 0, extreme ties). If the
# result is non-finite, we fall back to the median to keep downstream GA evaluation numerically stable.
mode_hsm_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch(modeest::hsm(x), error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

# Parzen window-based mode estimation depends on kernel density-like procedures that can fail or produce
# unexpected outputs for small samples or pathological distributions. This wrapper forces numeric output,
# catches errors, and falls back to the median when needed so the estimator remains usable as a GA component.
mode_parzen_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch(as.numeric(modeest::mlv(x, method = "parzen")),
                  error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}



# =============================================================================
# TRIMEAN AND M-ESTIMATORS
# =============================================================================

# Tukey’s trimean combines Q1, median, and Q3 with weights (1,2,1)/4. It is more robust than the mean
# (downweights extreme tails) but often more efficient than the median under clean data. It relies only
# on quantiles, making it stable under many contamination types used in this project.
trimean_safe <- function(x) {
  x <- as.numeric(x); q <- .qs3(x); (q[1] + 2*q[2] + q[3]) / 4
}

# Huber M-estimator: computed via robust regression (rlm) with Huber psi. This yields a location estimate
# that behaves like the mean near the center but limits the influence of extreme residuals. We use MAD
# scale estimation and cap iterations; failures revert to the sample median to avoid breaking the GA loop.
huber_mean_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch({
    fit <- MASS::rlm(x ~ 1, psi = MASS::psi.huber, scale.est = "MAD",
                     maxit = 50, na.action = na.omit)
    as.numeric(coef(fit)[1])
  }, error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

# Tukey biweight (bisquare) M-estimator: stronger downweighting of large residuals than Huber. This is
# useful under heavy contamination because extreme points can be almost fully ignored. As with huber_mean,
# we defensively catch failures and return the median when the robust regression cannot produce a finite fit.
biweight_mean_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch({
    fit <- MASS::rlm(x ~ 1, psi = MASS::psi.bisquare, scale.est = "MAD",
                     maxit = 50, na.action = na.omit)
    as.numeric(coef(fit)[1])
  }, error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}



# =============================================================================
# ESTIMATOR REGISTRY
# =============================================================================
# Maps human-readable estimator names to functions. This registry is the “single source of truth” for
# component definitions: the GA treats each entry as one feature in a convex combination, and downstream
# code relies on this ordering to interpret weights, export formulas, and ensure reproducibility across runs.
# =============================================================================

ESTIMATOR_REGISTRY <- list(
  mean        = function(x) base::mean(x, na.rm = TRUE),
  median      = function(x) stats::median(x, na.rm = TRUE),
  trimmed20   = function(x) base::mean(x, trim = 0.20, na.rm = TRUE),
  harmonic    = harmonic_mean_safe,
  geometric   = geometric_mean_safe,
  mode_hsm    = mode_hsm_safe,
  mode_parzen = mode_parzen_safe,
  trimean     = trimean_safe,
  huber       = huber_mean_safe,
  biweight    = biweight_mean_safe
)

# Convenience metadata used throughout the project: estimator names define column labels, and N_EST is
# the dimensionality of the GA weight vector. Many parts of the pipeline (simplex normalization, warm-start
# matrices, perturbation tests) assume length(weights) == N_EST, so this central definition is critical.
ESTIMATOR_NAMES <- names(ESTIMATOR_REGISTRY)
N_EST <- length(ESTIMATOR_NAMES)



# =============================================================================
# GA COMPOSITE ESTIMATOR
# =============================================================================
# Defines the composite estimator optimized by the GA: it computes each registered component on the same
# sample, then returns a convex combination under weights normalized to the simplex. This ensures the final
# estimator is interpretable as a weighted average of named components. If any component fails, it is replaced
# by the sample median so the composite remains well-defined and fitness can be computed without NA cascades.
# =============================================================================

custom_estimator <- function(sample, weights) {
  stopifnot(length(weights) == N_EST)
  w <- .normalize_simplex(as.numeric(weights))
  x <- as.numeric(sample)
  comps <- vapply(ESTIMATOR_REGISTRY, function(f) {
    out <- tryCatch(f(x), error = function(e) NA_real_)
    if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
  }, numeric(1))
  sum(w * comps)
}



# Converts weight vector into readable symbolic formula
# Produces a compact, human-readable representation of the GA solution by pairing each normalized weight
# with its estimator name (rounded for readability). This is used in logs and result tables so researchers
# can quickly interpret which components dominate and how the composite compares across distributions/stages.
weights_to_formula <- function(weights) {
  stopifnot(length(weights) == N_EST)
  w <- round(.normalize_simplex(as.numeric(weights)), 3)
  paste0(paste0(w, "*", ESTIMATOR_NAMES), collapse = " + ")
}



#Trigger for the module run tracker
mark_module_done("04_estimators_registry.R")
