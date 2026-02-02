# =============================================================================
# SCENARIO SAMPLING & STAGING UTILITIES
# =============================================================================
# This module defines how contamination scenarios are constructed, sampled, cached, and summarized
# across the experimental pipeline. A “scenario” is one synthetic stress-test environment (type/rate/
# scale) used to evaluate robustness. The utilities below enforce reproducible subset selection,
# balanced coverage across regimes, efficient reuse via caching, and diagnostics to quantify stage
# difficulty and diversity.
# =============================================================================



# ===== MODULE RUN TRACKING =====
# Tracks whether this module has already been sourced in the current session. The shared .MOD_STATUS
# environment acts as a lightweight “load ledger” across files, which helps debug partial sourcing or
# incorrect load order in a multi-module workflow without changing any experimental computations.
.MOD_STATUS <- if (exists(".MOD_STATUS", inherits = TRUE) && is.environment(.MOD_STATUS)) {
  .MOD_STATUS
} else {
  new.env(parent = emptyenv())
}

# Records successful sourcing of this module with a timestamp and optional metadata. This produces a
# standardized console line that makes long runs auditable
mark_module_done <- function(module_id, extra = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .MOD_STATUS[[module_id]] <- list(done = TRUE, time = ts, extra = extra)
  
  cat(sprintf("[MODULE DONE] %s | %s%s\n",
              ts, module_id,
              if (!is.null(extra)) paste0(" | ", extra) else ""))
  flush.console()
  
  invisible(TRUE)
}

# Checks if a module_id was previously marked as done in this R session. This is typically used for
# sanity checks in interactive work or orchestrated runs, ensuring prerequisite modules were sourced
# before executing dependent code, reducing silent failures from missing objects/functions.
is_module_done <- function(module_id) {
  x <- try(.MOD_STATUS[[module_id]], silent = TRUE)
  is.list(x) && isTRUE(x$done)
}



# =============================================================================
# SCENARIO GRID BUILDERS
# =============================================================================
# Two scenario universes are defined:
#   1) FULL  — high-resolution grid used for final evaluation and later stages where coverage matters.
#   2) LIGHT — smaller proxy grid used for fast iteration and early stages where compute is constrained.
# Both grids share the same schema so downstream code (prep, GA, diagnostics) can treat them uniformly.
# =============================================================================

# Builds the full cross-product of contamination rates, outlier severity scales, and contamination
# structures. This defines the “experimental universe” for robustness testing: many combinations yield
# a broad spectrum from clean data to extreme corruption, supporting stable conclusions across regimes.
build_scenarios_full <- function(
    rates  = c(0.00, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30),
    scales = c(1.5, 3, 4.5, 6, 9, 12, 15, 20, 30),
    types  = c(
      "upper_tail", "lower_tail", "symmetric_t", "point_mass",
      "clustered_upper", "clustered_symmetric",
      "mixture_bimodal_near", "mixture_bimodal_far"
    )
) {
  expand.grid(
    contamination_rate = rates,
    outlier_scale_mad  = scales,
    contamination_type = types,
    stringsAsFactors = FALSE
  )
}

# Builds a compact “proxy” scenario grid that preserves diversity (multiple rates/scales/types) while
# reducing the total scenario count dramatically. This is intended for quick runs, hyperparameter
# search, and early-stage optimization where relative ranking matters more than exhaustive coverage.
build_scenarios_light <- function() {
  expand.grid(
    contamination_rate = c(0.00, 0.02, 0.20),
    outlier_scale_mad  = c(3, 9, 20),
    contamination_type = c("upper_tail", "symmetric_t", "clustered_symmetric", "mixture_bimodal_near"),
    stringsAsFactors = FALSE
  )
}



# =============================================================================
# SCENARIO SUBSAMPLING
# =============================================================================

# In-memory cache for scenario subsets. This avoids re-sampling the same subset across repeated calls
# (e.g., CV folds or staged pipelines) when subset_tag + seed + scenario universe size match, ensuring
# deterministic reuse and faster iteration without persisting anything to disk.
.scenario_subset_cache <- new.env(parent = emptyenv())

# Selects a reproducible subset of scenarios from a scenario universe, with stratification across
# contamination_type and coarse bins of rate/scale. Supports either fraction-based sampling (scenario_frac)
# or fixed-size sampling (k). Also enforces a minimum size (min_n) to guarantee enough scenarios for the
# caller (e.g., at least k_folds scenarios for CV validation).
pick_scenario_subset <- function(scenarios,
                                 scenario_frac = 1.0,
                                 k = NULL,
                                 min_n = 1L,
                                 seed = 123,
                                 subset_tag = NULL) {
  stopifnot(is.data.frame(scenarios), nrow(scenarios) >= 1L)
  if (is.null(scenarios$scenario_id)) {
    scenarios <- add_scenario_ids(scenarios)
  }
  scenario_frac <- as.numeric(scenario_frac)
  if (!is.null(k)) k <- as.integer(k)
  min_n <- as.integer(min_n)
  
  # If subset_tag is provided, reuse a previously computed subset for the same tag/seed/universe size.
  # This guarantees that staged runs reference the exact same scenario rows, improving comparability
  # and making downstream exports/audits stable across repeated executions.
  if (!is.null(subset_tag)) {
    key <- paste0(subset_tag, "::", seed, "::", nrow(scenarios))
    if (exists(key, envir = .scenario_subset_cache, inherits = FALSE)) {
      return(get(key, envir = .scenario_subset_cache, inherits = FALSE))
    }
  }
  
  # Determine the target subset size. If k is specified, it dominates; otherwise we compute a size from
  # scenario_frac. Both paths are clamped to [min_n, nrow(scenarios)] to prevent degenerate subsets and
  # to ensure the caller always gets a usable subset even when fractions are tiny or inputs are small.
  target <- if (!is.null(k)) {
    max(min_n, min(nrow(scenarios), k))
  } else {
    if (!is.finite(scenario_frac) || scenario_frac <= 0) stop("scenario_frac must be > 0.")
    max(min_n, min(nrow(scenarios), as.integer(ceiling(nrow(scenarios) * scenario_frac))))
  }
  
  # Perform stratified sampling under a local seed scope so global RNG is not advanced. We build coarse
  # bins for contamination_rate and outlier_scale_mad using quantiles, then stratify on the interaction
  # (type, rate_bin, scale_bin). This reduces the chance of sampling only “easy” or only “hard” regimes.
  out <- .seed_scope(seed, {
    rate_bin  <- cut(scenarios$contamination_rate,
                     breaks = unique(stats::quantile(scenarios$contamination_rate,
                                                     probs = seq(0, 1, length.out = 4),
                                                     na.rm = TRUE)),
                     include.lowest = TRUE, labels = FALSE)
    scale_bin <- cut(scenarios$outlier_scale_mad,
                     breaks = unique(stats::quantile(scenarios$outlier_scale_mad,
                                                     probs = seq(0, 1, length.out = 4),
                                                     na.rm = TRUE)),
                     include.lowest = TRUE, labels = FALSE)
    g <- interaction(scenarios$contamination_type, rate_bin, scale_bin,
                     drop = TRUE, lex.order = TRUE)
    splits <- split(seq_len(nrow(scenarios)), g)
    
    # Allocate per-stratum quotas proportional to stratum size. This approximates “balanced coverage”
    # while still respecting the empirical scenario universe composition. We then sample within each
    # stratum without replacement to build the candidate index set.
    q <- ceiling(target * sapply(splits, length) / nrow(scenarios))
    take <- integer(0)
    for (nm in names(splits)) {
      pool <- splits[[nm]]
      if (!length(pool)) next
      kk <- min(q[[nm]], length(pool))
      take <- c(take, sample(pool, kk))
    }
    take <- unique(take)
    
    # Adjust for rounding effects: if we overshoot, downsample; if we undershoot, fill from the remaining
    # scenarios uniformly. This guarantees the final subset size equals target while staying close to the
    # stratified intent, and preserves determinism under the same seed.
    if (length(take) > target) take <- sample(take, target)
    if (length(take) < target) {
      rest <- setdiff(seq_len(nrow(scenarios)), take)
      need <- target - length(take)
      if (length(rest)) take <- c(take, sample(rest, min(need, length(rest))))
    }
    
    # Optional “coverage repair” across contamination_type: if the target is sufficiently large, we try
    # to guarantee at least min_per_type scenarios per contamination type. If a type is missing, we add
    # one scenario of that type and (when possible) replace a scenario from an overrepresented type so
    # subset size stays fixed and diversity is preserved.
    types_all <- unique(scenarios$contamination_type)
    min_per_type <- 1L
    
    if (length(types_all) > 1L && target >= (length(types_all) * min_per_type)) {
      counts <- table(scenarios$contamination_type[take])
      missing <- setdiff(types_all, names(counts))
      if (length(missing) > 0L) {
        for (tp in missing) {
          pool <- which(scenarios$contamination_type == tp)
          cand <- setdiff(pool, take)
          if (length(cand) == 0L) next
          add_one <- sample(cand, 1L)
          
          counts <- table(scenarios$contamination_type[take])
          replaceable_types <- names(counts[counts > min_per_type])
          if (length(replaceable_types) == 0L) break
          drop_type <- replaceable_types[which.max(as.numeric(counts[replaceable_types]))]
          drop_pos <- which(scenarios$contamination_type[take] == drop_type)
          if (length(drop_pos) == 0L) break
          take[sample(drop_pos, 1L)] <- add_one
        }
      }
    }
    
    scenarios[sort(unique(take)), , drop = FALSE]
  })
  
  # Persist subset into the in-memory cache keyed by tag/seed/universe size. This is especially useful
  # in staged workflows: callers can supply the same subset_tag to guarantee the same subset is reused
  # across multiple stages (HPF, halving, final), without requiring disk persistence.
  if (!is.null(subset_tag)) {
    key <- paste0(subset_tag, "::", seed, "::", nrow(scenarios))
    assign(key, out, envir = .scenario_subset_cache)
  }
  out
}



# Caches used across the broader framework. .prepped_cache stores expensive prepped scenario objects;
# .init_pop_cache stores warm-start GA populations per (family, seed). Both caches reduce repeated work
# during CV/staging loops while keeping results deterministic when inputs (keys) are identical.
.prepped_cache <- new.env(parent = emptyenv())
.init_pop_cache <- new.env(parent = emptyenv())

# Retrieves a stored warm-start population (if available) for a given distribution family and seed.
# Warm-starts accelerate convergence by seeding the GA with previously good candidates, which is most
# valuable in staged pipelines where later stages refine earlier discoveries on harder scenario sets.
get_warm_start <- function(fam, seed) {
  key <- paste(fam, seed, sep = "::")
  if (exists(key, envir = .init_pop_cache)) get(key, envir = .init_pop_cache) else NULL
}

# Stores a warm-start population under (family, seed). This acts like an in-session memory of good
# candidate solutions. Downstream calls can combine warm-starts with fresh random initialization to
# maintain exploration while still benefiting from prior information.
set_warm_start <- function(fam, seed, population) {
  key <- paste(fam, seed, sep = "::")
  assign(key, population, envir = .init_pop_cache)
}



# Implements early-stopping logic on validation history. If the best validation score in the last
# 'patience' checkpoints does not improve beyond a relative threshold (min_delta) compared to earlier
# checkpoints, the function returns TRUE, signaling the GA loop to stop and avoid wasted generations.
should_stop <- function(val_hist, patience = 3L, min_delta = 0.005) {
  L <- length(val_hist)
  if (L < patience + 1L) return(FALSE)
  best_prev   <- min(val_hist[1:(L - patience)])
  best_recent <- min(val_hist[(L - patience + 1L):L])
  improvement <- (best_prev - best_recent) / (abs(best_prev) + 1e-12)
  improvement < min_delta
}



# Normalizes a weight vector onto the probability simplex. It replaces non-finite or non-positive
# entries with eps, then rescales to sum to 1. This is critical because GA operators (mutation/crossover)
# can produce invalid weights; simplex normalization ensures every candidate is a valid convex combination
# of estimators, keeping fitness evaluation stable and interpretable.
.normalize_simplex <- function(w, eps = 1e-12) {
  w <- as.numeric(w)
  w[!is.finite(w)] <- eps
  w[w < eps] <- eps
  s <- sum(w)
  if (s <= eps) {
    w[] <- 1 / length(w)
  } else {
    w <- w / s
  }
  w
}



# Computes Shannon entropy for a positive probability vector. In this framework it serves as a diversity
# diagnostic and/or penalty component: higher entropy means weights are more spread (less dominance by a
# single estimator), while lower entropy indicates concentration, which can signal brittle solutions.
entropy1 <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# Produces a compact “difficulty” profile for a set of scenarios used in a stage. It summarizes typical
# contamination severity (rate*scale), upper-tail difficulty via 95th percentiles, and diversity via
# entropy across contamination types and coarse (type, rate_bin, scale_bin) groups. This helps compare
# how “hard” or “broad” each stage’s scenario subset is.
stage_difficulty <- function(sc_used) {
  sc_used <- as.data.frame(sc_used)
  
  sev <- sc_used$contamination_rate * sc_used$outlier_scale_mad
  
  tab_type <- table(sc_used$contamination_type)
  p_type <- as.numeric(tab_type) / sum(tab_type)
  
  rate_bin  <- cut(sc_used$contamination_rate, breaks = c(-Inf, 0, 0.02, 0.10, 0.20, Inf))
  scale_bin <- cut(sc_used$outlier_scale_mad,  breaks = c(-Inf, 3, 9, 20, Inf))
  g <- interaction(sc_used$contamination_type, rate_bin, scale_bin, drop = TRUE)
  tab_g <- table(g)
  p_g <- as.numeric(tab_g) / sum(tab_g)
  
  data.frame(
    n_scenarios = nrow(sc_used),
    rate_mean = mean(sc_used$contamination_rate),
    rate_q95  = as.numeric(quantile(sc_used$contamination_rate, 0.95, names = FALSE)),
    scale_mean = mean(sc_used$outlier_scale_mad),
    scale_q95  = as.numeric(quantile(sc_used$outlier_scale_mad, 0.95, names = FALSE)),
    severity_mean = mean(sev),
    severity_q95  = as.numeric(quantile(sev, 0.95, names = FALSE)),
    entropy_type = entropy1(p_type),
    entropy_groups = entropy1(p_g)
  )
}



# Samples a minibatch of configuration rows from cfg_grid, either by fraction (frac) or fixed count (k).
# If stratify_cols is provided, selection is approximately proportional across the interaction of those
# columns, improving coverage of the hyperparameter space. This is used in staged search to evaluate a
# representative subset of configs without running the full grid each time.
pick_minibatch_configs <- function(cfg_grid, frac = NULL, k = NULL, seed = 13, stratify_cols = NULL) {
  stopifnot(is.data.frame(cfg_grid), nrow(cfg_grid) >= 1)
  if (is.null(frac) && is.null(k)) stop("Specify 'frac' or 'k'.")
  seed <- .ensure_seed(seed, fallback = 101L)
  set.seed(seed)
  
  target <- if (!is.null(k)) max(1L, min(nrow(cfg_grid), as.integer(k))) else {
    max(1L, min(nrow(cfg_grid), as.integer(ceiling(nrow(cfg_grid) * frac))))
  }
  
  # If no stratification is requested, do a simple random sample without replacement. This is fastest
  # and sufficient when cfg_grid is already well-balanced or when only coarse exploration is needed.
  if (is.null(stratify_cols) || length(stratify_cols) == 0) {
    idx <- sample(seq_len(nrow(cfg_grid)), size = target, replace = FALSE)
    return(cfg_grid[idx, , drop = FALSE])
  }
  
  # Stratified sampling path: create groups via interaction of stratify_cols, allocate per-group quotas
  # proportional to group sizes, then sample within each group. Finally, fix rounding over/under-shoot
  # so the final output has exactly target rows while remaining approximately balanced across strata.
  g <- interaction(cfg_grid[, stratify_cols], drop = TRUE, lex.order = TRUE)
  splits <- split(seq_len(nrow(cfg_grid)), g)
  q <- ceiling(target * sapply(splits, length) / nrow(cfg_grid))
  
  take <- integer(0)
  for (nm in names(splits)) {
    pool <- splits[[nm]]
    if (length(pool) == 0) next
    kk <- min(q[[nm]], length(pool))
    take <- c(take, sample(pool, kk))
  }
  if (length(take) > target) take <- sample(take, target)
  if (length(take) < target) {
    rest <- setdiff(seq_len(nrow(cfg_grid)), take)
    need <- target - length(take)
    take <- c(take, sample(rest, min(need, length(rest))))
  }
  cfg_grid[sort(unique(take)), , drop = FALSE]
}



#Trigger for the module run tracker
mark_module_done("02_scenarios_sampling.R")
