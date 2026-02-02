# =============================================================================
# 06_DATA_PREP
# =============================================================================
# The purpose of this module is to serve a specific role within the simulation,
# optimization, and evaluation pipeline for robust location estimators under
# contamination. Comments focus on experimental design, reproducibility, and
# statistical meaning rather than user instructions.
# =============================================================================

# ===== MODULE RUN TRACKING =====
# Lightweight run-tracking utilities used across modules to confirm that a script was
# successfully sourced. The environment stores per-module status and timestamps, which helps debug
# partial loads and ordering issues in multi-file pipelines without affecting the experiment logic.
.MOD_STATUS <- if (exists(".MOD_STATUS", inherits = TRUE) && is.environment(.MOD_STATUS)) {
  .MOD_STATUS
} else {
  new.env(parent = emptyenv())
}

# Records that a module finished loading and prints a standardized confirmation line.
# This supports reproducibility and debugging by making it obvious which modules executed, when they
# executed, and any extra context (e.g., run mode) supplied by the caller.
mark_module_done <- function(module_id, extra = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .MOD_STATUS[[module_id]] <- list(done = TRUE, time = ts, extra = extra)
  
  cat(sprintf("[MODULE DONE] %s | %s%s\n",
              ts, module_id,
              if (!is.null(extra)) paste0(" | ", extra) else ""))
  flush.console()
  
  invisible(TRUE)
}

# Checks whether a given module_id has been marked as completed in the run-tracker.
# This is used for sanity checks in interactive sessions or orchestrated runs to verify that required
# dependencies were loaded before downstream modules are executed.
is_module_done <- function(module_id) {
  x <- try(.MOD_STATUS[[module_id]], silent = TRUE)
  is.list(x) && isTRUE(x$done)
}

# Injects controlled contamination into a sample by replacing k observations with
# synthetic outliers whose location/shape is defined by a contamination type. The magnitude is scaled
# by MAD (or a robust fallback) to be comparable across distributions, enabling fair robustness testing
# under multiple realistic contamination mechanisms.
inject_outliers_realistic <- function(sample,
                                      contamination_rate = 0.05,
                                      outlier_scale_mad = 12,
                                      type = c("upper_tail", "symmetric_t", "lower_tail", "point_mass",
                                               "clustered_upper", "clustered_symmetric",
                                               "mixture_bimodal_near", "mixture_bimodal_far"),
                                      df_t = 3,
                                      # Parameters controlling clustered 
                                      # contamination: instead of IID outlier locations,
                                      # outliers can appear in contiguous blocks
                                      # to mimic batch/systematic corruption (e.g., sensor bursts
                                      # or logging glitches). cluster_min/max 
                                      # set the typical block length and affect dependence structure.
                                      # clustered contamination knobs
                                      cluster_min = 3L,
                                      cluster_max = 12L,
                                      # Parameters controlling bimodal mixture 
                                      # contamination: generates outliers from two
                                      # modes around the median at ±delta, 
                                      # where delta scales with outlier_scale_mad*MAD. “near” vs “far”
                                      # only changes the distance multiplier, 
                                      # allowing controlled tests of moderate vs extreme bimodality.
                                      # mixture knobs (distance multiplier 
                                      #                relative to outlier_scale_mad * MAD)
                                      mixture_near_factor = 0.5,
                                      mixture_far_factor  = 1.5) {
  type <- match.arg(type)
  n <- length(sample)
  if (!n) return(sample)
  k <- ceiling(n * contamination_rate)
  if (k <= 0) return(sample)
  
  # Computes robust location/scale anchors (median and MAD) to parameterize outlier
  # magnitude in a distribution-agnostic way. If MAD is zero/invalid (e.g., constant samples), it falls
  # back to SD, and finally to 1, preventing degenerate scaling and ensuring contamination is applied.
  med  <- stats::median(sample)
  madv <- stats::mad(sample, constant = 1, na.rm = TRUE)
  if (!is.finite(madv) || madv <= 0) {
    # fallback: robust-ish SD (if sd=0, use 1 to avoid fixed-point scaling)
    sdv <- stats::sd(sample)
    madv <- if (is.finite(sdv) && sdv > 0) sdv else 1
  }
  
  # Helper to select which indices will be contaminated. In IID mode, it draws k unique
  # positions uniformly. In clustered mode, it builds one or more contiguous blocks of indices to mimic
  # temporally/spatially correlated corruption, then trims/fills to reach exactly k unique indices.
  # --- helper: pick k indices either IID or in contiguous clusters -------
  pick_indices <- function(n, k, clustered = FALSE,
                           cluster_min = 3L, cluster_max = 12L) {
    if (!clustered) {
      return(base::sample.int(n, size = k, replace = FALSE))
    }
    cluster_min <- max(1L, as.integer(cluster_min))
    cluster_max <- max(cluster_min, as.integer(cluster_max))
    take <- integer(0)
    guard <- 0L
    while (length(take) < k && guard < 2000L) {
      guard <- guard + 1L
      L <- base::sample(seq.int(cluster_min, cluster_max), size = 1L)
      start <- base::sample.int(n, size = 1L)
      idx <- start + seq.int(0L, L - 1L)
      idx <- idx[idx <= n]
      take <- unique(c(take, idx))
    }
    if (length(take) > k) take <- base::sample(take, k, replace = FALSE)
    if (length(take) < k) {
      rest <- setdiff(seq_len(n), take)
      if (length(rest)) take <- c(take, base::sample(rest, min(k - length(take), length(rest))))
    }
    take
  }
  
  # Determines whether clustered index selection is used and draws the final set of
  # contaminated positions. Clustered variants intentionally introduce dependence in the contamination
  # pattern; non-clustered variants keep contamination IID. This separation enables controlled studies
  # of estimator robustness under different corruption structures.
  clustered <- type %in% c("clustered_upper", "clustered_symmetric")
  idx <- pick_indices(n, k, clustered = clustered,
                      cluster_min = cluster_min, cluster_max = cluster_max)
  
  # Generates replacement values for the chosen indices according to the contamination
  # mechanism. Upper/lower tail use one-sided Gaussian magnitudes; symmetric_t uses heavy tails; point
  # mass forces identical extreme values; bimodal mixtures create two outlier modes around ±delta.
  # --- generate replacement values --------------------------------------
  if (type %in% c("upper_tail", "clustered_upper")) {
    bump <- abs(stats::rnorm(length(idx)))
    outlier_values <- med + outlier_scale_mad * madv * bump
    
  } else if (type %in% c("lower_tail")) {
    bump <- abs(stats::rnorm(length(idx)))
    outlier_values <- med - outlier_scale_mad * madv * bump
    
  } else if (type %in% c("symmetric_t", "clustered_symmetric")) {
    bump <- stats::rt(length(idx), df = df_t)
    outlier_values <- med + outlier_scale_mad * madv * bump
    
  } else if (type == "point_mass") {
    outlier_values <- rep(med + outlier_scale_mad * madv, length(idx))
    
  } else if (type %in% c("mixture_bimodal_near", "mixture_bimodal_far")) {
    # Two-mode mixture centered around +/- delta from the median.
    # "near" vs "far" only changes delta (distance between modes).
    factor <- if (type == "mixture_bimodal_near") mixture_near_factor else mixture_far_factor
    delta  <- outlier_scale_mad * madv * factor
    # random sign to pick component
    sgn <- ifelse(stats::runif(length(idx)) < 0.5, -1, 1)
    outlier_values <- stats::rnorm(length(idx), mean = med + sgn * delta, sd = madv)
    
  } else {
    stop(sprintf("inject_outliers_realistic: unsupported type '%s'", type))
  }
  
  # Applies the contamination by overwriting selected sample positions with the
  # generated outlier values. This preserves sample size while altering the empirical distribution in
  # a controlled manner, allowing consistent comparison of estimator behavior under fixed corruption
  # rates and mechanisms.
  sample[idx] <- outlier_values
  sample
}

# ============================ SCENARIOS & DATA ==========================
# Scenario-building and data-prep functions define the experimental design: 
# the grid of contamination settings (type/rate/scale/etc.) and the 
# simulation outputs (estimator components). These objects are the foundation
# of GA optimization and subsequent evaluation across distributions.

## Dispatcher light/full
# Routes scenario construction to either the comprehensive FULL design or the faster
# LIGHT design. This keeps downstream code consistent by providing a single entry point that returns
# a scenarios data frame with the same schema, while controlling compute via mode selection.
build_scenarios <- function(mode = c("full","light")) {
  mode <- match.arg(mode)
  if (mode == "full") build_scenarios_full() else build_scenarios_light()
}

# Computes the vector of estimator outputs for a single sample x using the global
# estimator registry. Ensures consistent ordering and naming across modules. Non-finite values are
# handled defensively: if an estimator fails, it falls back to the sample median to keep pipelines
# stable during large simulation sweeps.
compute_components_vector <- function(x) {
  # defensive: remove non-finite values if they appear
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, N_EST))
  # Use the registry to preserve consistent mapping: estimator names <-> component functions
  comps <- vapply(ESTIMATOR_REGISTRY, function(f) {
    out <- tryCatch(f(x), error = function(e) NA_real_)
    if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
  }, numeric(1))
  as.numeric(comps)
}

# The helpers below support reproducible simulation via local seeding and CRN
# (Common Random Numbers). CRN ensures different estimators are 
# compared on the same underlying random draws for each (family, n), 
# reducing Monte Carlo noise and improving cross-run comparability.

#Helper: evaluate an expression under a local seed without affecting global RNG -------

# --- PREP with CRN: uses sim-seeds by (family, n) for reproducible sampling
# 'crn_env' comes from make_crn_indices(); key = paste(family, n, sep="__")

# Prepares (and optionally caches) the scenario list for a distribution family by
# generating a pooled population from the parameter grid, then sampling 
# repeated datasets per scenario and sample size, computing estimator 
# components each time. Supports CRN for variance reduction and optional 
# scenario subsampling to enable CV and staged/halving workflows
# without regenerating data.
prep_scenarios <- function(dist_name, dist_param_grid,
                           sample_sizes = c(50, 100, 300, 500, 1000, 2000),
                           num_samples = 80,
                           scenario_mode = c("full","light"),
                           seed = 123,
                           crn_env = NULL,
                           fam_key = NULL,
                           # Optional scenario subset input: when provided, 
                           #prep_scenarios only simulates the
                           # requested rows instead of the full scenario universe. 
                           #This is used for CV folds, scenario_frac
                           # downsampling, and multi-stage pipelines. 
                           #subset_tag labels the subset for caching and auditing.
                           # scenario subsampling (prepare ONLY this subset)
                           scenario_subset = NULL,
                           subset_tag = NULL,
                           use_cache = TRUE) {
  scenario_mode <- match.arg(scenario_mode)
  
  # Builds the scenario table according to the requested mode unless a scenario_subset
  # is provided by the caller (e.g., CV or staged subsampling). This keeps the pipeline flexible: the
  # same prep function can serve FULL experiments, LIGHT quick runs, and targeted scenario subsets.
  # Build scenario table (full/light) or use provided subset
  scenarios <- if (!is.null(scenario_subset)) {
    as.data.frame(scenario_subset)
  } else {
    build_scenarios(scenario_mode)
  }
  
  # Constructs a stable cache key encoding the family name, parameter-grid signature,
  # sample sizes, replication count, and scenario subset tag. This enables deterministic reuse of the
  # expensive prepped object across runs without depending on external hashing packages, while keeping
  # the key sensitive to substantive changes in the experimental configuration.
  # Cache key (stable enough without extra hashing packages)
  grid_sig <- paste0(nrow(dist_param_grid), "x", ncol(dist_param_grid), "|",
                     paste(names(dist_param_grid), collapse = ","), "|",
                     paste(round(suppressWarnings(colMeans(as.data.frame(dist_param_grid), na.rm = TRUE)), 6),
                           collapse = ","))
  ss_sig <- if (!is.null(subset_tag)) subset_tag else paste0("ALL_", scenario_mode, "_", nrow(scenarios))
  cache_key <- paste(dist_name, grid_sig,
                     paste(sample_sizes, collapse = "_"),
                     paste0("ns", num_samples),
                     ss_sig,
                     sep = "::")
  
  # Returns a cached prepped object if available to avoid recomputation. Caching is
  # safe here because the cache_key captures the experimental configuration and subset identity. This
  # dramatically reduces runtime in multi-stage pipelines (CV, halving) while preserving reproducibility.
  if (isTRUE(use_cache) && exists(cache_key, envir = .prepped_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .prepped_cache, inherits = FALSE))
  }
  
  # Runs the full preparation inside a local seeding scope so the surrounding session’s
  # RNG stream is not permanently advanced. This block generates a pooled population from all parameter
  # grid rows, computes an analytic “true mean” benchmark, then simulates each scenario across sample
  # sizes and replications, applying contamination and computing estimator components into matrices.
  out <- .seed_scope(seed, {
    pool_total <- max(2000L, 3L * max(sample_sizes))
    per_row_n  <- max(1L, ceiling(pool_total / nrow(dist_param_grid)))
    
    combined_population <- unlist(lapply(seq_len(nrow(dist_param_grid)), function(i) {
      param <- as.list(dist_param_grid[i, , drop = FALSE])
      tryCatch(generate_population(dist_name, per_row_n, param), error = function(e) NULL)
    }))
    combined_population <- combined_population[is.finite(combined_population)]
    if (length(combined_population) == 0) stop("Failed to generate combined population.")
    
    row_means <- vapply(seq_len(nrow(dist_param_grid)), function(i) {
      param <- as.list(dist_param_grid[i, , drop = FALSE])
      analytic_mean_from_params(dist_name, param)
    }, numeric(1))
    
    true_mean <- base::mean(row_means)
    
    
    # For each scenario row and each sample size n, generates num_samples replications.
    # If CRN is enabled (crn_env + fam_key), each replication uses a deterministic sim seed so different
    # estimators and folds share identical base samples. Contamination is then injected per scenario, and
    # estimator component vectors are computed and stored as a num_samples x N_EST matrix for that n.
    # 2) For each scenario and each n: num_samples replications with CRN if available
    prepped <- lapply(seq_len(nrow(scenarios)), function(i) {
      sc <- scenarios[i, ]
      sc$contamination_type <- as.character(sc$contamination_type)
      
      components_by_size <- lapply(sample_sizes, function(n) {
        crn_key <- if (!is.null(fam_key)) paste(fam_key, n, sep = "__") else NULL
        sim_seeds <- NULL
        if (!is.null(crn_env) && !is.null(crn_key) && exists(crn_key, envir = crn_env)) {
          obj <- get(crn_key, envir = crn_env)
          sim_seeds <- obj$sim_idx
          if (length(sim_seeds) < num_samples) {
            sim_seeds <- rep(sim_seeds, length.out = num_samples) # deterministic recycling
          } else {
            sim_seeds <- sim_seeds[seq_len(num_samples)]
          }
        }
        
        comps_mat <- matrix(NA_real_, nrow = num_samples, ncol = N_EST)
        for (r in seq_len(num_samples)) {
          # Base sample draw (CRN-seeded if sim_seeds is available)
          if (!is.null(sim_seeds)) {
            replace_flag <- length(combined_population) < n
            s <- .seed_scope(sim_seeds[r],
                             base::sample(combined_population, n, replace = replace_flag))
          } else {
            replace_flag <- length(combined_population) < n
            s <- base::sample(combined_population, n, replace = replace_flag)
          }
          
          if (isTRUE(sc$contamination_rate > 0)) {
            s <- inject_outliers_realistic(
              s,
              contamination_rate = sc$contamination_rate,
              outlier_scale_mad  = sc$outlier_scale_mad,
              type               = sc$contamination_type
            )
          }
          comps_mat[r, ] <- compute_components_vector(s)
        }
        colnames(comps_mat) <- ESTIMATOR_NAMES
        comps_mat
      })
      
      names(components_by_size) <- as.character(sample_sizes)
      
      S_i <- list(
        components_by_size = components_by_size,
        true_mean = true_mean,
        scenario  = sc
      )
      attr(S_i, "scenario_mode") <- scenario_mode
      S_i
    })
    
    # Attaches metadata attributes used throughout the framework. sample_sizes records
    # which n values exist in each scenario’s components_by_size; scenario_mode records whether the
    # scenario table came from FULL or LIGHT design. These attributes let downstream modules adapt
    # behavior (e.g., minibatching) without re-deriving configuration.
    attr(prepped, "sample_sizes")  <- sample_sizes
    attr(prepped, "scenario_mode") <- scenario_mode
    prepped
  })
  
  # Stores the newly created prepped object in the module-level cache for future reuse.
  # Because prep_scenarios is often the most expensive step, caching is essential for iterative research,
  # CV loops, and staged pipelines. The returned object is identical whether computed fresh or loaded.
  if (isTRUE(use_cache)) {
    assign(cache_key, out, envir = .prepped_cache)
  }
  out
}

# ======== BALANCED MINIBATCH HELPERS =====================================
# Utilities to draw a stratified (balanced) minibatch of scenarios. The goal is to keep the minibatch
# representative of the full scenario pool by sampling across contamination “types/rates/scales”, which
# reduces training variance in LIGHT mode and prevents the GA from overfitting to a narrow scenario slice.

# Signature builder used to stratify sampling by scenario characteristics.
# It compresses each scenario into a stable string key based on contamination mechanism and severity.
.scenario_sig <- function(S) {
  sc <- S$scenario
  paste(sc$contamination_type,
        sprintf("r=%.3f", sc$contamination_rate),
        sprintf("s=%.2f", sc$outlier_scale_mad),
        sep = "|")
}

# Draws k scenario indices from pool_idx using balanced, round-robin sampling across strata (signatures).
# Uses a local seed scope for reproducibility without advancing the global RNG stream. If k is large
# relative to the pool, it returns the full pool to avoid unnecessary work and preserve determinism.
balanced_sample_idx <- function(prepped, pool_idx, k, seed) {
  .seed_scope(seed, {
    pool_idx <- unique(as.integer(pool_idx))
    k <- as.integer(k)
    if (length(pool_idx) == 0L || k <= 0L) return(integer(0))
    if (k >= length(pool_idx)) return(pool_idx)
    
    sigs <- vapply(prepped[pool_idx], .scenario_sig, FUN.VALUE = character(1))
    by_sig <- split(pool_idx, sigs)
    
    # Randomize bucket traversal order each round to avoid systematic ordering bias.
    bucket_names <- sample(names(by_sig))
    
    take <- integer(0)
    # Round-robin sampling: take 1 index per (shuffled) bucket repeatedly until reaching k.
    while (length(take) < k && length(by_sig) > 0) {
      for (nm in bucket_names) {
        b <- by_sig[[nm]]
        if (length(b) == 0L) next
        j <- base::sample(b, 1L)
        take <- c(take, j)
        # Remove the selected element so each bucket is sampled without replacement within the round.
        by_sig[[nm]] <- setdiff(b, j)
        if (length(take) >= k) break
      }
      by_sig <- Filter(length, by_sig)
      if (length(by_sig) == 0L) break
      bucket_names <- sample(names(by_sig))
    }
    unique(head(take, k))
  })
}

#Trigger for the module run tracker
mark_module_done("06_data_prep.R")




