
# =============================================================================
# 07_FITNESS_OBJECTIVES — RESEARCH COMMENTED EDITION
# =============================================================================
# The purpose of this module is to serve a specific role within the simulation,
# optimization, and evaluation pipeline for robust location estimators under
# contamination. Comments focus on experimental design, reproducibility, and
# statistical meaning rather than user instructions.
# =============================================================================

# ===== MODULE RUN TRACKING ===== To know if I ran or not this module.
.MOD_STATUS <- if (exists(".MOD_STATUS", inherits = TRUE) && is.environment(.MOD_STATUS)) {
  .MOD_STATUS
} else {
  new.env(parent = emptyenv())
}

mark_module_done <- function(module_id, extra = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .MOD_STATUS[[module_id]] <- list(done = TRUE, time = ts, extra = extra)
  
  cat(sprintf("[MODULE DONE] %s | %s%s\n",
              ts, module_id,
              if (!is.null(extra)) paste0(" | ", extra) else ""))
  flush.console()
  
  invisible(TRUE)
}

is_module_done <- function(module_id) {
  x <- try(.MOD_STATUS[[module_id]], silent = TRUE)
  is.list(x) && isTRUE(x$done)
}


# ===================== SENSITIVITY ====================

jitter_simplex <- function(w, eps = 0.02, seed = 1) {
  w <- .normalize_simplex(w)
  .seed_scope(seed, {
    noise <- stats::runif(length(w), min = 0, max = eps)
    .normalize_simplex(w + noise)
  })
}

run_sensitivity <- function(w_best, dist_name, dist_param_grid,
                            run_dir, scenario_mode,
                            sample_sizes, num_samples,
                            base_seed = 1001,
                            n_alt = 5) {
  
  rows <- list()
  
  for (i in seq_len(n_alt)) {
    alt_seed <- base_seed + i
    
    # (A) subset alternativo
    sc_u <- build_scenarios(scenario_mode)
    sc_u <- add_scenario_ids(sc_u)
    sc_alt <- pick_scenario_subset(sc_u, scenario_frac = 1.0, seed = alt_seed,
                                   subset_tag = paste0(dist_name, "::sens"))
    
    prepped_alt <- prep_scenarios(
      dist_name, dist_param_grid,
      sample_sizes = sample_sizes,
      num_samples = num_samples,
      scenario_mode = scenario_mode,
      seed = alt_seed,
      scenario_subset = sc_alt,
      subset_tag = paste0(dist_name, "::sens::", alt_seed),
      use_cache = FALSE
    )
    
    det <- fitness_universal(w_best, prepped_alt, return_details = TRUE,
                             run_id = dist_name, ga_seed = alt_seed)
    
    rows[[length(rows)+1]] <- data.frame(
      perturbation_type = "alt_subset_seed",
      delta = NA_real_,
      seed = alt_seed,
      fitness = det$fitness,
      base_obj = det$base_obj,
      mean_loss = det$loss_summary$mean,
      q95_loss  = det$loss_summary$q95,
      max_loss  = det$loss_summary$max,
      stringsAsFactors = FALSE
    )
    
    w_j <- jitter_simplex(w_best, eps = 0.02, seed = alt_seed)
    det2 <- fitness_universal(w_j, prepped_alt, return_details = TRUE,
                              run_id = dist_name, ga_seed = alt_seed)
    
    rows[[length(rows)+1]] <- data.frame(
      perturbation_type = "weight_jitter",
      delta = 0.02,
      seed = alt_seed,
      fitness = det2$fitness,
      base_obj = det2$base_obj,
      mean_loss = det2$loss_summary$mean,
      q95_loss  = det2$loss_summary$q95,
      max_loss  = det2$loss_summary$max,
      stringsAsFactors = FALSE
    )
  }
  
  out <- dplyr::bind_rows(rows)
  safe_write_csv(out, file.path(run_dir, sprintf("sensitivity_%s.csv", dist_name)))
  out
}



# =============================== FITNESS =================================

# --- Huber loss (robusto a outliers) -------------------------------------
huber_loss <- function(r, delta = 1.0) {
  r <- as.numeric(r)
  a <- base::abs(r)
  ifelse(a <= delta, 0.5 * r^2, delta * (a - 0.5 * delta))
}
scenario_weight <- function(sc_row) {
  cr  <- sc_row$contamination_rate
  osm <- sc_row$outlier_scale_mad
  typ <- as.character(sc_row$contamination_type)
  
  base_w  <- 1 + 2 * cr + 0.03 * osm
  type_adj <- switch(typ,
                     "upper_tail"            = 1.05,
                     "lower_tail"            = 1.05,
                     "symmetric_t"           = 1.10,
                     "point_mass"            = 1.00,
                     # clustered outliers (slightly higher risk weight)
                     "clustered_upper"       = 1.12,
                     "clustered_symmetric"   = 1.15,
                     # multimodal / mixture outliers (often challenging for location)
                     "mixture_bimodal_near"  = 1.10,
                     "mixture_bimodal_far"   = 1.20,
                     1.00)
  base_w * type_adj
}


# ===================== HIERARCHICAL GROUPIN =======================

# Bin contamination rate into stable buckets 
.rate_bin <- function(r) {
  # PSOCK safety: NULL / length-0 can happen if a scenario column is missing
  if (is.null(r) || length(r) == 0L || !is.finite(r)) return("r=NA")
  if (r <= 0)       return("r=0")
  if (r <= 0.02)    return("r<=0.02")
  if (r <= 0.10)    return("r<=0.10")
  if (r <= 0.20)    return("r<=0.20")
  "r>0.20"
}


# Bin outlier severity (MAD scale) into buckets
.scale_bin <- function(s) {
  # PSOCK safety: NULL / length-0 can happen if a scenario column is missing
  if (is.null(s) || length(s) == 0L || !is.finite(s)) return("s=NA")
  if (s <= 3)        return("s<=3")
  if (s <= 9)        return("s<=9")
  if (s <= 20)       return("s<=20")
  "s>20"
}


# Scenario -> group key (the “group of scenarios” concept)
scenario_group_key <- function(sc_row,
                               by = c("type","rate_bin","scale_bin")) {
  by <- match.arg(by, several.ok = TRUE)
  parts <- character(0)

  # Use safe extraction (missing columns in a scenario row -> NA)
  get_num <- function(x) if (is.null(x) || length(x) == 0L) NA_real_ else as.numeric(x)
  get_chr <- function(x) if (is.null(x) || length(x) == 0L) NA_character_ else as.character(x)

  if ("type" %in% by)      parts <- c(parts, paste0("type=", get_chr(sc_row$contamination_type)))
  if ("rate_bin" %in% by)  parts <- c(parts, .rate_bin(get_num(sc_row$contamination_rate)))
  if ("scale_bin" %in% by) parts <- c(parts, .scale_bin(get_num(sc_row$outlier_scale_mad)))

  paste(parts, collapse = "|")
}

# Aggregate losses inside each group (inner) then across groups (outer)
hierarchical_objective <- function(losses_by_scen, group_keys,
                                   objective = c("mean","q95","max","mixed"),
                                   inner = c("mean","q95","max"),
                                   inner_q = 0.95,
                                   mix_w_q95 = 0.7,
                                   mix_w_max = 0.3) {
  objective <- match.arg(objective)
  inner     <- match.arg(inner)
  ok <- is.finite(losses_by_scen) & !is.na(group_keys) & nzchar(group_keys)
  losses_by_scen <- losses_by_scen[ok]
  group_keys     <- group_keys[ok]
  if (!length(losses_by_scen)) return(Inf)
  
  by_group <- split(losses_by_scen, group_keys)
  
  # Inner summary: one number per group
  g_stat <- vapply(by_group, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(Inf)
    switch(inner,
           mean = mean(v),
           max  = max(v),
           q95  = as.numeric(stats::quantile(v, probs = inner_q, names = FALSE, type = 8)))
  }, numeric(1))
  
  g_stat <- g_stat[is.finite(g_stat)]
  if (!length(g_stat)) return(Inf)
  
  # Outer objective across groups
  q95_val  <- as.numeric(stats::quantile(g_stat, 0.95, names = FALSE, type = 8))
  max_val  <- max(g_stat)
  mean_val <- mean(g_stat)
  
  switch(objective,
         mean  = mean_val,
         max   = max_val,
         q95   = q95_val,
         mixed = mix_w_q95 * q95_val + mix_w_max * max_val)
}
# ===================== DIAGNOSTICS / REPORTING ======================

group_risk_stats <- function(losses_by_scen, group_keys,
                             inner = c("q95","mean","max"),
                             inner_q = 0.95) {
  inner <- match.arg(inner)
  
  ok <- is.finite(losses_by_scen) & !is.na(group_keys) & nzchar(group_keys)
  losses_by_scen <- losses_by_scen[ok]
  group_keys     <- group_keys[ok]
  if (!length(losses_by_scen)) return(setNames(numeric(0), character(0)))
  
  by_group <- split(losses_by_scen, group_keys)
  
  vapply(by_group, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(Inf)
    switch(inner,
           mean = mean(v),
           max  = max(v),
           q95  = as.numeric(stats::quantile(v, probs = inner_q, names = FALSE, type = 8))
    )
  }, numeric(1))
}

summarize_losses <- function(losses_by_scen) {
  losses_by_scen <- losses_by_scen[is.finite(losses_by_scen)]
  if (!length(losses_by_scen)) {
    return(list(mean = Inf, q95 = Inf, max = Inf))
  }
  list(
    mean = mean(losses_by_scen),
    q95  = as.numeric(stats::quantile(losses_by_scen, 0.95, names = FALSE, type = 8)),
    max  = max(losses_by_scen)
  )
}

# --- Bootstrap q95 con control de RNG y atajos en casos degenerados ------
q95_boot <- function(x, B = 300L, rng_seed = NULL, boot_indices = NULL) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  if (n < 5L)  return(as.numeric(stats::quantile(x, 0.95, names = FALSE, type = 8)))
  if (all(x == x[1])) return(x[1])
  
  s_bkp_exists <- exists(".Random.seed", inherits = FALSE)
  if (s_bkp_exists) s_bkp <- .Random.seed
  if (!is.null(rng_seed)) {
    set.seed(rng_seed)
  } else if (!s_bkp_exists) {
    stats::runif(1)  # garantiza .Random.seed
  }
  
  if (!is.null(boot_indices)) {
    B_eff <- base::min(B, length(boot_indices))
    qs <- vapply(seq_len(B_eff), function(b) {
      xb <- x[ boot_indices[[b]] ]
      as.numeric(stats::quantile(xb, 0.95, names = FALSE, type = 8))
    }, numeric(1))
  } else {
    qs <- replicate(B, {
      xb <- base::sample(x, replace = TRUE)
      as.numeric(stats::quantile(xb, 0.95, names = FALSE, type = 8))
    })
  }
  
  if (s_bkp_exists) .Random.seed <<- s_bkp
  stats::median(qs)
}

# --- Fitness universal: combina objetivo + penalizaciones -----------------
fitness_universal <- function(w, prepped,
                              objective = c("mean","q95","max","mixed"),
                              idx = NULL,
                              use_boot = TRUE, B = 300,
                              # penalizaciones / mixing
                              lambda_instab = 0.0,
                              lambda_entropy = 0.015,
                              max_w_dominance = 0.80,
                              lambda_dominance = 0.05,
                              mix_w_q95 = 0.7,
                              mix_w_max = 0.3,
                              lambda_bias = 0.05,
                              # knobs de robustez
                              huber_delta = 1.0,
                              q95_seed = NULL,
                              use_hierarchical = TRUE,
                              group_by = c("type","rate_bin","scale_bin"),
                              group_inner = c("q95","mean","max"),
                              group_inner_q = 0.95,
                              return_details = FALSE,
                              run_id = NA_character_,
                              ga_seed = NA_integer_,
                              q95_boot_indices = NULL) {
  w <- .normalize_simplex(as.numeric(w))
  objective <- match.arg(objective)
  
  # Subconjunto de escenarios
  Slist <- if (is.null(idx)) prepped else prepped[idx]
  if (length(Slist) == 0L) return(Inf)
  
  if (!is.null(attr(prepped, "scenario_mode")) &&
      identical(attr(prepped, "scenario_mode"), "light") &&
      missing(B)) {
    B <- 60L
  }
  
  losses_by_scen <- vapply(Slist, function(S) {
    ls_per_size <- vapply(S$components_by_size, function(C) {
      est <- as.vector(C %*% w)
      r   <- est - S$true_mean
      base::mean(huber_loss(r, delta = huber_delta), na.rm = TRUE)
    }, numeric(1))
    base::mean(ls_per_size) * scenario_weight(S$scenario)
  }, numeric(1))
  
  group_keys <- vapply(Slist, function(S) {
    scenario_group_key(S$scenario, by = group_by)
  }, character(1))
  
  
  
  # Filtra no finitos; si todo es no-finito, fitness = Inf
  finite_mask <- is.finite(losses_by_scen)
  if (!any(finite_mask)) return(Inf)
  
  if (!all(finite_mask)) {
    losses_by_scen <- losses_by_scen[finite_mask]
    group_keys     <- group_keys[finite_mask]
  }
  
  
  bias_by_scen <- vapply(Slist[finite_mask], function(S) {
    base::mean(vapply(S$components_by_size, function(C) {
      est <- as.vector(C %*% w)
      base::mean(est, na.rm = TRUE) - S$true_mean
    }, numeric(1)))
  }, numeric(1))
  mean_abs_bias <- base::mean(base::abs(bias_by_scen), na.rm = TRUE)
  
  if (isTRUE(use_hierarchical)) {
    # Hierarchical: inner summary per group, then objective across groups
    base_obj <- hierarchical_objective(
      losses_by_scen = losses_by_scen,
      group_keys     = group_keys,
      objective      = objective,
      inner          = match.arg(group_inner),
      inner_q        = group_inner_q,
      mix_w_q95      = mix_w_q95,
      mix_w_max      = mix_w_max
    )
  } else {
    # Flat (your original behavior)
    q95_val  <- if (isTRUE(use_boot)) {
      q95_boot(losses_by_scen, B = B, rng_seed = q95_seed, boot_indices = q95_boot_indices)
    } else {
      as.numeric(stats::quantile(losses_by_scen, 0.95, names = FALSE, type = 8))
    }
    max_val  <- base::max(losses_by_scen)
    mean_val <- base::mean(losses_by_scen)
    
    base_obj <- switch(objective,
                       mean  = mean_val,
                       max   = max_val,
                       q95   = q95_val,
                       mixed = mix_w_q95 * q95_val + mix_w_max * max_val)
  }
  
  # --- Penalizaciones -----------------------------------------------------
  penalty_instab <- lambda_instab * stats::sd(losses_by_scen, na.rm = TRUE)
  
  eps <- 1e-12
  entropy <- -sum(w * base::log(w + eps))
  entropy_norm <- entropy / base::log(length(w))
  pen_entropy <- lambda_entropy * (1 - entropy_norm)
  
  dom_excess <- base::max(0, base::max(w) - max_w_dominance)
  pen_dominance <- lambda_dominance * dom_excess
  
  pen_bias <- lambda_bias * mean_abs_bias
  
  # ----- FINAL FITNESS (same numeric as before) -----
  fitness_total <- base_obj + penalty_instab + pen_entropy + pen_dominance + pen_bias
  
  # Default behavior for GA: return a single scalar
  if (!isTRUE(return_details)) {
    return(fitness_total)
  }
  
  # ----- detailed diagnostics return -----
  loss_summary <- summarize_losses(losses_by_scen)
  
  g_stats <- if (isTRUE(use_hierarchical)) {
    group_risk_stats(losses_by_scen, group_keys,
                     inner = match.arg(group_inner),
                     inner_q = group_inner_q)
  } else {
    setNames(numeric(0), character(0))
  }
  
  return(list(
    fitness = fitness_total,
    base_obj = base_obj,
    penalties = list(
      instab    = penalty_instab,
      entropy   = pen_entropy,
      dominance = pen_dominance,
      bias      = pen_bias
    ),
    loss_summary = loss_summary,
    losses_by_scen = losses_by_scen,
    group_keys = group_keys,
    group_stats = g_stats,
    meta = list(
      run_id = run_id,
      ga_seed = ga_seed,
      use_hierarchical = use_hierarchical,
      objective = objective,
      group_by = group_by,
      group_inner = match.arg(group_inner),
      group_inner_q = group_inner_q
    ),
    w = w
  ))
  
  
}





# ===================== SENSITIVITY: WEIGHT PERTURBATION =====================

weight_perturbation_test <- function(w_star, prepped,
                                     sigma = 0.02, n = 200, seed = 1,
                                     ...) {
  W <- perturb_weights_simplex(w_star, sigma = sigma, n = n, seed = seed)
  
  fits <- apply(W, 1, function(w) {
    fitness_universal(
      w = w,
      prepped = prepped,
      return_details = FALSE,
      ...
    )
  })
  
  list(
    sigma = sigma,
    n = n,
    fitness_mean = mean(fits, na.rm = TRUE),
    fitness_sd   = stats::sd(fits, na.rm = TRUE),
    fitness_q95  = as.numeric(stats::quantile(fits, 0.95, na.rm = TRUE)),
    fitness_min  = min(fits, na.rm = TRUE),
    fitness_max  = max(fits, na.rm = TRUE)
  )
}


# ===================== BASELINE: RANDOM SEARCH =====================

random_search_baseline <- function(prepped, N = 5000,
                                   objective = "q95",
                                   use_hierarchical = TRUE,
                                   group_by = c("type","rate_bin","scale_bin"),
                                   group_inner = "q95",
                                   group_inner_q = 0.95,
                                   seed = 1,
                                   idx = NULL) {
  seed <- .ensure_seed(seed, fallback = 101L)
  set.seed(seed)
  
  n_estimators <- ncol(prepped[[1]]$components_by_size[[1]])  # number of estimators
  
  best <- list(fitness = Inf, w = NULL, details = NULL)
  
  for (i in seq_len(N)) {
    w <- as.numeric(stats::runif(n_estimators))
    w <- w / sum(w)
    
    out <- fitness_universal(
      w = w, prepped = prepped,
      objective = objective,
      idx = idx,
      use_boot = TRUE, B = 0,
      use_hierarchical = use_hierarchical,
      group_by = group_by,
      group_inner = group_inner,
      group_inner_q = group_inner_q,
      return_details = TRUE,
      run_id = paste0("rand_", seed, "_", i),
      ga_seed = seed
    )
    
    if (is.list(out) && is.finite(out$fitness) && out$fitness < best$fitness) {
      best <- list(fitness = out$fitness, w = out$w, details = out)
    }
  }
  
  best
}

# ===================== SENSITIVITY: WEIGHT PERTURBATION  =====================

perturb_weights_simplex <- function(w, sigma = 0.02, n = 200, seed = 1) {
  seed <- .ensure_seed(seed, fallback = 101L)
  set.seed(seed)
  n_dims <- length(w)
  W <- matrix(NA_real_, nrow = n, ncol = n_dims)
  for (i in seq_len(n)) {
    z  <- log(pmax(w, 1e-12)) + stats::rnorm(n_dims, 0, sigma)
    ww <- exp(z)
    ww <- ww / sum(ww)
    W[i, ] <- ww
  }
  W
}

weight_perturbation_test <- function(w_star, prepped,
                                     sigma = 0.02, n = 200, seed = 1,
                                     ...) {
  W <- perturb_weights_simplex(w_star, sigma = sigma, n = n, seed = seed)
  
  fits <- apply(W, 1, function(w) {
    fitness_universal(
      w = w,
      prepped = prepped,
      return_details = FALSE,
      ...
    )
  })
  
  list(
    sigma = sigma,
    n = n,
    fitness_mean = mean(fits, na.rm = TRUE),
    fitness_sd   = stats::sd(fits, na.rm = TRUE),
    fitness_q95  = as.numeric(stats::quantile(fits, 0.95, na.rm = TRUE)),
    fitness_min  = min(fits, na.rm = TRUE),
    fitness_max  = max(fits, na.rm = TRUE),
    fits = fits
  )
}




  # ============================== SAVE HELPERS ==============================
  
  # Utilitario: escribir CSV solo si el objeto existe y tiene filas
  .write_csv_safe <- function(obj, path) {
    if (is.null(obj)) return(invisible(FALSE))
    if (!is.data.frame(obj) || nrow(obj) == 0L) return(invisible(FALSE))
    readr::write_csv(obj, path)
    invisible(TRUE)
  }
  
  save_family_config_outputs <- function(root_out, family_tag, config_tag, res_cv) {
    if (is.null(res_cv) || is.null(res_cv$final)) return(invisible(FALSE))
    fam_dir   <- file.path(root_out, family_tag)
    folds_dir <- file.path(fam_dir, sprintf("FOLDS__%s", config_tag))
    mkdirp(fam_dir); mkdirp(folds_dir)
    
    # --- Final (retrain o mejor fold) ---
    .write_csv_safe(res_cv$final$overall,
                    file.path(fam_dir, sprintf("OVERALL__%s.csv",   config_tag)))
    .write_csv_safe(res_cv$final$scenario_table,
                    file.path(fam_dir, sprintf("SCENARIOS__%s.csv", config_tag)))
    .write_csv_safe(res_cv$final$convergence,
                    file.path(fam_dir, sprintf("CONVERGENCE__%s.csv", config_tag)))
    
    # --- Por fold ---
    fr_list <- res_cv$fold_results %||% list()
    if (length(fr_list) > 0) {
      for (k in seq_along(fr_list)) {
        fr <- fr_list[[k]]
        .write_csv_safe(fr$overall,        file.path(folds_dir, sprintf("fold%02d_OVERALL.csv",    k)))
        .write_csv_safe(fr$scenario_table, file.path(folds_dir, sprintf("fold%02d_SCENARIOS.csv",  k)))
        .write_csv_safe(fr$convergence,    file.path(folds_dir, sprintf("fold%02d_CONVERGENCE.csv", k)))
      }
    }
    invisible(TRUE)
  }
  
  # -------- BIG LOG HELPERS (unified CSV) ---------------------------------
  
  
  .safe_cbind <- function(df, extra) {
    if (is.null(extra) || length(extra) == 0) return(df)
    ex <- as.data.frame(extra, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(df) > 0 && nrow(ex) %in% c(0,1)) {
      ex <- ex[rep(1, nrow(df)), , drop = FALSE]
      rownames(ex) <- NULL
    }
    cbind(df, ex, stringsAsFactors = FALSE)
  }
  
  .parse_config <- function(config_tag) {
    as_num <- function(x) suppressWarnings(as.numeric(x))
    data.frame(
      ps    = as_num(sub(".*\\bps(\\d+).*", "\\1", config_tag)),
      mr    = as_num(sub(".*_mr([0-9.]+).*", "\\1", config_tag)),
      a0    = as_num(sub(".*_a0([0-9.]+).*", "\\1", config_tag)),
      am    = as_num(sub(".*_am([0-9.]+).*", "\\1", config_tag)),
      im    = as_num(sub(".*_im([0-9.]+).*", "\\1", config_tag)),
      t     = as_num(sub(".*_t(\\d+).*", "\\1", config_tag)),
      e     = as_num(sub(".*_e(\\d+).*", "\\1", config_tag)),
      seed  = as_num(sub(".*_seed(\\d+).*", "\\1", config_tag)),
      stage = as_num(sub(".*_STAGE(\\d+)$", "\\1", config_tag)),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  
  # Asegura que la tabla de escenarios tenga las columnas clave, rellenando con NA si faltan
  .ensure_scenario_cols <- function(df) {
    need <- c("distribution","sample_size",
              "contamination_rate","outlier_scale_mad","contamination_type",
              "true_mean","scenario_mode",
              "estimator","mse","mse_q95","mse_max","bias","variance","mad","iqr")
    for (nm in need) if (!nm %in% names(df)) df[[nm]] <- NA
    df
  }
  
  # Asegura que la tabla de convergencia tenga las columnas clave
  .ensure_convergence_cols <- function(df) {
    need <- c("gen","best_train","med_train","best_val","med_val")
    for (nm in need) if (!nm %in% names(df)) df[[nm]] <- NA
    df
  }
  
  # Construye un contexto enriquecido combinando:
  # - lo que viene en 'context' (presupuestos, cobertura, fase, etc.)
  # - lo que se pueda inferir del propio resultado (objective, scenario_mode) y del tag
  .build_context <- function(df, table, run_id, family, config_tag, is_winner, fold, context) {
    cfg <- .parse_config(config_tag)
    # Intento de deducir scenario_mode desde la propia tabla (si existe la col)
    scenemode_in_df <- df$scenario_mode[1] %||% NA
    # coverage_hint: por defecto "full" si scenario_mode == "full"; si no, "partial"
    inferred_cov <- if (!is.na(scenemode_in_df) && is.character(scenemode_in_df)) {
      ifelse(tolower(scenemode_in_df) == "full", "full", "partial")
    } else NA_character_
    base <- list(
      .table       = table,
      .run_id      = run_id,
      .family      = family,
      .config_tag  = config_tag,
      .is_winner   = is_winner,
      .fold        = as.integer(fold),
      coverage_hint = inferred_cov
    )
    # El 'context' que pase el llamador (budget_gens, budget_pop, used_num_samples, used_bootstrapB,
    c(base, as.list(cfg), context %||% list())
  }
  
  .add_block <- function(biglog, df, table, run_id, family, config_tag,
                         is_winner = FALSE, fold = NA_integer_, context = NULL) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(biglog)
    
    if (identical(table, "final_scenarios") || identical(table, "fold_scenarios")) {
      df <- .ensure_scenario_cols(df)
    } else if (identical(table, "final_convergence") || identical(table, "fold_convergence")) {
      df <- .ensure_convergence_cols(df)
    }
    
    ctx <- .build_context(df, table, run_id, family, config_tag, is_winner, fold, context)
    df  <- .safe_cbind(df, ctx)
    
    biglog[[length(biglog) + 1L]] <- df
    biglog
  }
  
  # Extraer todo de un resultado de CV a filas largas
  # 'context' permite estampar presupuestos/escenario/cobertura/fase:
  #   list(
  #     objective=..., scenario_mode=..., coverage_hint=..., phase="HPF1|HPF2|HALVING|FINAL",
  #     grind_stage=1L, budget_gens=..., budget_pop=..., used_num_samples=..., used_bootstrapB=...
  #   )
  .collect_from_cv <- function(biglog, res_cv, run_id, family, config_tag,
                               is_winner = FALSE, context = NULL) {
    if (is.null(res_cv) || is.null(res_cv$final)) return(biglog)
    
    # Intenta extraer objective/scenario_mode del overall final si no vinieron en 'context'
    ov <- res_cv$final$overall
    if (is.data.frame(ov) && nrow(ov)) {
      context <- modifyList(
        list(
          objective    = ov$objective[1] %||% context$objective,
          scenario_mode= ov$scenario_mode[1] %||% context$scenario_mode
        ),
        context %||% list()
      )
    }
    
    biglog <- .add_block(biglog, res_cv$final$overall,        "final_overall",
                         run_id, family, config_tag, is_winner, fold = NA_integer_, context = context)
    biglog <- .add_block(biglog, res_cv$final$scenario_table, "final_scenarios",
                         run_id, family, config_tag, is_winner, fold = NA_integer_, context = context)
    biglog <- .add_block(biglog, res_cv$final$convergence,    "final_convergence",
                         run_id, family, config_tag, is_winner, fold = NA_integer_, context = context)
    
    fr_list <- res_cv$fold_results %||% list()
    for (k in seq_along(fr_list)) {
      fr <- fr_list[[k]]
      biglog <- .add_block(biglog, fr$overall,        "fold_overall",
                           run_id, family, config_tag, is_winner, fold = k, context = context)
      biglog <- .add_block(biglog, fr$scenario_table, "fold_scenarios",
                           run_id, family, config_tag, is_winner, fold = k, context = context)
      biglog <- .add_block(biglog, fr$convergence,    "fold_convergence",
                           run_id, family, config_tag, is_winner, fold = k, context = context)
    }
    biglog
  }
  

# ====================== EVALUATION: fixed weights (NO GA) ======================

.evaluate_fixed_weights <- function(dist_name, dist_param_grid, weights,
                                    sample_sizes, num_samples,
                                    scenario_mode = c("full","light"),
                                    seed = 1,
                                    q95_B = 300L,
                                    crn_env = NULL,
                                    fam_key = NULL,
                                    # allow optional subset (default NULL = full grid)
                                    scenario_subset = NULL,
                                    subset_tag = NULL,
                                    use_cache = TRUE,
                                    # perturbation at EVAL time (independent from GA perturbation)
                                    do_perturb = TRUE,
                                    pert_sigma = 0.02,
                                    pert_n = 200L,
                                    pert_seed_offset = 7777L) {
  
  scenario_mode <- match.arg(scenario_mode)
  w <- .normalize_simplex(as.numeric(weights))
  
  prepped_full <- prep_scenarios(
    dist_name, dist_param_grid,
    sample_sizes    = sample_sizes,
    num_samples     = num_samples,
    scenario_mode   = scenario_mode,
    seed            = seed,
    crn_env         = crn_env,
    fam_key         = if (is.null(fam_key)) dist_name else fam_key,
    scenario_subset = scenario_subset,
    subset_tag      = subset_tag,
    use_cache       = use_cache
  )
  
  scen_tbl <- dplyr::bind_rows(lapply(prepped_full, function(S) {
    dplyr::bind_rows(lapply(sample_sizes, function(n) {
      summarize_scenario(
        S = S,
        weights = w,
        distribution = dist_name,
        sample_size = n,
        q95_B = as.integer(q95_B),
        q95_seed = seed
      )
    }))
  }))
  
  rob <- dplyr::filter(scen_tbl, estimator == "robust")
  robust_mean_mse <- mean(rob$mse, na.rm = TRUE)
  robust_q95_mse  <- mean(rob$mse_q95, na.rm = TRUE)
  robust_max_mse  <- mean(rob$mse_max, na.rm = TRUE)
  
  base <- dplyr::filter(scen_tbl, estimator != "robust")
  best_by_q95 <- base |>
    dplyr::group_by(estimator) |>
    dplyr::summarise(score = mean(mse_q95, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(score) |>
    dplyr::slice(1)
  
  best_est <- if (nrow(best_by_q95)) best_by_q95$estimator[1] else NA_character_
  
  overall <- tibble::tibble(
    distribution    = dist_name,
    estimator       = best_est,
    objective       = "q95",
    scenario_mode   = scenario_mode,
    robust_mean_mse = robust_mean_mse,
    robust_q95_mse  = robust_q95_mse,
    robust_max_mse  = robust_max_mse
  )
  
  # ---- perturbation at EVAL time (NO GA) ----
  pert <- NULL
  if (isTRUE(do_perturb)) {
    pert <- weight_perturbation_test(
      w_star = w,
      prepped = prepped_full,
      sigma = pert_sigma,
      n = as.integer(pert_n),
      seed = as.integer(seed) + as.integer(pert_seed_offset),
      objective = "q95",
      idx = seq_along(prepped_full),
      use_boot = TRUE,
      B = 0,
      lambda_instab    = 0,
      lambda_entropy   = 0,
      max_w_dominance  = 1,
      lambda_dominance = 0,
      mix_w_q95        = 0.7,
      mix_w_max        = 0.3,
      lambda_bias      = 0,
      use_hierarchical = TRUE,
      group_by         = c("type","rate_bin","scale_bin"),
      group_inner      = c("q95","mean","max"),
      group_inner_q    = 0.95
    )
    
    cat(sprintf(
      "[PERTURB_EVAL] %s | mean=%.6g sd=%.6g q95=%.6g min=%.6g max=%.6g\n",
      dist_name,
      pert$fitness_mean, pert$fitness_sd, pert$fitness_q95,
      pert$fitness_min, pert$fitness_max
    ))
  }
  
  list(
    weights = w,
    estimator_str = weights_to_formula(w),
    overall = overall,
    scenario_table = scen_tbl,
    perturbation = pert
  )
}

# ===== EVALUATION: top-k fixed weights (NO GA) ======================
#This stage is the final one, its a full evaluation outside of the GA, the
# idea is to identity the estimators that perform better in each family under all scenarios.
  # the complete grind runs here. 
.evaluate_topk_fixed_weights <- function(dist_name, dist_param_grid, weights_mat,
                                         sample_sizes, num_samples,
                                         scenario_mode = "full",
                                         seed = 1,
                                         q95_B = 300L,
                                         crn_env = NULL,
                                         fam_key = NULL,
                                         rank_metric = c("robust_q95_mse","robust_mean_mse","robust_max_mse"),
                                         # control perturbation at EVAL time
                                         do_perturb = TRUE,
                                         pert_sigma = 0.02,
                                         pert_n = 200L,
                                         pert_seed_offset = 7777L) {
  
  rank_metric <- match.arg(rank_metric)
  if (is.null(weights_mat) || length(weights_mat) == 0) {
    stop("weights_mat está vacío. No hay top-k para evaluar.")
  }
  weights_mat <- as.matrix(weights_mat)
  
  eval_list <- vector("list", nrow(weights_mat))
  for (k in seq_len(nrow(weights_mat))) {
    eval_list[[k]] <- .evaluate_fixed_weights(
      dist_name = dist_name,
      dist_param_grid = dist_param_grid,
      weights = weights_mat[k, ],
      sample_sizes = sample_sizes,
      num_samples  = num_samples,
      scenario_mode = scenario_mode,
      seed = seed + k,
      q95_B = q95_B,
      crn_env = crn_env,
      fam_key = fam_key,
      scenario_subset = NULL,   # IMPORTANT: NULL => 100% FULL grid
      subset_tag = paste0("EVAL_TOPK__", dist_name),
      use_cache = TRUE,
      do_perturb = do_perturb,
      pert_sigma = pert_sigma,
      pert_n = pert_n,
      pert_seed_offset = pert_seed_offset
    )
  }
  
  # Build ranking table (includes perturbation metrics)
  rank_df <- dplyr::bind_rows(lapply(seq_along(eval_list), function(k) {
    ov <- eval_list[[k]]$overall
    pt <- eval_list[[k]]$perturbation
    
    tibble::tibble(
      candidate_k     = k,
      estimator_str   = eval_list[[k]]$estimator_str,
      robust_mean_mse = .safe_scalar1(ov$robust_mean_mse),
      robust_q95_mse  = .safe_scalar1(ov$robust_q95_mse),
      robust_max_mse  = .safe_scalar1(ov$robust_max_mse),
      
      # perturbation summary (EVAL time)
      pert_mean = pt$fitness_mean %||% NA_real_,
      pert_sd   = pt$fitness_sd   %||% NA_real_,
      pert_q95  = pt$fitness_q95  %||% NA_real_,
      pert_min  = pt$fitness_min  %||% NA_real_,
      pert_max  = pt$fitness_max  %||% NA_real_
    )
  }))
  
  rank_df <- rank_df |>
    dplyr::arrange(.data[[rank_metric]])
  
  
  list(
    rank = rank_df,
    candidates = eval_list
  )
}

evolve_universal_estimator_per_family_cv_two_phase <- function(
    dist_name, dist_param_grid,
    
    # --- PHASE 1 knobs (search) ---
    sample_sizes_light = c(300, 1000),
    num_samples_light  = 30,
    scenario_mode_light = "light",
    
    # --- PHASE 2 knobs (final eval) ---
    sample_sizes_full = c(300, 500, 1000, 2000, 5000),
    num_samples_full  = 80,
    scenario_mode_full = "full",
    
    # --- GA/CV knobs (same as your CV function) ---
    pop_size = 100,
    generations_per_fold = 50,
    seed = 101,
    objective = c("mean","q95","max","mixed"),
    use_parallel = TRUE,
    k_folds = 3L,
    lambda_instab = 0.15,
    bootstrap_B = 200L,
    t_size = 2L,
    elitism = 2L,
    immigrant_rate = 0.10,
    mutation_rate_init = 0.18,
    init_alpha = 1.0,
    alpha_mut = 1.0,
    check_every = 5L,
    patience = 3L,
    min_delta = 0.005,
    final_retrain = TRUE,
    mix_w_q95 = 0.7,
    mix_w_max = 0.3,
    minibatch_frac = NULL,
    minibatch_min = 12L,
    crn_env = NULL,
    fam_key = NULL,
    
    # perturbation controls for PHASE 2 evaluation
    phase2_do_perturb = TRUE,
    phase2_pert_sigma = 0.02,
    phase2_pert_n = 200L,
    phase2_pert_seed_offset = 7777L
) {
  
  # ---------------- PHASE 1: cheap GA search (LIGHT) ----------------
  res_light <- evolve_universal_estimator_per_family_cv(
    dist_name            = dist_name,
    dist_param_grid      = dist_param_grid,
    sample_sizes         = sample_sizes_light,
    num_samples          = num_samples_light,
    scenario_mode        = scenario_mode_light,
    force_scenario_mode  = "light",
    
    pop_size             = pop_size,
    generations_per_fold = generations_per_fold,
    seed                 = seed,
    objective            = objective,
    use_parallel         = use_parallel,
    k_folds              = k_folds,
    lambda_instab        = lambda_instab,
    bootstrap_B          = bootstrap_B,
    t_size               = t_size,
    elitism              = elitism,
    immigrant_rate       = immigrant_rate,
    mutation_rate_init   = mutation_rate_init,
    init_alpha           = init_alpha,
    alpha_mut            = alpha_mut,
    check_every          = check_every,
    patience             = patience,
    min_delta            = min_delta,
    final_retrain        = final_retrain,
    mix_w_q95            = mix_w_q95,
    mix_w_max            = mix_w_max,
    minibatch_frac       = minibatch_frac,
    minibatch_min        = minibatch_min,
    crn_env              = crn_env,
    fam_key              = if (is.null(fam_key)) dist_name else fam_key
  )
  
  w_star <- res_light$final$weights
  
  # ---------------- PHASE 2: evaluate weights on FULL grid (no GA) ----------------
  eval_full <- .evaluate_fixed_weights(
    dist_name       = dist_name,
    dist_param_grid = dist_param_grid,
    weights         = w_star,
    sample_sizes    = sample_sizes_full,
    num_samples     = num_samples_full,
    scenario_mode   = scenario_mode_full,
    seed            = seed,
    q95_B           = bootstrap_B,
    crn_env         = crn_env,
    fam_key         = if (is.null(fam_key)) dist_name else fam_key,
    do_perturb      = phase2_do_perturb,
    pert_sigma      = phase2_pert_sigma,
    pert_n          = phase2_pert_n,
    pert_seed_offset= phase2_pert_seed_offset
  )
  
  # return in the same "shape" 
  list(
    fold_results = res_light$fold_results,
    final = list(
      weights = eval_full$weights,
      overall = eval_full$overall,
      scenario_table = eval_full$scenario_table,
      convergence = res_light$final$convergence,
      perturbation = eval_full$perturbation,
      best_val_score = eval_full$overall$robust_q95_mse[1]
    ),
    phase1_light = res_light
  )
}


#Trigger for the module run tracker
mark_module_done("07_fitness_objectives.R")
