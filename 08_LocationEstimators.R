
# =============================================================================
# 08_LOCATIONESTIMATORS — RESEARCH COMMENTED EDITION
# =============================================================================
# The purpose of this module is to serve a specific role within the simulation,
# optimization, and evaluation pipeline for robust location estimators under
# contamination. Comments focus on experimental design, reproducibility, and
# statistical meaning rather than user instructions.
# =============================================================================

this_dir = setwd("C:/Users/AlvaroRivera-Eraso/Documents/simulating distributions and estimators/LocationEstimators11_modular_v1")



# ====================== AUTO-ENABLE DEBUG FLAGS ======================
# Turn on verbose debug logs for MAIN + TopK diagnostics
Sys.setenv(DEBUG_MAIN = "1")

# ======================================================================

catf("\n[DEFENSIVE] this_dir = %s\n", this_dir)
catf("[DEFENSIVE] getwd()  = %s\n\n", getwd())

# 1) Confirm module files exist where you think they do
module_files <- c(
  "00_utils_debug_io.R",
  "01_paths_repro.R",
  "02_scenarios_sampling.R",
  "03_distributions_params.R",
  "04_estimators_registry.R",
  "05_ga_core.R",
  "06_data_prep.R",
  "07_fitness_objectives.R"
)

missing_files <- module_files[!file.exists(file.path(this_dir, module_files))]
if (length(missing_files) > 0L) {
  msg <- paste0(
    "[DEFENSIVE] Missing module files under this_dir:\n  - ",
    paste(missing_files, collapse = "\n  - "),
    "\n\nFix: check this_dir / working directory, or rename module_files list."
  )
  stop(msg, call. = FALSE)
} else {
  catf("[DEFENSIVE] All module files exist under this_dir.\n")
}

# 2) Confirm critical symbols are loaded into the environment
#    (1–3 per module is enough; tweak as you evolve.)
required <- list(
  "00_utils" = list(
    fun = c("safe_write_csv", "safe_save_rds", "safe_write_lines"),
    obj = c()
  ),
  "01_paths_repro" = list(
    fun = c("make_crn_indices", ".ensure_seed", ".seed_scope"),
    obj = c("IS_WINDOWS", "OUT_ROOT")
  ),
  "02_scenarios_sampling" = list(
    fun = c("build_scenarios_full", "build_scenarios_light", "pick_scenario_subset"),
    obj = c(".scenario_subset_cache")
  ),
  "03_distributions_params" = list(
    fun = c("generate_population", "analytic_mean_from_params"),
    obj = c("param_grids")
  ),
  "04_estimators_registry" = list(
    fun = c("custom_estimator", "weights_to_formula"),
    obj = c("ESTIMATOR_REGISTRY", "ESTIMATOR_NAMES", "N_EST")
  ),
  "05_ga_core" = list(
    fun = c("init_population", "crossover", "mutate_weights", "tournament_select", "make_folds"),
    obj = c(".ga_default_ctrl")
  ),
  "06_data_prep" = list(
    fun = c("build_scenarios", "prep_scenarios", "inject_outliers_realistic"),
    obj = c(".prepped_cache")
  ),
  "07_fitness_objectives" = list(
    fun = c("fitness_universal", "random_search_baseline"),
    obj = c()
  )
)

missing_syms <- character(0)

for (mod in names(required)) {
  # functions
  for (fn in required[[mod]]$fun) {
    if (!exists(fn, mode = "function", inherits = TRUE)) {
      missing_syms <- c(missing_syms, paste0(mod, "::", fn, " (function)"))
    }
  }
  # objects
  for (ob in required[[mod]]$obj) {
    if (!exists(ob, inherits = TRUE)) {
      missing_syms <- c(missing_syms, paste0(mod, "::", ob, " (object)"))
    }
  }
}

if (length(missing_syms) > 0L) {
  msg <- paste0(
    "[DEFENSIVE] Module load check FAILED. Missing symbols:\n  - ",
    paste(missing_syms, collapse = "\n  - "),
    "\n\nMost common cause: source() paths are wrong (this_dir/getwd mismatch).\n",
    "Fix: ensure the MAIN file is sourced from its own folder, or setwd(this_dir) before source()."
  )
  stop(msg, call. = FALSE)
} else {
  catf("[DEFENSIVE] Module load check OK: all required symbols are present.\n\n")
}
# ===============================================================================




# ====================== DEBUG SWITCH (MAIN) ======================
# Set environment variable DEBUG_MAIN=1 to enable verbose console logs from MAIN pipeline.
.dbg_main <- function(msg) {
  if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
    message(sprintf("[DEBUG_MAIN %s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  }
}
# Convenience: print short object summary (names, dim, class)
.dbg_obj <- function(name, x) {
  if (!identical(Sys.getenv("DEBUG_MAIN","0"), "1")) return(invisible(NULL))
  cls <- paste(class(x), collapse=",")
  dn  <- tryCatch(paste(dim(x), collapse="x"), error=function(e) NA_character_)
  nms <- tryCatch(paste(utils::head(names(x), 20), collapse=","), error=function(e) "<no-names>")
  message(sprintf("[DEBUG_MAIN %s] %s | class=%s | dim=%s | names(head)=%s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), name, cls, dn, nms))
  invisible(NULL)
}

# ===== Original launcher / main code =====
# ====================== LAUNCHER: TODO EN UNA TIRADA ======================

run_all_one_shot <- function(
    families_to_run = c(
      "normal",
      "lognormal",
      "weibull",
      "exgaussian",
      "exwald"
    ),
    pop_sizes       = c(80, 100),
    mutation_rates  = c(0.12, 0.18, 0.24),
    init_alphas     = c(0.5, 1.0),
    alpha_mut_set   = c(0.5, 1.0),
    immigrant_rates = c(0.05, 0.10),
    t_sizes         = c(2L, 3L),
    elitism_set     = c(1L, 2L),
    seeds           = c(101, 202),
    
    # GA/scenarios 
    sample_sizes = c(300, 500, 1000, 2000, 5000),
    num_samples  = 80,
    generations_per_fold = 60,
    k_folds = 3,
    final_retrain = TRUE,
    objective    = "q95",
    lambda_instab_default = 0.15,
    bootstrap_B  = 200L,
    check_every  = 5L,
    patience     = 3L,
    min_delta    = 0.005,
    use_parallel = TRUE,
    mix_w_q95 = 0.7, mix_w_max = 0.3,
    out_root = OUT_ROOT,
    minibatch_frac = 0.5,
    minibatch_min  = 8L,
    pick_single_winner_per_family = TRUE,
    suppress_intermediate_saves   = FALSE,
    write_all_in_one = FALSE,
    winner_metric = c("robust_q95_mse","robust_mean_mse","robust_max_mse"),
    
    search_strategy = c("grid","random","halving"),
    random_n = 40,
    eta = 3,                   # halving reduction factor
    finalists_per_family = 2,  # cuántos empujar a “full grind”
    
    # >>> Schedules para HALVING (multi-fidelity, FASE 3)
    # Early stages = menor presupuesto y menos muestras para velocidad
    stage_gens_frac = c(0.25, 0.60, 1.00),
    stage_pop_size  = c(40,   60,   NA),  # NA -> usa max(pop_sizes)
    stage_num_samples = c(30,  60,   NA), # NA -> usa full num_samples
    stage_bootstrap_B = c(60L, 100L, NA), # NA -> usa full bootstrap_B
    
    # FASES ESTRICTAS DEL GRIND DE HP (FASE 1 y FASE 2)
    hp_strict = TRUE,            # activa Fase 1 + Fase 2 strict
    # Fase 1 (minibatch de configs + minibatch de escenarios)
    hp1_frac = 0.25,             # o usa hp1_k exacto
    hp1_k    = NULL,
    hp1_gens = 10L,
    hp1_num_samples = NULL,      # por defecto: round(0.4 * num_samples)
    hp1_bootstrap_B = NULL,      # por defecto: round(0.33 * bootstrap_B)
    hp1_minibatch_frac = 0.50,   # minibatch de escenarios (solo si modo=light)
    hp1_minibatch_min  = 10L,
    hp2_k    = 20L,              # top-K exacto tras F1
    hp2_gens = 20L,
    hp2_num_samples = NULL,      # por defecto: round(0.6 * num_samples)
    hp2_bootstrap_B = NULL,      # por defecto: round(0.6 * bootstrap_B)
    hp2_minibatch_frac = 0.50,
    hp2_minibatch_min  = 12L,
    
    # scenario coverage per phase (subsampling)
    hp1_scenario_frac   = 0.20,
    hp2_scenario_frac   = 0.40,
    stage_scenario_frac = c(0.60, 0.80, 1.00),
    # ensure fair comparison (same scenario subset across configs in a phase/stage)
    fair_scenario_subsets = TRUE,
    
    # Pool que pasa a Fase 3 (HALVING)
    hp3_pool_k = 10L          
) {
  # Accumulate per-family winners (data frames) to return at end
  winners_acc <- list()

  run_id  <- paste0("GA_UNIFIED_", RUN_TS())
  root_out <- file.path(out_root, run_id); mkdirp(root_out)
  
  # acumulador para la mega-CSV ===
  biglog <- list() #This is not in use anymore, I can take this out. 
  
  # --- MANIFEST enriquecido ---
  search_strategy <- match.arg(search_strategy)
  winner_metric   <- match.arg(winner_metric)
  # -- scientific scenario universe exports ---
  index_init()
  run_dir <- root_out
  sc_universe_full  <- write_scenarios_universe(run_dir, scenario_mode = "full")
  sc_universe_light <- write_scenarios_universe(run_dir, scenario_mode = "light")
  meta <- list(
    run_id = run_id,
    families = families_to_run,
    k_folds = k_folds,
    objective = objective,
    seeds = paste(seeds, collapse = ","),
    pop_sizes = pop_sizes,
    generations_per_fold = generations_per_fold,
    final_retrain = final_retrain,
    bootstrap_B = bootstrap_B,
    lambda_instab_default = lambda_instab_default,
    mutation_rates = mutation_rates,
    init_alphas = init_alphas,
    alpha_mut_set = alpha_mut_set,
    immigrant_rates = immigrant_rates,
    t_sizes = t_sizes,
    elitism_set = elitism_set,
    minibatch_frac = minibatch_frac,
    minibatch_min = minibatch_min,
    scenario_mode = "both",
    search_strategy = search_strategy,
    random_n = random_n,
    eta = eta,
    finalists_per_family = finalists_per_family,
    stage_gens_frac = stage_gens_frac,
    stage_pop_size  = stage_pop_size,
    stage_num_samples = stage_num_samples,
    stage_bootstrap_B = stage_bootstrap_B,
    winner_metric = winner_metric,
    pick_single_winner_per_family = pick_single_winner_per_family,
    # HP strict manifest
    hp_strict = hp_strict,
    hp1_frac = hp1_frac, hp1_k = hp1_k, hp1_gens = hp1_gens,
    hp1_num_samples = hp1_num_samples, hp1_bootstrap_B = hp1_bootstrap_B,
    hp1_minibatch_frac = hp1_minibatch_frac, hp1_minibatch_min = hp1_minibatch_min,
    hp2_k = hp2_k, hp2_gens = hp2_gens,
    hp2_num_samples = hp2_num_samples, hp2_bootstrap_B = hp2_bootstrap_B,
    hp2_minibatch_frac = hp2_minibatch_frac, hp2_minibatch_min = hp2_minibatch_min,
    hp3_pool_k = hp3_pool_k
  )
  # # Bind winners across families (if any)
  # if (length(winners_acc) > 0L) {
  #   winners_df <- dplyr::bind_rows(winners_acc)
  # } else {
  #   winners_df <- NULL
  # } # useless, need to be taken out *****


  write_manifest(root_out, meta)
  hp1_gens <- max(5L, as.integer(hp1_gens))
  hp2_gens <- max(5L, as.integer(hp2_gens))
  
  
  
  # ---- CRN ----
  crn_env <- make_crn_indices(
    families     = families_to_run,
    sample_sizes = sample_sizes,
    B_boot       = bootstrap_B,
    num_samples  = num_samples
  )
  
  cfg_grid <- expand.grid(
    pop_size           = pop_sizes,
    mutation_rate_init = mutation_rates,
    init_alpha         = init_alphas,
    alpha_mut          = alpha_mut_set,
    immigrant_rate     = immigrant_rates,
    t_size             = t_sizes,
    elitism            = elitism_set,
    seed               = seeds,
    lambda_instab      = lambda_instab_default,
    stringsAsFactors   = FALSE
  )
  
  # ---- Planificador de pasos para mostrar 1-100% --------------------------
  .plan_family_steps <- function(n_cfg) {
    #if (search_strategy != "halving") return(as.integer(n_cfg))  # 1 eval por config
    
    # F1/F2 si hp_strict
    n1 <- if (isTRUE(hp_strict)) {
      if (!is.null(hp1_k)) min(hp1_k, n_cfg) else ceiling(n_cfg * hp1_frac)
    } else 0
    n1 <- as.integer(n1)
    
    n2 <- if (isTRUE(hp_strict)) {
      if (!is.null(hp2_k)) min(hp2_k, n1) else ceiling(n1 * 0.33)
    } else n_cfg
    n2 <- as.integer(n2)
    
    start_halv <- if (isTRUE(hp_strict)) min(hp3_pool_k, n2) else n2
    start_halv <- as.integer(start_halv)
    
    S <- length(stage_gens_frac)
    if (length(stage_scenario_frac) != S) {
      stop(sprintf("stage_scenario_frac debe tener longitud %d (como stage_gens_frac).", S))
    }
    cur <- start_halv
    n_halv <- 0L
    for (s in seq_len(S)) {
      n_halv <- n_halv + max(0L, as.integer(cur))
      cur <- if (s < S) max(1L, as.integer(ceiling(cur / eta))) else as.integer(finalists_per_family)
    }
    n_final <- as.integer(finalists_per_family)
    
    as.integer(n1 + n2 + n_halv + n_final)
  }
  
  
  total_steps <- sum(vapply(families_to_run, function(fam) .plan_family_steps(nrow(cfg_grid)), integer(1)))
  progress_init(total_steps)
  
  # ======== util interno: score de etapa (q95 val con fallback) =========
  .stage_score <- function(res_cv) {
    s <- res_cv$final$best_val_score
    if (is.null(s) || length(s) != 1L || !is.finite(s)) {
      ov <- res_cv$final$overall
      s <- .safe_scalar1(ov$robust_q95_mse)
    }
    if (is.null(s) || length(s) != 1L || !is.finite(s)) s <- Inf
    s
  }
  .safe_scalar1 <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) < 1L) return(default)
    v <- suppressWarnings(as.numeric(x[[1]]))
    if (!is.finite(v)) return(default)
    v
  }

  # ======== INTERNAL: CONFIG TAG FALLBACK (TRACEABILITY) ========================
  # Some pipelines keep the hyperparameter pool ("fam_pool") as a pure numeric grid
  # without carrying the original human-readable tag column (e.g., HPF1__/HPF2__).
  # When that happens, downstream winner-writing can lose traceability (tag becomes NA).
  # This helper rebuilds a deterministic tag from the config row *without affecting*
  # any scoring or optimization logic. It is used only for naming/logging outputs.
  .make_cfg_tag_fallback <- function(cfg, prefix = "CFG") {
    # Defensive: accept either data.frame row or list-like object.
    ps <- as.integer(cfg$pop_size)
    mr <- as.numeric(cfg$mutation_rate_init)
    a  <- as.numeric(cfg$init_alpha)
    am <- as.numeric(cfg$alpha_mut)
    im <- as.numeric(cfg$immigrant_rate)
    tt <- as.integer(cfg$t_size)
    ee <- as.integer(cfg$elitism)
    li <- as.numeric(cfg$lambda_instab)
    sd <- as.integer(cfg$seed)
    sprintf("%s__ps%d_mr%.2f_a%.2f_am%.2f_im%.2f_t%d_e%d_li%.2f_seed%d",
            prefix, ps, mr, a, am, im, tt, ee, li, sd)
  }
  # ============================================================================
  if (search_strategy == "random") {
    set.seed(13)
    if (nrow(cfg_grid) > random_n) {
      cfg_grid <- cfg_grid[sample(seq_len(nrow(cfg_grid)), random_n), , drop = FALSE]
    }
  }
  if (search_strategy == "halving") {
    S <- length(stage_gens_frac)
    if (length(stage_scenario_frac) != S) {
      stop(sprintf("stage_scenario_frac debe tener longitud %d (como stage_gens_frac).", S))
    }
    stopifnot(length(stage_pop_size)   == S,
              length(stage_num_samples)== S,
              length(stage_bootstrap_B)== S)
  }
  
  # ======== Acumuladores =================================================
  big_rows <- list()   # para incluir el rollup dentro de la mega-CSV

  # winners_acc initialized at function start
  # ============================== LOOP FAMILIAS ==========================
  for (dist in families_to_run) {
    message(sprintf("\n========== FAMILY: %s ==========\n", dist))
    grid <- param_grids[[dist]]
    fam_pool <- cfg_grid  # candidatos iniciales
    fam_dir <- file.path(root_out, dist)
    mkdirp(fam_dir)
    
    # Universe con IDs (para trazabilidad)
    sc_u_light <- add_scenario_ids(build_scenarios("light"))
    sc_u_full  <- add_scenario_ids(build_scenarios("full"))
    
    # ===================== FASE 1: SIMPLE (LIGHT) ===================
    if (isTRUE(hp_strict)) {
      hp1_ns <- hp1_num_samples %||% max(1L, round(0.4 * num_samples))
      hp1_B  <- hp1_bootstrap_B %||% max(20L, round(0.33 * bootstrap_B))
      
      fam_pool_phase1 <- if (!is.null(hp1_k)) {
        pick_minibatch_configs(fam_pool, k = hp1_k, seed = 13,
                               stratify_cols = c("pop_size","t_size"))
      } else {
        pick_minibatch_configs(fam_pool, frac = hp1_frac, seed = 13,
                               stratify_cols = c("pop_size","t_size"))
      }
      
      message(sprintf("[HP-F1 | SIMPLE] %s | configs=%d | gens=%d | ns=%d | B=%d",
                      dist, nrow(fam_pool_phase1), hp1_gens, hp1_ns, hp1_B))
      
      sc_seed_hpf1 <- if (isTRUE(fair_scenario_subsets)) .stable_int_seed(paste0(dist,"::HPF1"), base_seed = 10001L) else 13L
      sc_used_hpf1 <- pick_scenario_subset(sc_u_light, scenario_frac = hp1_scenario_frac, min_n = k_folds, seed = sc_seed_hpf1,
                                           subset_tag = paste0(dist, "::HPF1"))
      safe_write_csv(sc_used_hpf1, file.path(fam_dir, sprintf("scenarios_used_%s_HPF1.csv", dist)))
      safe_write_csv(stage_difficulty(sc_used_hpf1), file.path(fam_dir, sprintf("difficulty_%s_HPF1.csv", dist)))
      index_add(dist, "HPF1", file.path(fam_dir, sprintf("scenarios_used_%s_HPF1.csv", dist)), n_candidates = nrow(fam_pool_phase1), n_scenarios = nrow(sc_used_hpf1), seed = sc_seed_hpf1)
      index_add(dist, "HPF1", file.path(fam_dir, sprintf("difficulty_%s_HPF1.csv", dist)), n_candidates = nrow(fam_pool_phase1), n_scenarios = nrow(sc_used_hpf1), seed = sc_seed_hpf1)
      
      phase1_rows <- vector("list", nrow(fam_pool_phase1))
      
      for (i in seq_len(nrow(fam_pool_phase1))) {
        cfg <- fam_pool_phase1[i, ]
        tag <- sprintf("HPF1__ps%d_mr%.2f_a%.2f_am%.2f_im%.2f_t%d_e%d_li%.2f_B%d_ns%d_sf%.2f_sc%d_seed%d",
                       cfg$pop_size, cfg$mutation_rate_init, cfg$init_alpha, cfg$alpha_mut,
                       cfg$immigrant_rate, cfg$t_size, cfg$elitism, cfg$lambda_instab,
                       hp1_B, hp1_ns, hp1_scenario_frac, sc_seed_hpf1, cfg$seed)
        
        res1 <- do.call(
          evolve_universal_estimator_per_family_cv,
          list(
            dist_name            = dist,
            dist_param_grid      = grid,
            pop_size             = cfg$pop_size,
            generations_per_fold = hp1_gens,
            seed                 = cfg$seed,
            objective            = objective,
            use_parallel         = use_parallel,
            k_folds              = k_folds,
            lambda_instab        = cfg$lambda_instab,
            bootstrap_B          = hp1_B,
            t_size               = cfg$t_size,
            elitism              = cfg$elitism,
            immigrant_rate       = cfg$immigrant_rate,
            mutation_rate_init   = cfg$mutation_rate_init,
            init_alpha           = cfg$init_alpha,
            alpha_mut            = cfg$alpha_mut,
            check_every          = check_every,
            patience             = patience,
            min_delta            = min_delta,
            final_retrain        = FALSE,
            mix_w_q95            = mix_w_q95,
            mix_w_max            = mix_w_max,
            minibatch_frac       = hp1_minibatch_frac,
            minibatch_min        = hp1_minibatch_min,
            crn_env              = crn_env,
            fam_key              = dist,
            sample_sizes         = c(300, 1000),
            num_samples          = hp1_ns,
            force_scenario_mode  = "light",
            scenario_frac        = hp1_scenario_frac,
            scenario_seed        = if (isTRUE(fair_scenario_subsets)) .stable_int_seed(paste0(dist,"::HPF1"), base_seed = 10001L) else cfg$seed,
            subset_tag           = paste0("HPF1__", dist),
            scenario_universe    = sc_u_light,
            scenario_subset_override = sc_used_hpf1,
            run_dir              = fam_dir,
            stage_name           = "HPF1"
          )
        )
        
        progress_step(sprintf("HPF1 %s | %s", dist, tag))
        
        ov <- res1$final$overall
        sc <- .stage_score(res1)
        
        phase1_rows[[i]] <- data.frame(
          family = dist, config_tag = tag, stage = 1L,
          is_finalist = FALSE, score_stage = sc,
          robust_mean_mse_stage = .safe_scalar1(ov$robust_mean_mse),
          robust_q95_mse_stage  = .safe_scalar1(ov$robust_q95_mse),
          robust_max_mse_stage  = .safe_scalar1(ov$robust_max_mse),
          stringsAsFactors = FALSE
        )
      }
      
      big_rows <- c(big_rows, phase1_rows)
      
      df1  <- dplyr::bind_rows(phase1_rows)
      ord1 <- order(df1$score_stage)
      keep1 <- max(1L, min(hp2_k, nrow(df1)))
      fam_pool <- fam_pool_phase1[ord1[seq_len(keep1)], , drop = FALSE]

      .dbg_main(sprintf("POST-HPF1 pool built | %s | keep1=%d | fam_pool_n=%d | has_cfg_tag_col=%s",
                        dist, keep1, nrow(fam_pool), ("cfg_tag" %in% names(fam_pool))))
      if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
        .dbg_main(sprintf("POST-HPF1 fam_pool colnames: %s", paste(names(fam_pool), collapse=", ")))
        .dbg_main(sprintf("POST-HPF1 df1 tag sample: %s", paste(utils::head(df1$config_tag, 3), collapse=" || ")))
      }
    }
    
    # ===================== FASE 2: INTERMEDIA (LIGHT+) ===================
    if (isTRUE(hp_strict)) {
      hp2_ns <- hp2_num_samples %||% max(1L, round(0.6 * num_samples))
      hp2_B  <- hp2_bootstrap_B %||% max(40L, round(0.6 * bootstrap_B))
      
      message(sprintf("[HP-F2 | INTERMEDIA] %s | configs=%d | gens=%d | ns=%d | B=%d",
                      dist, nrow(fam_pool), hp2_gens, hp2_ns, hp2_B))
      
      sc_seed_hpf2 <- if (isTRUE(fair_scenario_subsets)) .stable_int_seed(paste0(dist,"::HPF2"), base_seed = 20002L) else 13L
      sc_used_hpf2 <- pick_scenario_subset(sc_u_full, scenario_frac = hp2_scenario_frac, min_n = k_folds, seed = sc_seed_hpf2,
                                           subset_tag = paste0(dist, "::HPF2"))
      safe_write_csv(sc_used_hpf2, file.path(fam_dir, sprintf("scenarios_used_%s_HPF2.csv", dist)))
      safe_write_csv(stage_difficulty(sc_used_hpf2), file.path(fam_dir, sprintf("difficulty_%s_HPF2.csv", dist)))
      index_add(dist, "HPF2", file.path(fam_dir, sprintf("scenarios_used_%s_HPF2.csv", dist)), n_candidates = nrow(fam_pool), n_scenarios = nrow(sc_used_hpf2), seed = sc_seed_hpf2)
      index_add(dist, "HPF2", file.path(fam_dir, sprintf("difficulty_%s_HPF2.csv", dist)), n_candidates = nrow(fam_pool), n_scenarios = nrow(sc_used_hpf2), seed = sc_seed_hpf2)
      
      phase2_rows <- vector("list", nrow(fam_pool))
      phase2_res_list  <- vector("list", nrow(fam_pool))
      phase2_topk_list <- vector("list", nrow(fam_pool))
      
      for (i in seq_len(nrow(fam_pool))) {
        cfg <- fam_pool[i, ]
        tag <- sprintf("HPF2__ps%d_mr%.2f_a%.2f_am%.2f_im%.2f_t%d_e%d_li%.2f_B%d_ns%d_sf%.2f_sc%d_seed%d",
                       cfg$pop_size, cfg$mutation_rate_init, cfg$init_alpha, cfg$alpha_mut,
                       cfg$immigrant_rate, cfg$t_size, cfg$elitism, cfg$lambda_instab,
                       hp2_B, hp2_ns, hp2_scenario_frac, sc_seed_hpf2, cfg$seed)
        
        res2 <- do.call(
          evolve_universal_estimator_per_family_cv,
          list(
            dist_name            = dist,
            dist_param_grid      = grid,
            pop_size             = cfg$pop_size,
            generations_per_fold = hp2_gens,
            seed                 = cfg$seed,
            objective            = objective,
            use_parallel         = use_parallel,
            k_folds              = k_folds,
            lambda_instab        = cfg$lambda_instab,
            bootstrap_B          = hp2_B,
            t_size               = cfg$t_size,
            elitism              = cfg$elitism,
            immigrant_rate       = cfg$immigrant_rate,
            mutation_rate_init   = cfg$mutation_rate_init,
            init_alpha           = cfg$init_alpha,
            alpha_mut            = cfg$alpha_mut,
            check_every          = check_every,
            patience             = patience,
            min_delta            = min_delta,
            final_retrain        = FALSE,
            mix_w_q95            = mix_w_q95,
            mix_w_max            = mix_w_max,
            minibatch_frac       = hp2_minibatch_frac,
            minibatch_min        = hp2_minibatch_min,
            crn_env              = crn_env,
            fam_key              = dist,
            sample_sizes         = c(300, 500, 1000),
            num_samples          = hp2_ns,
            force_scenario_mode  = "full",
            scenario_frac        = hp2_scenario_frac,
            scenario_seed        = if (isTRUE(fair_scenario_subsets)) .stable_int_seed(paste0(dist,"::HPF2"), base_seed = 20002L) else cfg$seed,
            subset_tag           = paste0("HPF2__", dist),
            scenario_universe    = sc_u_full,
            scenario_subset_override = sc_used_hpf2,
            run_dir              = fam_dir,
            stage_name           = "HPF2"
          )
        )
        
        progress_step(sprintf("HPF2 %s | %s", dist, tag))
        
        phase2_res_list[[i]]  <- res2
        
        phase2_topk_list[[i]] <- .extract_topk_from_res(res2, K = finalists_per_family)
        
        ov <- res2$final$overall
        sc <- .stage_score(res2)
        
        phase2_rows[[i]] <- data.frame(
          family = dist, config_tag = tag, stage = 2L,
          is_finalist = FALSE, score_stage = sc,
          robust_mean_mse_stage = .safe_scalar1(ov$robust_mean_mse),
          robust_q95_mse_stage  = .safe_scalar1(ov$robust_q95_mse),
          robust_max_mse_stage  = .safe_scalar1(ov$robust_max_mse),
          stringsAsFactors = FALSE
        )
      }
      
      big_rows <- c(big_rows, phase2_rows)
      
      df2  <- dplyr::bind_rows(phase2_rows)
      ord2 <- order(df2$score_stage)
      sel2 <- ord2[seq_len(min(hp3_pool_k, nrow(df2)))]
      fam_pool <- fam_pool[sel2, , drop = FALSE]

      .dbg_main(sprintf("POST-HPF2 pool built | %s | sel2_n=%d | fam_pool_n=%d | has_cfg_tag_col=%s",
                        dist, length(sel2), nrow(fam_pool), ("cfg_tag" %in% names(fam_pool))))
      if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
        .dbg_main(sprintf("POST-HPF2 fam_pool colnames: %s", paste(names(fam_pool), collapse=", ")))
        .dbg_main(sprintf("POST-HPF2 df2 tag sample: %s", paste(utils::head(df2$config_tag, 3), collapse=" || ")))
      }
      # inicializar herencia para HALVING desde HPF2 (res + topk)
      fam_pool_res  <- phase2_res_list[sel2]
      fam_pool_topk <- phase2_topk_list[sel2]
    }
    
    # ===================== FASE 3: HALVING (FULL) ===================
    S <- length(stage_gens_frac)
    if (length(stage_scenario_frac) != S) {
      stop(sprintf("stage_scenario_frac debe tener longitud %d (como stage_gens_frac).", S))
    }
    
    
# IMPORTANT: fam_pool_res / fam_pool_topk se usan para heredar fallback entre stages.
.dbg(sprintf("HALVING init | %s | incoming pool=%d | have_inherited_res=%s | have_inherited_topk=%s",
             dist, nrow(fam_pool),
             (!is.null(fam_pool_res) && length(fam_pool_res) == nrow(fam_pool)),
             (!is.null(fam_pool_topk) && length(fam_pool_topk) == nrow(fam_pool))))

    for (s in seq_len(S)) {
      gens_s <- max(30L, ceiling(generations_per_fold * stage_gens_frac[s]))
      pop_s  <- if (is.na(stage_pop_size[s])) max(pop_sizes) else stage_pop_size[s]
      num_s  <- if (is.na(stage_num_samples[s])) num_samples else stage_num_samples[s]
      B_s    <- if (is.na(stage_bootstrap_B[s])) bootstrap_B else stage_bootstrap_B[s]
      
      sc_seed_stage <- if (isTRUE(fair_scenario_subsets)) .stable_int_seed(paste0(dist,"::STAGE",s), base_seed = 30003L) else 13L
      sc_used_stage <- pick_scenario_subset(sc_u_full, scenario_frac = stage_scenario_frac[s], min_n = k_folds, seed = sc_seed_stage,
                                            subset_tag = paste0(dist, "::HALVING::", s))
      stname <- paste0("HALVING_s", s)
      safe_write_csv(sc_used_stage, file.path(fam_dir, sprintf("scenarios_used_%s_%s.csv", dist, stname)))
      safe_write_csv(stage_difficulty(sc_used_stage), file.path(fam_dir, sprintf("difficulty_%s_%s.csv", dist, stname)))
      index_add(dist, stname, file.path(fam_dir, sprintf("scenarios_used_%s_%s.csv", dist, stname)), n_candidates = nrow(fam_pool), n_scenarios = nrow(sc_used_stage), seed = sc_seed_stage)
      index_add(dist, stname, file.path(fam_dir, sprintf("difficulty_%s_%s.csv", dist, stname)), n_candidates = nrow(fam_pool), n_scenarios = nrow(sc_used_stage), seed = sc_seed_stage)
      
      # progressive sample_sizes (reduce hot-loop early)
      ss_s <- if (s == 1L) {
        intersect(sample_sizes, c(300, 1000))
      } else if (s == 2L) {
        intersect(sample_sizes, c(300, 500, 1000, 2000))
      } else {
        sample_sizes
      }
      if (!length(ss_s)) ss_s <- sample_sizes
      
      stage_eval <- vector("list", nrow(fam_pool))
      
      for (i in seq_len(nrow(fam_pool))) {
        cfg <- fam_pool[i, ]
        
        res_cv <- tryCatch(
          do.call(
            evolve_universal_estimator_per_family_cv,
          list(
            dist_name            = dist,
            dist_param_grid      = grid,
            pop_size             = pop_s,
            generations_per_fold = gens_s,
            seed                 = cfg$seed,
            objective            = objective,
            use_parallel         = use_parallel,
            k_folds              = k_folds,
            lambda_instab        = cfg$lambda_instab,
            bootstrap_B          = B_s,
            t_size               = cfg$t_size,
            elitism              = cfg$elitism,
            immigrant_rate       = cfg$immigrant_rate,
            mutation_rate_init   = cfg$mutation_rate_init,
            init_alpha           = cfg$init_alpha,
            alpha_mut            = cfg$alpha_mut,
            final_retrain        = FALSE,
            crn_env              = crn_env,
            fam_key              = dist,
            sample_sizes         = ss_s,
            num_samples          = num_s,
            force_scenario_mode  = "full",
            scenario_frac        = stage_scenario_frac[s],
            scenario_seed        = if (isTRUE(fair_scenario_subsets))
              .stable_int_seed(paste0(dist,"::STAGE",s), base_seed = 30003L)
            else cfg$seed,
            subset_tag           = paste0("HALVING__", dist, "__s", s),
            scenario_universe    = sc_u_full,
            scenario_subset_override = sc_used_stage,
            run_dir              = fam_dir,
            stage_name           = stname
          )
        )
        , error = function(e) {
          
message(sprintf("[HALVING] %s | stage %d/%d | cfg %d/%d | ERROR: %s",
                dist, s, S, i, nrow(fam_pool), conditionMessage(e)))
.dbg(sprintf("HALVING error stack | %s | stage %d/%d | cfg %d/%d | tag=%s | stack=%s",
             dist, s, S, i, nrow(fam_pool),
             if (!is.null(cfg$tag)) cfg$tag else "<no-tag>",
             .debug_stack()))
NULL
        }
        )

        topk_payload <- .extract_topk_from_res(res_cv, K = finalists_per_family)
        if (is.null(topk_payload) || length(topk_payload) == 0L) {
          .debug_topk_diagnose(res_cv, label = sprintf("HALVING topk missing (pre-fallback) | %s | stage %d/%d | cfg %d/%d | ", dist, s, S, i, nrow(fam_pool)))
        }
        if (is.null(topk_payload) || length(topk_payload) == 0L) {
          if (!is.null(fam_pool_topk) && length(fam_pool_topk) >= i &&
              !is.null(fam_pool_topk[[i]]) && length(fam_pool_topk[[i]]) > 0L) {
            topk_payload <- fam_pool_topk[[i]]
          } else if (!is.null(fam_pool_res) && length(fam_pool_res) >= i) {
            topk_payload <- .extract_topk_from_res(fam_pool_res[[i]], K = finalists_per_family)
          }
        }

if (is.null(topk_payload) || length(topk_payload) == 0L) {
  .dbg(sprintf("HALVING topk still missing (post-fallback) | %s | stage %d/%d | cfg %d/%d | tag=%s | have_prev_topk=%s | have_prev_res=%s",
               dist, s, S, i, nrow(fam_pool),
               if (!is.null(cfg$tag)) cfg$tag else "<no-tag>",
               (!is.null(fam_pool_topk) && length(fam_pool_topk) >= i && !is.null(fam_pool_topk[[i]]) && length(fam_pool_topk[[i]]) > 0L),
               (!is.null(fam_pool_res) && length(fam_pool_res) >= i && !is.null(fam_pool_res[[i]]))))
  .debug_topk_diagnose(res_cv, label = sprintf("HALVING topk missing (post-fallback) | %s | stage %d/%d | cfg %d/%d | ", dist, s, S, i, nrow(fam_pool)))
}

stage_eval[[i]] <- list(
          score  = .stage_score(res_cv),
          cfg    = cfg,
          res    = res_cv,
          topk_w = topk_payload
        )
        progress_step(sprintf("HALVING %s | stage %d/%d | cfg %d/%d",
                              dist, s, S, i, nrow(fam_pool)))
      }
      
      
keep_n <- if (s < S) max(1L, ceiling(nrow(fam_pool) / eta)) else finalists_per_family
keep_n <- min(keep_n, nrow(fam_pool))

# --- Selection rule ---
# We always rank by score, but at the FINAL halving stage (s==S) we prefer configs that
# actually produced a non-empty Top-K payload, otherwise the "FINAL EVAL (NO GA)" cannot run.
if (s == S) {
  ok_topk <- vapply(stage_eval, function(z) !is.null(z$topk_w) && length(z$topk_w) > 0L, logical(1))

  if (any(ok_topk)) {
    ord_all <- order(vapply(stage_eval, `[[`, numeric(1), "score"))
    ord_ok  <- ord_all[ord_all %in% which(ok_topk)]
    keep_ix <- ord_ok[seq_len(min(keep_n, length(ord_ok)))]

    # Backfill (diagnostic only) if we still need more configs than valid-topk ones
    if (length(keep_ix) < keep_n) {
      remaining <- setdiff(ord_all, keep_ix)
      keep_ix <- c(keep_ix, remaining[seq_len(keep_n - length(keep_ix))])
      .dbg(sprintf("HALVING selection backfill | %s | stage %d/%d | keep_n=%d | ok_topk=%d | backfilled=%d",
                   dist, s, S, keep_n, sum(ok_topk), keep_n - length(ord_ok)))
    }
  } else {
    keep_ix <- order(vapply(stage_eval, `[[`, numeric(1), "score"))[seq_len(keep_n)]
    .dbg(sprintf("HALVING selection WARNING: no candidates with valid topk at final stage | %s | stage %d/%d",
                 dist, s, S))
  }
} else {
  keep_ix <- order(vapply(stage_eval, `[[`, numeric(1), "score"))[seq_len(keep_n)]
}

      
      # >>> FIX #1: ACTUALIZAR EL POOL PARA EL SIGUIENTE STAGE
      fam_pool <- fam_pool[keep_ix, , drop = FALSE]
      
      # >>> FIX #1b: GUARDAR SOLO LOS RESULTADOS QUE SOBREVIVEN
      fam_pool_res  <- lapply(stage_eval[keep_ix], `[[`, "res")
      fam_pool_topk <- lapply(stage_eval[keep_ix], `[[`, "topk_w")
      message(sprintf("[HALVING] %s | stage %d/%d | kept=%d",
                      dist, s, S, nrow(fam_pool)))
    }
    
    # ===================== FINAL: EVALUAR TOP5 (NO GA) SOBRE 100% FULL GRID ===================
    if (is.null(fam_pool_res)  || length(fam_pool_res)  != nrow(fam_pool) ||
        is.null(fam_pool_topk) || length(fam_pool_topk) != nrow(fam_pool)) {
      stop("HALVING produjo desalineación: pool_res/topk no coincide con fam_pool. Revisa keep/update.")
    }


    .dbg_main(sprintf("ENTER FINAL EVAL | %s | finalists=%d | res_len=%d | topk_len=%d | has_cfg_tag_col=%s",
                      dist, nrow(fam_pool),
                      if (!is.null(fam_pool_res)) length(fam_pool_res) else -1L,
                      if (!is.null(fam_pool_topk)) length(fam_pool_topk) else -1L,
                      ("cfg_tag" %in% names(fam_pool))))
    if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
      .dbg_main(sprintf("FINAL fam_pool colnames: %s", paste(names(fam_pool), collapse=", ")))
      .dbg_main(sprintf("FINAL fam_pool first row summary: %s",
                        paste(sprintf("%s=%s", names(fam_pool)[1:min(8, ncol(fam_pool))],
                                      as.character(fam_pool[1, 1:min(8, ncol(fam_pool))])),
                              collapse=", ")))
    }

    for (i in seq_len(nrow(fam_pool))) {
      cfg <- fam_pool[i, ]


    # ---------- TRACEABILITY: ensure we always have a non-NA config tag ----------
    # Downstream writers/loggers may rely on a tag string to name files and rows.
    # If the pool does not carry cfg_tag/config_tag, rebuild one deterministically.
    tag_val <- NA_character_
    if ("cfg_tag" %in% names(cfg)) tag_val <- as.character(cfg$cfg_tag)
    if (!isTRUE(nzchar(tag_val)) && "config_tag" %in% names(cfg)) tag_val <- as.character(cfg$config_tag)
    if (!isTRUE(nzchar(tag_val)) || is.na(tag_val)) tag_val <- .make_cfg_tag_fallback(cfg, prefix = "FINAL")
    # ---------------------------------------------------------------------------
      if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
        has_tag <- ("cfg_tag" %in% names(cfg))
        .dbg_main(sprintf("FINAL loop | %s | i=%d/%d | seed=%s | has_cfg_tag=%s | cfg_tag=%s",
                          dist, i, nrow(fam_pool), as.character(cfg$seed), has_tag, tag_val))
      }
      top5_w <- fam_pool_topk[[i]]

      
if (is.null(top5_w) || length(top5_w) == 0) {
  .dbg(sprintf("FINAL EVAL missing topk | %s | finalist i=%d/%d | tag=%s | fam_pool_topk_len=%d",
               dist, i, nrow(fam_pool),
               if (!is.null(cfg$tag)) cfg$tag else "<no-tag>",
               if (!is.null(fam_pool_topk)) length(fam_pool_topk) else -1L))
  if (!is.null(fam_pool_res) && length(fam_pool_res) >= i) {
    .debug_topk_diagnose(fam_pool_res[[i]], label = sprintf("FINAL EVAL diag from fam_pool_res | %s | i=%d | ", dist, i))
  }
  stop("HALVING: no se pudieron extraer pesos (top-k) para el final. Revisa el stage que produjo este resultado.")
}

      # >>> FIX #3: normaliza top5_w (matriz) y limita a 5
      top5_w <- as.matrix(top5_w)
      if (nrow(top5_w) > 5L) top5_w <- top5_w[1:5, , drop = FALSE]
      
      cat(sprintf("\n[%s] FINAL EVAL (NO GA): evaluando Top-%d sobre FULL grid...\n",
                  dist, nrow(top5_w)))
      
      eval_topk <- .evaluate_topk_fixed_weights(
        dist_name       = dist,
        dist_param_grid = grid,
        weights_mat     = top5_w,
        sample_sizes    = sample_sizes,
        num_samples     = num_samples,
        scenario_mode   = "full",
        seed            = cfg$seed + 50000L,
        q95_B           = bootstrap_B,
        crn_env         = crn_env,
        fam_key         = dist,
        rank_metric     = "robust_q95_mse",
        do_perturb      = TRUE
      )
      
      # Seleccionar ganador Top-1
      rank_df <- eval_topk$rank |>
        dplyr::mutate(
          rank_pos = dplyr::row_number(),
          is_top5  = rank_pos <= 5L,
          is_top1  = rank_pos == 1L
        )
      
      winner_row <- dplyr::filter(rank_df, is_top1) |> dplyr::slice(1)
      if (nrow(winner_row) == 0L) {
        if (identical(Sys.getenv("DEBUG_TOPK","0"), "1")) {
          message(sprintf("[DEBUG %s] WINNER row missing | %s | seed=%d | cfg_tag=%s | rank_cols=%s", 
                          format(Sys.time(), "%Y-%m-%d %H:%M:%S"), dist, cfg$seed, tag_val, paste(names(rank_df), collapse=",")))
        }
        stop("FINAL EVAL: winner_row vacío; eval_topk$rank no produjo top-1.")
      }

      winner_k   <- winner_row$candidate_k[1]
      winner_str <- winner_row$estimator_str[1]
      

      # --- DEBUG/EXPORT: store winner summary for this family/config ---


      # Robust cfg tag: keep previously computed deterministic tag_val (with fallback)
      tag_val <- as.character(tag_val)[1]

if (identical(Sys.getenv("DEBUG_MAIN","0"), "1")) {
  .dbg_main(sprintf("ABOUT TO BUILD win_one | %s | seed=%s | tag_val=%s",
                    dist, as.character(cfg$seed),
                    ifelse(is.na(tag_val), "<NA>", tag_val)))
}
      win_one <- tryCatch(data.frame(
        family = dist,
        seed   = cfg$seed,
        cfg_tag = tag_val,
        winner_k = winner_k,
        winner_estimator_str = winner_str,
        winner_robust_q95_mse  = as.numeric(winner_row$robust_q95_mse[1]),
        winner_robust_mean_mse = if ("robust_mean_mse" %in% names(winner_row)) as.numeric(winner_row$robust_mean_mse[1]) else NA_real_,
        winner_robust_max_mse  = if ("robust_max_mse" %in% names(winner_row))  as.numeric(winner_row$robust_max_mse[1])  else NA_real_,
        stringsAsFactors = FALSE
      ), error = function(e) {
        .dbg_main(sprintf("win_one ERROR | %s | message=%s", dist, conditionMessage(e)))
        .dbg_main(sprintf("win_one ERROR diagnostics | cfg_names=%s", paste(names(cfg), collapse=",")))
        stop(e)
      })
      winners_acc[[length(winners_acc) + 1L]] <- win_one
      if (identical(Sys.getenv("DEBUG_TOPK","0"), "1")) {
        message(sprintf("[DEBUG %s] WINNER stored | %s | seed=%d | cfg_tag=%s | k=%d",
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                        dist, cfg$seed, tag_val, winner_k))
      }

      # sensitivity analysis on final winner weights ---
      w_best <- as.numeric(eval_topk$candidates[[winner_k]]$weights)
      run_sensitivity(
        w_best = w_best,
        dist_name = dist,
        dist_param_grid = grid,
        run_dir = fam_dir,
        scenario_mode = "full",
        sample_sizes = sample_sizes,
        num_samples = num_samples,
        base_seed = 9000,
        n_alt = 6
      )
      
      index_add(dist, "SENS", file.path(fam_dir, sprintf("sensitivity_%s.csv", dist)), n_candidates = NA_integer_, n_scenarios = NA_integer_, seed = cfg$seed)
      
      cat(sprintf("[%s] WINNER (Top-1): k=%d | %s | robust_q95_mse=%.6g\n",
                  dist, winner_k, winner_str, winner_row$robust_q95_mse[1]))
      
      # Adjuntar info del ganador a cada fila del ranking
      rank_df <- rank_df |>
        dplyr::mutate(
          winner_k = winner_k,
          winner_estimator_str = winner_str
        )
      
      # ------------------ Guardar outputs ------------------
      if (!isTRUE(suppress_intermediate_saves)) {
        
        # >>> FIX #2: usar root_out (incluye run_id)
        fam_dir <- file.path(root_out, dist)
        mkdirp(fam_dir)
        
        # Ranking Top-k + flags + winner info
        readr::write_csv(rank_df, file.path(fam_dir, sprintf("FINAL_TOP5_RANK__seed%d.csv", cfg$seed)))
        
        # Escenarios del ganador
        best_tbl <- eval_topk$candidates[[winner_k]]$scenario_table
        best_tbl <- best_tbl |>
          dplyr::mutate(
            winner_k = winner_k,
            winner_estimator_str = winner_str,
            winner_robust_q95_mse  = winner_row$robust_q95_mse[1],
            winner_robust_mean_mse = winner_row$robust_mean_mse[1],
            winner_robust_max_mse  = winner_row$robust_max_mse[1]
          )
        
        readr::write_csv(best_tbl, file.path(fam_dir, sprintf("FINAL_WINNER_SCENARIOS__seed%d.csv", cfg$seed)))
        
        # Resumen 1-fila del ganador
        winner_summary <- winner_row |>
          dplyr::select(
            winner_k,
            winner_estimator_str = estimator_str,
            winner_robust_mean_mse = robust_mean_mse,
            winner_robust_q95_mse  = robust_q95_mse,
            winner_robust_max_mse  = robust_max_mse,
            pert_mean, pert_sd, pert_q95, pert_min, pert_max
          )
        
        readr::write_csv(winner_summary, file.path(fam_dir, sprintf("FINAL_WINNER_SUMMARY__seed%d.csv", cfg$seed)))
      }
    }
    
  } # <-- fin loop familias
  
  
  # ======================= CIERRE (SIN ROLLUP / SIN MEGA-CSV) =======================
  
  write_manifest(root_out, meta)
  index_flush(run_dir)
  progress_finalize()
  
  # ======================= WINNERS (ROLLUP) =======================
  #this builds the winners list
  winners_df <- if (length(winners_acc) > 0L) { 
    dplyr::bind_rows(winners_acc)
  } else {
    NULL
  }
  
  #Creates a unique CSV with all winners
  if (!is.null(winners_df) && nrow(winners_df) > 0L) {
    safe_write_csv(winners_df, file.path(root_out, "WINNERS__ALL_FAMILIES.csv"))
  }
  invisible(list(out_dir = root_out, winners = winners_df))
}


# ====================== MULTI-NODE WRAPPER ======================
# Runs the same pipeline, but distributes families across multiple machines.
# Design:
#   - One PSOCK worker per host (master+workers)
#   - Each worker receives a subset of families (capacity-aware)
#   - Inside each node, existing use_parallel=TRUE uses local cores

run_all_one_shot_multinode <- function(
    hosts,
    families_to_run,
    script_path = NULL,
    bench_seconds = 0.25,
    seed = 123,
    outfile = "",
    family_weights = NULL,
    ...
) {
  stopifnot(is.character(hosts), length(hosts) >= 1L)
  stopifnot(is.character(families_to_run), length(families_to_run) >= 1L)
  
  # 1) Build controller (1 worker per node) and benchmark capacity
  ctrl <- setup_multinode_controller(hosts, bench_seconds = bench_seconds, seed = seed, outfile = outfile)
  cl   <- ctrl$cl
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  # 2) Ensure workers have the full function environment
  if (is.null(script_path)) {
    script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
    if (!is.character(script_path) || is.na(script_path) || !nzchar(script_path)) script_path <- NA_character_
  }
  if (!is.character(script_path) || is.na(script_path) || !nzchar(script_path)) {
    stop("run_all_one_shot_multinode requires script_path that exists on ALL nodes (shared filesystem recommended).")
  }
  # Source the script on each worker node.
  parallel::clusterCall(cl, function(p) { source(p, local = FALSE); NULL }, script_path)
  
  # 3) Capacity-aware family assignment
  fams <- families_to_run
  
  if (is.null(family_weights)) {
    if (exists("param_grids", mode = "list")) {
      family_weights <- vapply(fams, function(d) {
        g <- param_grids[[d]]
        if (is.null(g)) 1 else nrow(g)
      }, numeric(1))
      names(family_weights) <- fams
    } else {
      family_weights <- setNames(rep(1, length(fams)), fams)
    }
  }
  
  buckets <- .assign_families_by_capacity(fams, node_scores = ctrl$scores, family_weights = family_weights)
  
  # 4) Run each bucket on its node. Each node will create its own out_dir.
  # Make output paths node-specific to avoid collisions.
  res <- parallel::parLapply(cl, seq_along(buckets), function(i) {
    my_fams <- buckets[[i]]
    if (length(my_fams) == 0L) return(NULL)
    run_all_one_shot(families_to_run = my_fams, ...)
  })
  
  list(
    hosts = hosts,
    node_caps = ctrl$caps,
    node_scores = ctrl$scores,
    family_buckets = buckets,
    results = res
  )
}


# # =============================== FINAL RUN ===============================
# # Objetivo: ejecutar la corrida completa sobre TODAS las familias principales.
# # Resultados se guardan en OUT_ROOT/GA_UNIFIED_<timestamp>
# #
# # REGLA L2: etapa final SIN GA (solo evaluar Top-K fijo en FULL grid).
# # => Por eso: final_retrain = FALSE
# # => Y empujamos K finalistas a full-budget: finalists_per_family = 5
# 
# if (identical(Sys.getenv("RUN_FINAL", "0"), "1")) {
# out_final <- run_all_one_shot(
#   families_to_run = c(
#     "normal",
#     "lognormal",
#     "weibull",
#     "exgaussian",
#     "exwald"
#   ),
#   
#   pop_sizes       = c(80, 100),
#   mutation_rates  = c(0.12, 0.18),
#   init_alphas     = c(0.5, 1.0),
#   alpha_mut_set   = c(0.5, 1.0),
#   immigrant_rates = c(0.05, 0.10),
#   t_sizes         = c(2L, 3L),
#   elitism_set     = c(1L, 2L),
#   seeds           = c(101, 202),
#   
#   # escenarios / CV
#   sample_sizes    = c(300, 500, 1000, 2000, 5000),
#   num_samples     = 100,
#   generations_per_fold = 60,
#   k_folds         = 3,
#   
#   # Final = evaluar pesos finalistas en FULL, SIN GA/retrain.
#   final_retrain   = FALSE,
#   
#   # objetivo y penalizaciones
#   objective       = "mixed",       # mezcla q95 + max (ver mix_w_* abajo)
#   mix_w_q95       = 0.7,
#   mix_w_max       = 0.3,
#   lambda_instab_default = 0.15,
#   bootstrap_B     = 300L,
#   check_every     = 5L,
#   patience        = 3L,
#   min_delta       = 0.005,
#   
#   use_parallel    = TRUE,
#   out_root        = OUT_ROOT,
#   
#   pick_single_winner_per_family = TRUE,
#   winner_metric   = "robust_q95_mse",  # o "robust_mean_mse" / "robust_max_mse"
#   suppress_intermediate_saves   = FALSE, # FALSE = CSV por cada config/stage
#   
#   # ===== FASES ESTRICTAS DEL GRIND DE HP =====
#   hp_strict = TRUE,
#   
#   hp1_frac = 0.25,      # ~25% del pool (o usa hp1_k)
#   hp1_k    = NULL,
#   hp1_gens = 10L,
#   hp1_num_samples = NULL,      # por defecto 40% de num_samples
#   hp1_bootstrap_B = NULL,      # por defecto ~33% de bootstrap_B
#   hp1_minibatch_frac = 0.50,   # minibatch de escenarios (si modo=light)
#   hp1_minibatch_min  = 10L,
#   
#   hp2_k    = 20L,
#   hp2_gens = 20L,
#   hp2_num_samples = NULL,      # por defecto 60% de num_samples
#   hp2_bootstrap_B = NULL,      # por defecto 60% de bootstrap_B
#   hp2_minibatch_frac = 0.50,
#   hp2_minibatch_min  = 12L,
#   
#   # Pool que pasa a Fase 3 (Halving)
#   hp3_pool_k = 10L,
#   
#   # ===== Fase 3 — HALVING (multi-fidelity, sin minibatch de escenarios) =====
#   search_strategy      = "halving",
#   eta                  = 3,                 # reduce a ~1/3 en cada etapa
#   
#   finalists_per_family = 5,                 # Top-5 (cumple L2)
#   
#   stage_gens_frac      = c(0.25, 0.60, 1.00),
#   stage_pop_size       = c(40,   60,   NA), # NA -> usa max(pop_sizes)
#   stage_num_samples    = c(30,   60,   NA), # NA -> usa num_samples completo
#   stage_bootstrap_B    = c(60L,  100L, NA),
#   
#   # (lo siguiente solo aplica si el modo light se activara en Fases 1/2)
#   minibatch_frac = 0.40,
#   minibatch_min  = 12L)
# }






# =============================== MAIN (AUTO-RUN) ===============================
# By default:
# - it runs a QUICK smoke run so you always get outputs.
# - In non-interactiv: does NOTHING unless set RUN_MODE.

RUN_MODE <- Sys.getenv("RUN_MODE", unset = if (interactive()) "quick" else "none")
RUN_MODE <- tolower(trimws(RUN_MODE))

if (RUN_MODE %in% c("quick", "final")) {
  catf("RUN_MODE = %s", RUN_MODE)
  catf("OUT_ROOT  = %s", OUT_ROOT)
  
  if (RUN_MODE == "quick") {
    
    # Ultra-fast smoke-real run -> should produce CSVs/logs
    families_quick <- c("normal") #For now takes only normal family
    
    out_quick <- run_all_one_shot(
      families_to_run = families_quick,
      
      # ===== GA hyperparams: micro pool (1 config) =====
      pop_sizes       = c(20),
      mutation_rates  = c(0.15),
      init_alphas     = c(0.7),
      alpha_mut_set   = c(0.7),
      immigrant_rates = c(0.08),
      t_sizes         = c(2L),
      elitism_set     = c(1L),
      seeds           = c(101),
      
      sample_sizes         = c(300),  # 1 tamaño para speed
      num_samples          = 10,      # muy bajo pero imprime resultados
      generations_per_fold = 8,       # pocas generaciones
      k_folds              = 2,       # CV real mínimo
      
      # ===== Regla L2 =====
      final_retrain         = FALSE,
      finalists_per_family  = 2,
      
      # ===== objetivo y penalizaciones =====
      objective             = "mixed",
      mix_w_q95             = 0.7,
      mix_w_max             = 0.3,
      lambda_instab_default = 0.15,
      bootstrap_B           = 30L,    # bajo
      check_every           = 2L,
      patience              = 1L,
      min_delta             = 0.02,
      
      use_parallel                 = TRUE,
      out_root                     = OUT_ROOT,
      suppress_intermediate_saves  = FALSE,  
      pick_single_winner_per_family = TRUE,
      winner_metric                = "robust_q95_mse",
      
      # ===== HP strict: activo, pero micro =====
      hp_strict = TRUE,
      
      # HPF1: corto
      hp1_frac            = 0.50,
      hp1_k               = NULL,
      hp1_gens            = 3L,
      hp1_num_samples     = 6,
      hp1_bootstrap_B     = 20L,
      hp1_minibatch_frac  = 0.60,
      hp1_minibatch_min   = 6L,
      
      # HPF2: corto
      hp2_k               = 4L,
      hp2_gens            = 4L,
      hp2_num_samples     = 8,
      hp2_bootstrap_B     = 25L,
      hp2_minibatch_frac  = 0.60,
      hp2_minibatch_min   = 8L,
      
      # Pool -> halving: mini
      hp3_pool_k = 3L,
      
      # ===== HALVING: 2 stages (todo alineado) =====
      search_strategy     = "halving",
      eta                 = 2,
      
      stage_gens_frac      = c(0.50, 1.00),
      stage_scenario_frac  = c(0.25, 1.00),  
      
      stage_pop_size       = c(20, NA),
      stage_num_samples    = c(8,  NA),
      stage_bootstrap_B    = c(25L, NA),
      
    
      minibatch_frac = 0.50,
      minibatch_min  = 8L
    )
    
    print(out_quick)
    catf("Quick run complete. Check: %s",
         normalizePath(file.path(OUT_ROOT), winslash = "/", mustWork = FALSE))
    
  } else if (RUN_MODE == "final") {
    
    # Full run: heavier budgets. 
    out_final <- run_all_one_shot(
      families_to_run = c("normal","lognormal","weibull","exgaussian","exwald"),
      pop_sizes       = c(80, 100),
      mutation_rates  = c(0.12, 0.18),
      init_alphas     = c(0.5, 1.0),
      alpha_mut_set   = c(0.5, 1.0),
      immigrant_rates = c(0.05, 0.10),
      t_sizes         = c(2L, 3L),
      elitism_set     = c(1L, 2L),
      seeds           = c(101, 202),
      
      sample_sizes    = c(300, 500, 1000, 2000, 5000),
      num_samples     = 100,
      generations_per_fold = 60,
      k_folds         = 3,
      
      # REGLA L2:
      final_retrain   = FALSE,
      finalists_per_family = 5,
      
      objective       = "mixed",
      mix_w_q95       = 0.7,
      mix_w_max       = 0.3,
      lambda_instab_default = 0.15,
      bootstrap_B     = 300L,
      check_every     = 5L,
      patience        = 3L,
      min_delta       = 0.005,
      
      use_parallel    = TRUE,
      out_root        = OUT_ROOT,
      suppress_intermediate_saves = FALSE,
      
      pick_single_winner_per_family = TRUE,
      winner_metric   = "robust_q95_mse"
    )
    
    catf("Final run complete. Check: %s",
         normalizePath(file.path(OUT_ROOT), winslash = "/", mustWork = FALSE))
  }
  
} else {
  invisible(NULL)
}