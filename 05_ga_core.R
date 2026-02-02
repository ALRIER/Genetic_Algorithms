
# =============================================================================
# 05_GA_CORE — RESEARCH COMMENTED EDITION
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




init_population <- function(size, n_estimators, alpha = 1, warm_start = NULL, jitter = 0.15) {
  stopifnot(size >= 1L, n_estimators >= 1L, is.finite(alpha), alpha > 0)
  # Base Dirichlet
  pop <- matrix(gtools::rdirichlet(size, rep(alpha, n_estimators)),
                nrow = size, byrow = TRUE)
  
  if (!is.null(warm_start)) {
    ws <- as.matrix(warm_start)
    if (ncol(ws) != n_estimators) {
      if (ncol(ws) > n_estimators) ws <- ws[, seq_len(n_estimators), drop = FALSE]
      if (ncol(ws) < n_estimators) {
        add <- matrix(0, nrow = nrow(ws), ncol = n_estimators - ncol(ws))
        ws <- cbind(ws, add)
      }
    }
    ws <- t(apply(ws, 1, .normalize_simplex))
    k <- min(nrow(ws), size)
    if (k > 0) {
      # jitter suave for diversity
      eps <- matrix(gtools::rdirichlet(k, rep(alpha + jitter, n_estimators)),
                    nrow = k, byrow = TRUE)
      pop[seq_len(k), ] <- t(apply((ws[seq_len(k), , drop = FALSE] + eps) / 2, 1, .normalize_simplex))
    }
  }
  # Normalizing
  pop <- t(apply(pop, 1, .normalize_simplex))
  pop
}

crossover <- function(parent1, parent2) {
  stopifnot(length(parent1) == length(parent2))
  a <- stats::runif(1)  # mezcla aleatoria
  child <- a * as.numeric(parent1) + (1 - a) * as.numeric(parent2)
  .normalize_simplex(child)
}

mutate_weights <- function(weights, mutation_rate = 0.10, alpha_mut = 1) {
  w <- .normalize_simplex(as.numeric(weights))
  mutation_rate <- min(max(mutation_rate, 0), 1)  # clip a [0,1]
  if (stats::runif(1) < mutation_rate) {
    alpha_mut <- max(alpha_mut, 1e-6)
    perturb <- as.numeric(gtools::rdirichlet(1, rep(alpha_mut, length(w))))
    w <- (w + perturb) / 2
  }
  .normalize_simplex(w)
}




tournament_select <- function(scores, pop_size, t_size = 2L, replace = TRUE) {
  stopifnot(is.numeric(scores), length(scores) >= 2L, pop_size >= 2L, t_size >= 2L)
  n <- length(scores)
  t_size  <- min(t_size, n)             # security if t_size > n
  n_pairs <- max(1L, pop_size %/% 2L) #size of mating pool
  winners <- integer(n_pairs)
  
  for (i in seq_len(n_pairs)) {
    cand <- base::sample.int(n, size = t_size, replace = FALSE)
    local_best <- cand[which.min(scores[cand])]
    winners[i] <- local_best
  }
  
  if (!replace) winners <- unique(winners)
  if (!replace && length(winners) < n_pairs) {
    ord  <- base::order(scores) # complete with bests... scores
    need <- n_pairs - length(winners)
    winners <- c(winners, setdiff(ord, winners)[seq_len(need)])
  }
  winners
}

make_folds <- function(n, k, seed = NULL) {
  stopifnot(k >= 2L, n >= k)
  if (!is.null(seed)) {
    s_bkp <- if (exists(".Random.seed", inherits = FALSE)) .Random.seed else NULL
    on.exit({
      if (!is.null(s_bkp)) {
        .Random.seed <<- s_bkp
      } else {
        # Avoid noisy warning in fresh sessions / PSOCK workers
        if (exists(".Random.seed", inherits = FALSE)) rm(".Random.seed", inherits = FALSE)
      }
    }, add = TRUE)
    seed <- .ensure_seed(seed, fallback = 101L)
    set.seed(seed)
  }
  ids <- base::sample(seq_len(n))
  split(ids, rep(seq_len(k), length.out = n))
}

build_warm_starts <- function() {
  I   <- diag(N_EST)
  avg <- matrix(rep(1 / N_EST, N_EST), nrow = 1)
  
  z <- function(names) {
    v <- numeric(N_EST)
    present <- ESTIMATOR_NAMES %in% names
    k <- sum(present)
    if (k > 0) v[present] <- 1 / k
    v
  }
  
  combs <- rbind(
    z(c("mean","median")),
    z(c("median","trimmed20")),
    z(c("median","trimean","biweight")),
    z(c("huber","biweight")),
    z(c("geometric","harmonic"))
  )
  
  W <- unique(rbind(I, avg, combs))
  W <- t(apply(W, 1, .normalize_simplex))
  rs <- rowSums(W)
  if (any(rs == 0)) {
    W[rs == 0, ] <- matrix(rep(1 / N_EST, N_EST), nrow = sum(rs == 0), byrow = TRUE)
  }
  unique(W)
}


summarize_scenario <- function(S, weights, distribution, sample_size,
                               q95_B = 300L, q95_seed = NULL, q95_boot_indices = NULL) {
  true_mu <- S$true_mean
  C <- S$components_by_size[[as.character(sample_size)]]
  
  if (is.null(C) || !is.matrix(C) || ncol(C) != N_EST) {
    base_row <- tibble::tibble(
      distribution = distribution,
      sample_size = as.integer(sample_size),
      contamination_rate = S$scenario$contamination_rate,
      outlier_scale_mad  = S$scenario$outlier_scale_mad,
      contamination_type = S$scenario$contamination_type,
      true_mean = true_mu,
      scenario_mode = { sm <- attr(S, "scenario_mode", exact = TRUE); if (is.null(sm)) "unknown" else sm }
    )
    
    na_metrics <- list(
      mse=NA_real_, mse_q95=NA_real_, mse_max=NA_real_,
      bias=NA_real_, abs_bias=NA_real_,
      variance=NA_real_, mad=NA_real_, iqr=NA_real_
    )
    
    make_row <- function(est_name) {
      dplyr::mutate(base_row, estimator = est_name, !!!na_metrics)
    }
    
    return(dplyr::bind_rows(
      make_row("robust"),
      make_row("mean"),
      make_row("median"),
      make_row("trimmed20"),
      make_row("harmonic"),
      make_row("geometric"),
      make_row("mode_hsm"),
      make_row("mode_parzen"),
      make_row("trimean"),
      make_row("huber"),
      make_row("biweight")
    ))
  }
  
  # robust column selector: if it does not exist it returns NA 
  #NEed to resolve... the NA could not be useful - contaminations. 
  col_safe <- function(nm) {
    j <- match(nm, ESTIMATOR_NAMES)
    if (is.na(j)) rep(NA_real_, nrow(C)) else C[, j, drop = TRUE]
  }
  
  est_mean    <- col_safe("mean")
  est_median  <- col_safe("median")
  est_trim20  <- col_safe("trimmed20")
  est_harm    <- col_safe("harmonic")
  est_geom    <- col_safe("geometric")
  est_mode_h  <- col_safe("mode_hsm")
  est_mode_p  <- col_safe("mode_parzen")
  est_trimean <- col_safe("trimean")
  est_huber   <- col_safe("huber")
  est_biwt    <- col_safe("biweight")
  
  w_robust   <- .normalize_simplex(as.numeric(weights))
  est_robust <- as.vector(C %*% w_robust)
  
  agg <- function(z) {
    z <- as.numeric(z); z <- z[is.finite(z)]
    if (length(z) == 0L) {
      return(c(mse=NA, mse_q95=NA, mse_max=NA,
               bias=NA, abs_bias=NA,
               variance=NA, mad=NA, iqr=NA))
    }
    errs2 <- (z - true_mu)^2
    c(
      mse      = base::mean(errs2, na.rm = TRUE),
      mse_q95  = q95_boot(errs2, B = q95_B, rng_seed = q95_seed, boot_indices = q95_boot_indices),
      mse_max  = base::max(errs2, na.rm = TRUE),
      bias     = base::mean(z, na.rm = TRUE) - true_mu,
      abs_bias = base::abs(base::mean(z, na.rm = TRUE) - true_mu),
      variance = stats::var(z, na.rm = TRUE),
      mad      = stats::mad(z, constant = 1, na.rm = TRUE),
      iqr      = stats::IQR(z, na.rm = TRUE)
    )
  }
  
  m_mean    <- agg(est_mean)
  m_median  <- agg(est_median)
  m_trim20  <- agg(est_trim20)
  m_harm    <- agg(est_harm)
  m_geom    <- agg(est_geom)
  m_mode_h  <- agg(est_mode_h)
  m_mode_p  <- agg(est_mode_p)
  m_trimean <- agg(est_trimean)
  m_huber   <- agg(est_huber)
  m_biwt    <- agg(est_biwt)
  m_robust  <- agg(est_robust)
  
  base_row <- tibble::tibble(
    distribution = distribution,
    sample_size = as.integer(sample_size),
    contamination_rate = S$scenario$contamination_rate,
    outlier_scale_mad  = S$scenario$outlier_scale_mad,
    contamination_type = S$scenario$contamination_type,
    true_mean = true_mu,
    scenario_mode = { sm <- attr(S, "scenario_mode", exact = TRUE); if (is.null(sm)) "unknown" else sm }
  )
  
  make_row2 <- function(est_name, m) {
    dplyr::mutate(base_row, estimator = est_name, !!!as.list(m))
  }
  
  dplyr::bind_rows(
    make_row2("robust",      m_robust),
    make_row2("mean",        m_mean),
    make_row2("median",      m_median),
    make_row2("trimmed20",   m_trim20),
    make_row2("harmonic",    m_harm),
    make_row2("geometric",   m_geom),
    make_row2("mode_hsm",    m_mode_h),
    make_row2("mode_parzen", m_mode_p),
    make_row2("trimean",     m_trimean),
    make_row2("huber",       m_huber),
    make_row2("biweight",    m_biwt)
  )
}

.ga_default_ctrl <- list(
  objective           = "q95",
  use_hierarchical = TRUE,
  group_by         = c("type","rate_bin","scale_bin"),
  group_inner      = "q95",
  group_inner_q    = 0.95,
  lambda_instab       = 0.15,
  bootstrap_B         = 200L,
  t_size              = 2L,
  elitism             = 2L,
  immigrant_rate      = 0.10,
  mutation_rate_init  = 0.18,
  init_alpha          = 1.0,
  alpha_mut           = 1.0,
  check_every         = 5L,
  patience            = 3L,
  min_delta           = 0.005,
  mix_w_q95           = 0.7,
  mix_w_max           = 0.3,
  lambda_entropy      = 0.015,
  max_w_dominance     = 0.80,
  lambda_dominance    = 0.05,
  lambda_bias         = 0.05,
  minibatch_frac      = NULL,#fraction of train scenarios -only in lightmode
  minibatch_min       = 12L #minibatch/ minimum number scenarios
)
#Controls scenarios picks up light/full 
.ga_pick_mode <- function(sample_sizes, num_samples, pop_size, generations) {
  n_scen_full <- nrow(build_scenarios_full())
  work_score  <- n_scen_full * length(sample_sizes) * num_samples * pop_size * generations
  if (work_score >= 5e7) "light" else "full"
}

# # Cluster PSOCK con RNG reproducible
#Creates a PSOCK cluster (separate R worker processes) 
#and assigns each worker an independent, reproducible 
#RNG (random number generator) stream using a seed. 
#Ensures parallel simulations and genetic algorithm runs 
#produce consistent yet non-overlapping random sequences 
#across cores.
.ga_setup_cluster <- function(use_parallel, seed) {
  if (!use_parallel) return(NULL)
  n_cores <- suppressWarnings(max(1L, parallel::detectCores() - 1L))
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::clusterEvalQ(cl, { runif(1); NULL })
  cl
}

# Exports required functions, data, and settings to PSOCK workers 
# so each parallel R process can run the GA fitness and diagnostics. 
# Then loads needed libraries on each worker to ensure all dependencies 
# are available during parallel evaluation.
.ga_export_cluster <- function(cl, env) {
  if (is.null(cl)) return(invisible(NULL))
  parallel::clusterExport(
    cl,
    varlist = c(
      # funciones de fitness y dependencias
      "fitness_universal","q95_boot","huber_loss","scenario_weight",".normalize_simplex",
      # ---hierarchical grouping helpers (needed on PSOCK workers) ---
      "scenario_group_key",".rate_bin",".scale_bin","hierarchical_objective",
      # ---diagnostics helpers (only used if return_details=TRUE) ---
      "summarize_losses","group_risk_stats",
      "prepped","N_EST","ESTIMATOR_NAMES",".ga_eval_split","ctrl"
    ),
    envir = env
  )
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(dplyr); library(tibble); library(MASS) })
    NULL
  })
}
# ====================== MULTI-NODE ORCHESTRATION ====================

#This is important sice the complete training in a normal computer can take
#up a 3 months, so I produce a clustered solution based on 3 computers
#I have at home>>>>>>>>
#Orchestrates 2-level parallelism safely: one PSOCK worker per remote node (Level A),
# while each node can still run its own local-core parallelism (Level B) via .ga_setup_cluster().
# This avoids oversubscribing cores and lets the master distribute work by distribution-family.
# Goal: 2-level parallelism without oversubscription
#   - Level A (multi-node): master assigns *families* to nodes
#       -> we create 1 PSOCK worker per node
#   - Level B (intra-node): inside each node-process we keep the existing
#       use_parallel=TRUE logic, which will fan out across local cores via
#       .ga_setup_cluster().


#Benchmarks each node’s compute capacity using a short, CPU-heavy loop (matrix ops).
# Returns cores detected, repetitions completed, elapsed time, and a capacity score used to load-balance
# work across nodes (faster nodes get more work; slower nodes get less).
.node_benchmark <- function(seconds = 0.25) {
  cores <- suppressWarnings(parallel::detectCores())
  if (!is.finite(cores) || is.na(cores)) cores <- 1L
  t0 <- proc.time()[[3]]
  reps <- 0L
  while ((proc.time()[[3]] - t0) < seconds) {
    x <- matrix(runif(2500), nrow = 250)
    y <- crossprod(x)
    reps <- reps + 1L
  }
  elapsed <- proc.time()[[3]] - t0
  score <- (reps / max(elapsed, 1e-6)) * as.numeric(cores)
  list(cores = as.integer(cores), reps = reps, elapsed = elapsed, score = as.numeric(score))
}

#Assigns distribution families to nodes based on node capacity scores (and optional
# family weights). Uses greedy bin-packing on “effective load” = weight/score so each node receives
# a balanced share of work, minimizing total wall-clock time across heterogeneous machines.
.assign_families_by_capacity <- function(families, node_scores, family_weights = NULL) {
  stopifnot(is.character(families), length(families) >= 1L)
  stopifnot(is.numeric(node_scores), length(node_scores) >= 1L)
  node_scores <- pmax(node_scores, 1e-6)
  
  if (is.null(family_weights)) {
    family_weights <- rep(1, length(families))
    names(family_weights) <- families
  } else {
    # ensure named weights
    if (is.null(names(family_weights))) {
      if (length(family_weights) != length(families)) {
        stop("family_weights must be NULL, or have same length as families, or be a named vector.")
      }
      names(family_weights) <- families
    }
    family_weights <- family_weights[families]
    family_weights[!is.finite(family_weights) | is.na(family_weights)] <- 1
    family_weights <- pmax(family_weights, 1e-6)
  }
  
  # Greedy bin-packing on "effective load" = weight / score
  ord <- order(family_weights, decreasing = TRUE)
  fams <- families[ord]
  wts  <- family_weights[fams]
  
  k <- length(node_scores)
  loads <- rep(0, k)
  buckets <- vector("list", k)
  for (i in seq_along(fams)) {
    j <- which.min(loads)
    buckets[[j]] <- c(buckets[[j]], fams[[i]])
    loads[[j]] <- loads[[j]] + (wts[[i]] / node_scores[[j]])
  }
  buckets
}

#Creates a PSOCK cluster across the provided hostnames (one worker per host),
# initializes reproducible RNG streams on workers using a seed, “warms up” RNG state, benchmarks
# each node with .node_benchmark(), and returns the cluster plus per-node capacity scores for scheduling.
setup_multinode_controller <- function(hosts,
                                       bench_seconds = 0.25,
                                       seed = 123,
                                       outfile = "") {
  stopifnot(is.character(hosts), length(hosts) >= 1L)
  cl_nodes <- parallel::makePSOCKcluster(hosts, outfile = outfile)
  parallel::clusterSetRNGStream(cl_nodes, iseed = seed)
  parallel::clusterEvalQ(cl_nodes, { runif(1); NULL })
  caps <- parallel::clusterEvalQ(cl_nodes, .node_benchmark(seconds = bench_seconds))
  scores <- vapply(caps, function(z) z$score, numeric(1))
  scores[!is.finite(scores)] <- 1
  list(cl = cl_nodes, caps = caps, scores = scores)
}


#Computes a smooth mutation-rate schedule that decays over generations: starts near
# m0 early to explore broadly, then gradually approaches mmin as g→gmax to stabilize convergence.
# Uses a log-based curve to avoid overly fast decay in early generations.
.ga_decay_mut <- function(g, gmax, m0, mmin = 0.05) {
  g <- max(1L, as.integer(g)); gmax <- max(1L, as.integer(gmax))
  mmin + (m0 - mmin) * (log(gmax + 1) - log(g + 1)) / log(gmax + 1)
}

#Wrapper for evaluating a candidate weight vector w on a subset of scenarios idx.
# It forwards GA control settings (penalties, mixing weights, bootstrap size, hierarchical grouping)
# into fitness_universal(), and sets q95_seed so bootstrap/q95 components remain reproducible.
.ga_eval_split <- function(w, prepped, idx, ctrl, q95_seed) {
  fitness_universal(
    w, prepped, objective = ctrl$objective,
    idx = idx, use_boot = TRUE, B = ctrl$bootstrap_B,
    lambda_instab    = ctrl$lambda_instab,
    lambda_entropy   = ctrl$lambda_entropy,
    max_w_dominance  = ctrl$max_w_dominance,
    lambda_dominance = ctrl$lambda_dominance,
    mix_w_q95        = ctrl$mix_w_q95,
    mix_w_max        = ctrl$mix_w_max,
    lambda_bias      = ctrl$lambda_bias,
    use_hierarchical = ctrl$use_hierarchical,
    group_by         = ctrl$group_by,
    group_inner      = ctrl$group_inner,
    group_inner_q    = ctrl$group_inner_q,
    q95_seed         = q95_seed
  )
}

# ===================== SENSITIVITY: WEIGHT PERTURBATION =====================
# Evaluates how stable an optimal weight vector is by adding small random
# perturbations in log-space (so weights remain positive and normalized). It re-runs the
# fitness on many nearby weight vectors and summarizes how much performance degrades,
# revealing local robustness vs fragility of the solution.

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
# Provides a non-evolutionary baseline by sampling many random weight
# vectors from the simplex and evaluating their fitness. This helps quantify how much
# the GA improves over chance. The best random solution (weights + diagnostics) is
# retained for fair comparison against GA-optimized results.

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
  
 
  n_estimators <- ncol(prepped[[1]]$components_by_size[[1]])# number of estimators
  
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


# ===================== SENSITIVITY: WEIGHT PERTURBATION =====================
#Generates nearby alternative weight vectors on the probability simplex
# by applying Gaussian noise in log-space and renormalizing. This preserves positivity
# and sum-to-one constraints while exploring the local neighborhood around w, enabling
# smooth robustness and sensitivity analyses of estimator combinations.

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

#It Re-runs the perturbation-based robustness test but also returns the
# full vector of fitness values. This allows deeper diagnostics such as distribution
# plots, tail-risk analysis, and identifying asymmetric sensitivity where some nearby
# weight directions cause much larger degradation than others.

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

# =========================== GA (single split) ===========================
#Main GA driver for one distribution family using a single train/validation split.
# Builds/uses a scenario set, initializes a simplex-constrained population 
# (optionally warm-started), evolves weights with selection/crossover/mutation, 
# evaluates on train/val, early-stops, and returns best weights plus detailed
# per-scenario and aggregate outputs.

evolve_universal_estimator_per_family <- function(dist_name,
                                                  dist_param_grid,
                                                  sample_sizes = c(300, 500, 1000, 2000, 5000),
                                                  num_samples = 80,
                                                  pop_size = 100,
                                                  generations = 50,
                                                  seed = 101,
                                                  objective = c("mean","q95","max","mixed"),
                                                  record_convergence = TRUE,
                                                  use_parallel = TRUE,
                                                  train_idx = NULL, val_idx = NULL,
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
                                                  warm_start_matrix = NULL,
                                                  return_val_best = TRUE,
                                                  mix_w_q95 = 0.7,
                                                  mix_w_max = 0.3,
                                                  lambda_entropy = 0.015,
                                                  max_w_dominance = 0.80,
                                                  lambda_dominance = 0.05,
                                                  lambda_bias = 0.05,
                                                  # minibatch only in LIGHT mode ----
                                                  minibatch_frac = NULL,
                                                  minibatch_min  = 12L,
                                                  # ---- CRN ----
                                                  crn_env = NULL,
                                                  fam_key = NULL,
                                                  force_scenario_mode = NULL,
                                                  #permite inyectar el prepped ya subseteado desde CV
                                                  prepped_override = NULL,
                                                  diagnostics = FALSE){
  
  #Validates key inputs and establishes a reproducible RNG state for the entire run.
  # This ensures deterministic behavior across GA initialization, scenario generation, folds/splits,
  # and any stochastic evaluation components that depend on the seed, including parallel RNG streams.
  stopifnot(pop_size >= 2L, generations >= 1L, length(sample_sizes) >= 1L)
  seed <- .ensure_seed(seed, fallback = 101L)
  set.seed(seed)
  
  #Consolidates all GA hyperparameters and objective settings into a single control
  # list (ctrl). This centralizes configuration used throughout selection/mutation/evaluation and
  # enables consistent logging/diagnostics and safe parameter passing to worker processes.
  ctrl <- modifyList(.ga_default_ctrl, list(
    objective          = match.arg(objective),
    lambda_instab      = lambda_instab,
    bootstrap_B        = as.integer(bootstrap_B),
    t_size             = as.integer(t_size),
    elitism            = as.integer(elitism),
    immigrant_rate     = immigrant_rate,
    mutation_rate_init = mutation_rate_init,
    init_alpha         = init_alpha,
    alpha_mut          = alpha_mut,
    check_every        = as.integer(check_every),
    patience           = as.integer(patience),
    min_delta          = min_delta,
    mix_w_q95          = mix_w_q95,
    mix_w_max          = mix_w_max,
    lambda_entropy     = lambda_entropy,
    max_w_dominance    = max_w_dominance,
    lambda_dominance   = lambda_dominance,
    lambda_bias        = lambda_bias,
    minibatch_frac     = minibatch_frac,
    minibatch_min      = as.integer(minibatch_min)
  ))
  ctrl$diagnostics <- isTRUE(diagnostics)
  ctrl$run_id  <- ctrl$run_id %||% NA_character_
  ctrl$ga_seed <- seed
  best_details_hist <- list()
  
  #Chooses the scenario universe mode (full vs light). FULL is heavier but more
  # comprehensive; LIGHT is faster and can use minibatching. If force_scenario_mode is provided,
  # it overrides auto-selection to ensure the caller controls computational cost and coverage.
  scenario_mode <- if (!is.null(force_scenario_mode)) {
    match.arg(force_scenario_mode, c("full","light"))
  } else {
    .ga_pick_mode(sample_sizes, num_samples, pop_size, generations)
  }
  
  #Builds or reuses the prepared scenario list (prepped) containing simulations and
  # estimator components for each scenario. In CV, prepped_override avoids regenerating scenarios
  # and guarantees train/val subsets refer to the same scenario universe and CRN configuration.
  if (!is.null(prepped_override)) {
    prepped <- prepped_override
  } else {
    prepped <- prep_scenarios(
      dist_name, dist_param_grid,
      sample_sizes  = sample_sizes,
      num_samples   = num_samples,
      scenario_mode = scenario_mode,
      seed          = seed,
      crn_env       = crn_env,
      fam_key       = if (is.null(fam_key)) dist_name else fam_key
    )
  }
  
  #Determines the effective scenario mode used by the prepped object (attribute wins)
  # so downstream logic (e.g., minibatching) matches how scenarios were actually generated. This is
  # defensive when a caller provides overrides or cached objects with their own scenario_mode.
  prepped_mode <- attr(prepped, "scenario_mode") %||% scenario_mode
  
  #Validates that train/val indices are within bounds of the scenario list. This
  # prevents silent recycling or out-of-range errors when the caller passes subsets (e.g., from CV)
  # and ensures the GA evaluates exactly the intended scenario subset for training and validation.
  if (!is.null(train_idx) && length(train_idx)) {
    if (max(train_idx) > length(prepped)) stop("train_idx excede length(prepped). Revisa el prepped_override/subset.")
  }
  if (!is.null(val_idx) && length(val_idx)) {
    if (max(val_idx) > length(prepped)) stop("val_idx excede length(prepped). Revisa el prepped_override/subset.")
  }
  
  
  #Computes a random-search baseline for context: how good can we do without GA,
  # just by sampling many random simplex weight vectors? This provides a fair comparison point to
  # quantify GA benefit and helps detect if the objective surface is easy (baseline already strong).
  # ---- RANDOM BASELINE (optional) ----
  rand_base <- random_search_baseline(
    prepped = prepped,
    N = 5000,
    objective = ctrl$objective,
    use_hierarchical = ctrl$use_hierarchical,
    group_by = ctrl$group_by,
    group_inner = ctrl$group_inner,
    group_inner_q = ctrl$group_inner_q,
    seed = seed,
    idx = val_idx
  )
  
  cat(sprintf("[BASELINE] %s | random_search | best=%.6g\n", dist_name, rand_base$fitness))
  
  #Establishes train/validation scenario splits. If indices are not provided, a
  # deterministic K=3 split is created and the last fold is used for validation. 
  # This keeps the run reproducible and ensures validation is held out
  # from training selection pressure.
  all_idx <- seq_along(prepped)
  if (is.null(val_idx) || is.null(train_idx)) {
    folds <- make_folds(length(all_idx), k = 3, seed = seed)
    val_idx   <- sort(unlist(folds[[3]]))
    train_idx <- sort(setdiff(all_idx, val_idx))
  } else {
    train_idx <- sort(unique(train_idx))
    val_idx   <- sort(unique(val_idx))
  }
  
  #Constructs the initial GA population on the simplex. Combines deterministic
  # warm-start candidates (built-in + cached + user-provided) with randomly generated individuals
  # to reach pop_size, then normalizes all rows to enforce positivity and sum-to-one constraints.
  # ---- Warm-starts ------------------------------------------------------
  warm_cache <- get_warm_start(dist_name, seed)
  warm <- build_warm_starts()
  if (!is.null(warm_start_matrix)) warm <- unique(rbind(warm_start_matrix, warm))
  if (!is.null(warm_cache))        warm <- unique(rbind(warm_cache, warm))
  rest <- max(0L, pop_size - nrow(warm))
  ga_population <- if (rest > 0L) rbind(warm, init_population(rest, N_EST, alpha = ctrl$init_alpha)) else warm
  ga_population <- t(apply(ga_population, 1, .normalize_simplex))
  
  #Sets up optional parallel evaluation using a PSOCK cluster with reproducible RNG
  # streams derived from seed. Exports required functions/data to workers and ensures the cluster is
  # properly stopped on exit, preventing orphaned processes and non-deterministic RNG behavior.
  cl <- .ga_setup_cluster(use_parallel, seed)
  on.exit({ if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  .ga_export_cluster(cl, environment())
  
  #Initializes convergence tracking and early-stopping state. Stores best/median
  # train and validation metrics over time, maintains a validation history, and tracks the best score
  # seen so far to trigger patience-based stopping and adaptive mutation boosts when stagnation occurs.
  # ---- Tracking ---------------------------------------------------------
  conv <- if (record_convergence)
    data.frame(gen=integer(0), best_train=numeric(0), med_train=numeric(0),
               best_val=numeric(0),   med_val=numeric(0)) else NULL
  val_hist <- numeric(0)
  best_val_so_far <- Inf
  no_improve <- 0L
  
  #Core evolutionary loop: per generation, compute a decayed mutation rate and a
  # generation-specific seed for q95/boot components, optionally minibatch train scenarios in LIGHT
  # mode, evaluate population fitness (parallel if available), then apply elitism/selection/crossover/
  # mutation and occasional immigrants to maintain diversity and avoid premature convergence.
  # ======================== Bucle evolutivo ==============================
  for (gen in seq_len(generations)) {
    mut_rate     <- .ga_decay_mut(gen, generations, ctrl$mutation_rate_init)
    q95_seed_gen <- 1e6 + seed*1000 + gen
    
    train_sub <- train_idx
    if (identical(prepped_mode, "light") && !is.null(ctrl$minibatch_frac)) {
      k <- ceiling(length(train_idx) * ctrl$minibatch_frac)
      k <- max(ctrl$minibatch_min, min(k, length(train_idx)))
      seed_mb <- 9e6 + seed*1000 + gen
      train_sub <- balanced_sample_idx(prepped, train_idx, k, seed_mb)
    }
    
    eval_fun <- function(w)
      .ga_eval_split(w, prepped, train_sub, ctrl, q95_seed = q95_seed_gen)
    
    fitness_scores <- if (!is.null(cl)) {
      unlist(parallel::parLapply(
        cl,
        X = seq_len(nrow(ga_population)),
        fun = function(i, pop) eval_fun(pop[i, ]),
        ga_population
      ))
    } else {
      apply(ga_population, 1, eval_fun)
    }
    
    ord_fit <- order(fitness_scores)
    n_elite <- min(ctrl$elitism, nrow(ga_population))
    elites  <- ga_population[ord_fit[seq_len(n_elite)], , drop = FALSE]
    
    sel_idx  <- tournament_select(fitness_scores, nrow(ga_population),
                                  t_size = ctrl$t_size, replace = TRUE)
    selected <- ga_population[sel_idx, , drop = FALSE]
    
    n_off <- nrow(ga_population) - nrow(elites)
    offspring <- if (n_off > 0L) t(replicate(n_off, {
      p1 <- selected[base::sample(nrow(selected), 1L), ]
      p2 <- selected[base::sample(nrow(selected), 1L), ]
      mutate_weights(crossover(p1, p2), mutation_rate = mut_rate, alpha_mut = ctrl$alpha_mut)
    })) else matrix(numeric(0), 0, ncol(ga_population))
    
    ga_population <- rbind(elites, offspring)
    
    if (gen %% 15L == 0L && ctrl$immigrant_rate > 0) {
      m <- max(1L, floor(nrow(ga_population) * ctrl$immigrant_rate))
      worst <- order(fitness_scores, decreasing = TRUE)[seq_len(m)]
      ga_population[worst, ] <- init_population(length(worst), N_EST, alpha = ctrl$init_alpha)
    }
    
    
    #Periodic evaluation checkpoint: select current best individual on training,
    # re-evaluate on full train and held-out validation (never minibatched), record convergence stats,
    # store rich diagnostics for the best solution, adapt mutation if no improvement, and optionally
    # early-stop when validation has plateaued beyond patience/min_delta thresholds.
    if (gen %% ctrl$check_every == 0L || gen == generations) {
      best_idx <- which.min(fitness_scores)
      w_best   <- ga_population[best_idx, ]
      
      eval_on <- function(w, idx)
        .ga_eval_split(w, prepped, idx, ctrl, q95_seed = q95_seed_gen)
      
      best_train <- eval_on(w_best, train_idx)   # reporte (train completo)
      best_val   <- eval_on(w_best,   val_idx)   # val completo (nunca minibatch)
      
      #capture rich diagnostics for the best individual ----
      best_details <- fitness_universal(
        w = w_best,
        prepped = prepped,
        objective = ctrl$objective,
        idx = val_idx,
        use_boot = TRUE,
        B = ctrl$bootstrap_B,
        lambda_instab = ctrl$lambda_instab,
        lambda_entropy = ctrl$lambda_entropy,
        max_w_dominance = ctrl$max_w_dominance,
        lambda_dominance = ctrl$lambda_dominance,
        mix_w_q95 = ctrl$mix_w_q95,
        mix_w_max = ctrl$mix_w_max,
        lambda_bias = ctrl$lambda_bias,
        q95_seed = q95_seed_gen,
        use_hierarchical = ctrl$use_hierarchical,
        group_by = ctrl$group_by,
        group_inner = ctrl$group_inner,
        group_inner_q = ctrl$group_inner_q,
        return_details = TRUE,
        run_id = ctrl$run_id %||% NA_character_,
        ga_seed = ctrl$ga_seed %||% NA_integer_
      )
      best_details_hist[[as.character(gen)]] <- best_details
      if (record_convergence) {
        med_train <- stats::median(fitness_scores)
        probe_val <- sapply(base::sample(seq_len(nrow(ga_population)),
                                         min(20L, nrow(ga_population))),
                            function(i) eval_on(ga_population[i, ], val_idx))
        conv <- rbind(conv, data.frame(gen=gen, best_train=best_train, med_train=med_train,
                                       best_val=best_val,   med_val=stats::median(probe_val)))
      }
      
      val_hist <- c(val_hist, best_val)
      if (best_val < best_val_so_far - 1e-12) {
        best_val_so_far <- best_val
        no_improve <- 0L
      } else {
        no_improve <- no_improve + 1L
      }
      if (no_improve > 0L) mut_rate <- min(0.30, mut_rate * (1 + 0.10 * no_improve))
      
      cat(sprintf("[%-10s] Gen %3d | Train %.6g | Val %.6g | mut=%.3f | mode=%s%s\n",
                  dist_name, gen, best_train, best_val, mut_rate, prepped_mode,
                  if (!is.null(ctrl$minibatch_frac) && identical(prepped_mode, "light"))
                    sprintf(" (mb=%d/%d)", length(train_sub), length(train_idx)) else ""))
      
      
      # prints interpretable summaries of the current best
      # weights (max weight, entropy, sparsity via nnz) and attempts to surface penalty components
      # (instability/dominance/entropy/bias) from the detailed fitness output for transparent monitoring.
      if (isTRUE(ctrl$diagnostics)) {
        ww <- as.numeric(w_best)
        mx <- suppressWarnings(max(ww))
        ent <- suppressWarnings(entropy1(ww))
        nnz <- sum(ww > 1e-4)
        # Try to surface penalty components if available
        pen <- best_details$penalties %||% NULL
        instab <- if (is.list(pen) && !is.null(pen$instability)) pen$instability else NA_real_
        domp   <- if (is.list(pen) && !is.null(pen$dominance))   pen$dominance   else NA_real_
        entp   <- if (is.list(pen) && !is.null(pen$entropy))     pen$entropy     else NA_real_
        biasp  <- if (is.list(pen) && !is.null(pen$bias))        pen$bias        else NA_real_
        cat(sprintf("[DIAG] %s | gen=%d | max_w=%.4f | entropy=%.4f | nnz=%d | pen(instab=%.4g dom=%.4g ent=%.4g bias=%.4g)\n",
                    dist_name, gen, mx, ent, nnz, instab, domp, entp, biasp))
      }
      
      if (gen >= 30L && should_stop(val_hist, patience = ctrl$patience, min_delta = ctrl$min_delta)) {
        message(sprintf("[%s] Early stopping at gen %d.", dist_name, gen))
        break
      }
    }
    
    
  }
  
  #Final selection stage: re-evaluates the entire population on validation using a
  # fixed final seed to ensure deterministic ranking, chooses the best individual, and records its
  # validation score. This avoids relying on minibatch or checkpoint noise when declaring the winner.
  q95_seed_final <- 2e6 + seed*1000 + 777
  eval_val <- function(w)
    .ga_eval_split(w, prepped, val_idx, ctrl, q95_seed = q95_seed_final)
  
  final_scores_val <- if (!is.null(cl)) {
    unlist(parallel::parLapply(
      cl,
      X = seq_len(nrow(ga_population)),
      fun = function(i, pop) eval_val(pop[i, ]),
      ga_population
    ))
  } else {
    apply(ga_population, 1, eval_val)
  }
  
  best_idx       <- which.min(final_scores_val)
  best_weights   <- ga_population[best_idx, , drop = FALSE]
  best_val_score <- min(final_scores_val)
  
  
  #Postcheck step (diagnostics only): re-evaluates the chosen best solution on full
  # train and full validation to expose potential overfitting or logging inconsistencies. Reports the
  # gap (val - train) as a quick sanity check that selection/generalization behavior matches expectations.
  # ---- POSTCHECK: re-evaluate chosen best_weights on full train/val (debug overfit vs logging) ----
  if (isTRUE(ctrl$diagnostics)) {
    tr2 <- .ga_eval_split(as.numeric(best_weights), prepped, train_idx, ctrl, q95_seed = q95_seed_final)
    va2 <- .ga_eval_split(as.numeric(best_weights), prepped,   val_idx, ctrl, q95_seed = q95_seed_final)
    cat(sprintf("[POSTCHECK] %s | Train %.6g | Val %.6g | gap=%.6g\n",
                dist_name, tr2, va2, (va2 - tr2)))
  }
  
  #Robustness test around the final GA weights: perturbs the selected best_weights
  # within the simplex neighborhood and measures the distribution of
  # resulting fitness values. This quantifies local sensitivity 
  # (fragility vs stability) and provides summary statistics for reporting.
  # ======= PERTURBATION TEST (around final GA weights) =====================
  pert_ga <- weight_perturbation_test(
    w_star = as.numeric(best_weights),
    prepped = prepped,
    sigma = 0.02,
    n = 200,
    seed = seed,
    objective = ctrl$objective,
    idx = val_idx,
    use_boot = TRUE,
    B = 0,
    lambda_instab    = ctrl$lambda_instab,
    lambda_entropy   = ctrl$lambda_entropy,
    max_w_dominance  = ctrl$max_w_dominance,
    lambda_dominance = ctrl$lambda_dominance,
    mix_w_q95        = ctrl$mix_w_q95,
    mix_w_max        = ctrl$mix_w_max,
    lambda_bias      = ctrl$lambda_bias,
    use_hierarchical = ctrl$use_hierarchical,
    group_by         = ctrl$group_by,
    group_inner      = ctrl$group_inner,
    group_inner_q    = ctrl$group_inner_q
  )
  
  cat(sprintf(
    "[PERTURB] %s | mean=%.6g sd=%.6g q95=%.6g min=%.6g max=%.6g\n",
    dist_name,
    pert_ga$fitness_mean, pert_ga$fitness_sd, pert_ga$fitness_q95,
    pert_ga$fitness_min, pert_ga$fitness_max
  ))
  
  #Warm-start persistence: stores the top-K individuals (best validation-ranked)
  # so future runs for the same family/seed can initialize from strong candidates. 
  # This accelerates convergence, improves reproducibility of best solutions, 
  # and supports multi-stage/halving pipelines.
  keep_k <- min(nrow(ga_population), max(10L, as.integer(ctrl$elitism)))
  set_warm_start(dist_name, seed, ga_population[order(final_scores_val)[seq_len(keep_k)], , drop = FALSE])
  
  #Builds a per-scenario results table by re-evaluating the chosen best weights across
  # every scenario and each sample size, extracting summary metrics (e.g., MSE, q95, max) in a tidy form.
  # This enables downstream aggregation, plotting, and reviewer-friendly reporting of scenario behavior.
  # ======================== Resultados por escenario =====================
  scenario_rows <- lapply(seq_along(prepped), function(i) {
    S <- prepped[[i]]
    dplyr::bind_rows(lapply(sample_sizes, function(n)
      summarize_scenario(S, as.numeric(best_weights), dist_name, n,
                         q95_B = ctrl$bootstrap_B, q95_seed = q95_seed_final)))
  })
  scenario_df <- dplyr::bind_rows(scenario_rows)
  
  #Aggregates scenario-level results by estimator label to report average and tail
  # performance components. Extracts the “robust” row as the primary summary for the optimized mixture,
  # capturing mean/q95/max MSE and mean bias for inclusion in the overall output table.
  agg_by_est <- scenario_df %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(mean_mse      = base::mean(mse, na.rm = TRUE),
                     q95_mse       = base::mean(mse_q95, na.rm = TRUE),
                     max_mse       = base::max(mse_max, na.rm = TRUE),
                     mean_bias     = base::mean(bias, na.rm = TRUE),
                     mean_variance = base::mean(variance, na.rm = TRUE),
                     .groups = "drop")
  robust_row <- dplyr::filter(agg_by_est, estimator == "robust")
  
  #Constructs the “overall” summary row combining metadata (distribution, objective,
  # formula strings, scenario mode), the optimized weight vector columns (w_*), and robust estimator
  # performance metrics. Also records top-K formulas for interpretability and comparative reporting.
  # =============================== OVERALL ===============================
  wcols <- as.list(round(as.numeric(best_weights), 6))
  names(wcols) <- paste0("w_", ESTIMATOR_NAMES)
  wdf <- as.data.frame(wcols, check.names = FALSE, stringsAsFactors = FALSE)
  
  topk_idx <- head(order(final_scores_val), 5L)
  base_df <- data.frame(
    distribution  = dist_name,
    objective     = ctrl$objective,
    estimator     = weights_to_formula(as.numeric(best_weights)),
    topk_formulas = paste(apply(ga_population[topk_idx, , drop = FALSE],
                                1, function(w) weights_to_formula(as.numeric(w))),
                          collapse = " || "),
    scenario_mode = prepped_mode,
    stringsAsFactors = FALSE
  )
  metrics_df <- data.frame(robust_mean_mse  = robust_row$mean_mse,
                           robust_q95_mse   = robust_row$q95_mse,
                           robust_max_mse   = robust_row$max_mse,
                           robust_mean_bias = robust_row$mean_bias,
                           stringsAsFactors = FALSE)
  overall <- cbind(base_df, wdf, metrics_df)
  
  #Final return object: exposes optimized weights and readable formula strings,
  # top-K candidate weights, full scenario table, compact overall summary, convergence history,
  # perturbation robustness stats, and best validation score (optional). Designed for downstream logging,
  # paper-ready tables, and reproducible re-analysis.
  list(weights        = as.numeric(best_weights),
       estimator_str  = weights_to_formula(as.numeric(best_weights)),
       topk_weights   = ga_population[topk_idx, , drop = FALSE],
       scenario_table = scenario_df,
       overall        = overall,
       convergence    = conv,
       perturbation   = pert_ga,
       best_val_score = if (return_val_best) best_val_score else NA_real_)
  
}


# =========================== GA con K-fold CV (K=3) ======================
# K-fold cross-validation wrapper (default K=3) around the single-split GA.
# Builds a shared scenario universe (optionally subsampled/exported), runs GA per fold with warm-start
# accumulation, selects the best fold by validation score, and optionally performs a final retrain on
# all scenarios using the CV warm bank for a stronger final estimator.

evolve_universal_estimator_per_family_cv <- function(dist_name,
                                                     dist_param_grid,
                                                     sample_sizes = c(300, 500, 1000, 2000, 5000),
                                                     num_samples = 80,
                                                     pop_size = 100,
                                                     generations_per_fold = 35,
                                                     seed = 101,
                                                     objective = "q95",
                                                     use_parallel = TRUE,
                                                     k_folds = 3,
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
                                                     lambda_entropy = 0.015,
                                                     max_w_dominance = 0.80,
                                                     lambda_dominance = 0.05,
                                                     lambda_bias = 0.05,
                                                     minibatch_frac = NULL,
                                                     minibatch_min  = 12L,
                                                     crn_env = NULL,
                                                     fam_key = NULL,
                                                     #scenario subsampling (works for FULL and LIGHT)
                                                     scenario_frac = 1.0,
                                                     scenario_seed = NULL,
                                                     subset_tag = NULL,
                                                     #forzar full/light desde fuera
                                                     force_scenario_mode = NULL,
                                                     #override scenario universe/subset from launcher
                                                     scenario_universe = NULL,
                                                     scenario_subset_override = NULL,
                                                     run_dir = NULL,
                                                     stage_name = NULL) {
  # Enables lightweight debug instrumentation via a global option without changing
  # GA selection logic. This allows turning on richer prints and postchecks in development runs while
  # keeping default behavior identical in production/cluster runs for reproducible comparisons.
  # Extra: debug instrumentation (does not change selection)
  diagnostics <- isTRUE(getOption("GA_DIAGNOSTICS", FALSE))
  
  #Validates CV/GA configuration and seeds the run deterministically. This ensures
  # folds, subsampling, and all stochastic components remain reproducible. Also computes scenario_mode
  # (full/light), either forced by caller or auto-selected based on computational load expectations.
  stopifnot(length(sample_sizes) >= 1L, pop_size >= 2L, generations_per_fold >= 1L, k_folds >= 2L)
  seed <- .ensure_seed(seed, fallback = 101L)
  set.seed(seed)
  
  scenario_mode <- if (!is.null(force_scenario_mode)) {
    match.arg(force_scenario_mode, c("full","light"))
  } else {
    .ga_pick_mode(sample_sizes, num_samples, pop_size, generations_per_fold)
  }
  
  
  
  
  
  # Defines the scenario subset used for CV. If an override is provided, 
  # it is used directly (with scenario_id ensured). Otherwise,
  # it builds/uses the full scenario universe and selects a reproducible
  # subset (scenario_frac, min_n>=k_folds).
  # Produces tags/seeds for tracking.
  # ---- scenario subsampling (FULL or LIGHT) -----------------
  if (!is.null(scenario_subset_override)) {
    scenario_subset <- as.data.frame(scenario_subset_override)
    if (is.null(scenario_subset$scenario_id)) scenario_subset <- add_scenario_ids(scenario_subset)
    sc_seed <- scenario_seed %||% seed
    ss_tag  <- subset_tag %||% paste0(dist_name, "__", scenario_mode, "__override")
  } else {
    sc_all <- if (!is.null(scenario_universe)) as.data.frame(scenario_universe) else build_scenarios(scenario_mode)
    if (is.null(sc_all$scenario_id)) sc_all <- add_scenario_ids(sc_all)
    sc_seed <- scenario_seed %||% seed
    ss_tag  <- subset_tag %||% paste0(dist_name, "__", scenario_mode, "__frac", format(scenario_frac, digits = 3))
    scenario_subset <- pick_scenario_subset(
      sc_all,
      scenario_frac = scenario_frac,
      min_n = k_folds,
      seed = sc_seed,
      subset_tag = ss_tag
    )
  }
  
  #  writes the 
  # exact scenario subset used (scenarios_used_*.csv) and its 
  # difficulty summary to run_dir, then registers both in an index.
  # # This makes the CV experiment auditable and enables re-running
  # identical scenario universes later.
  # # ----per-stage exports (scenarios_used + difficulty) ---------------
  if (!is.null(run_dir)) {
    st <- stage_name %||% ss_tag
    st <- gsub("[^A-Za-z0-9]+", "_", st)
    out_sc <- file.path(run_dir, sprintf("scenarios_used_%s_%s.csv", dist_name, st))
    safe_write_csv(scenario_subset, out_sc)
    out_df <- file.path(run_dir, sprintf("difficulty_%s_%s.csv", dist_name, st))
    safe_write_csv(stage_difficulty(scenario_subset), out_df)
    # Index entries (n_candidates: NA here unless caller sets it via stage_name conventions)
    index_add(dist_name, st, out_sc, n_candidates = NA_integer_, n_scenarios = nrow(scenario_subset), seed = sc_seed)
    index_add(dist_name, st, out_df, n_candidates = NA_integer_, n_scenarios = nrow(scenario_subset), seed = sc_seed)
  }
  
  # #Prepares (and caches) the reference scenario simulations
  # for this family using the selected subset. This prepped_ref 
  # is shared across folds so CV compares apples-to-apples: 
  # each fold uses consistent scenario data and CRN settings, 
  # with only index splits changing between train/val.
  prepped_ref <- prep_scenarios(
    dist_name, dist_param_grid,
    sample_sizes  = sample_sizes,
    num_samples   = num_samples,
    scenario_mode = scenario_mode,
    seed          = seed,
    crn_env       = crn_env,
    fam_key         = if (is.null(fam_key)) dist_name else fam_key,
    scenario_subset = scenario_subset,
    subset_tag      = ss_tag,
    use_cache       = TRUE
  )
  all_idx <- seq_along(prepped_ref)
  if (length(all_idx) < k_folds) {
    stop(sprintf("k_folds=%d es mayor que #escenarios=%d.", k_folds, length(all_idx)))
  }
  
  #Creates deterministic fold assignments and initializes fold result containers.
  #warm_bank accumulates strong solutions from previous folds to 
  # warm-start later folds, which reduces compute and improves stability. 
  # Each fold is evaluated on its own held-out scenario subset.
  folds <- make_folds(length(all_idx), k = k_folds, seed = seed)
  fold_results <- vector("list", k_folds)
  warm_bank <- NULL
  
  # Helper that runs the single-split GA for a given fold (train_idx/val_idx) with a
  # fold-specific seed and optional warm-start matrix. It forwards all key hyperparameters and injects
  # prepped_ref to guarantee scenarios are identical across folds; retrain_all controls scoring output.
  run_fold <- function(train_idx, val_idx, seed_fold, warm_mat, gens, pop_sz, retrain_all = FALSE) {
    evolve_universal_estimator_per_family(
      dist_name          = dist_name,
      dist_param_grid    = dist_param_grid,
      sample_sizes       = sample_sizes,
      num_samples        = num_samples,
      pop_size           = pop_sz,
      generations        = gens,
      seed               = seed_fold,
      objective          = objective,
      use_parallel       = use_parallel,
      train_idx          = train_idx,
      val_idx            = val_idx,
      lambda_instab      = lambda_instab,
      bootstrap_B        = bootstrap_B,
      t_size             = t_size,
      elitism            = elitism,
      immigrant_rate     = immigrant_rate,
      mutation_rate_init = mutation_rate_init,
      init_alpha         = init_alpha,
      alpha_mut          = alpha_mut,
      check_every        = check_every,
      patience           = patience,
      min_delta          = min_delta,
      warm_start_matrix  = warm_mat,
      return_val_best    = !retrain_all,
      mix_w_q95          = mix_w_q95,
      mix_w_max          = mix_w_max,
      lambda_entropy     = lambda_entropy,
      max_w_dominance    = max_w_dominance,
      lambda_dominance   = lambda_dominance,
      lambda_bias        = lambda_bias,
      minibatch_frac     = minibatch_frac,
      minibatch_min      = minibatch_min,
      crn_env            = crn_env,
      fam_key            = if (is.null(fam_key)) dist_name else fam_key,
      
      force_scenario_mode = force_scenario_mode,
      
      prepped_override = prepped_ref,
      diagnostics = diagnostics
    )
  }
  
  #Executes the K folds: for each fold, defines validation indices and complementary
  # training indices, runs the GA, stores fold outputs, and appends the fold’s optimized weights to the
  # warm_bank for subsequent folds. This progressively seeds later folds with already-good candidates.
  for (k in seq_len(k_folds)) {
    val_idx   <- sort(unlist(folds[[k]]))
    train_idx <- sort(setdiff(all_idx, val_idx))
    cat(sprintf("\n[%s] Fold %d/%d | Train=%d | Val=%d | mode=%s\n",
                dist_name, k, k_folds, length(train_idx), length(val_idx), scenario_mode))
    res_k <- run_fold(train_idx, val_idx, seed + k, warm_bank,
                      gens = generations_per_fold, pop_sz = pop_size)
    fold_results[[k]] <- res_k
    warm_bank <- rbind(warm_bank, res_k$weights)
  }
  
  # Selects the best-performing fold by its validation score and extracts the
  # corresponding weights as the CV champion. This provides a principled choice when folds vary in
  # difficulty, and sets up the final retrain step (if enabled) with strong warm-start candidates.
  val_scores <- vapply(fold_results, function(x) x$best_val_score %||% Inf, numeric(1))
  best_fold  <- which.min(val_scores)
  best_weights_cv <- fold_results[[best_fold]]$weights
  
  #Final retrain option: runs one more GA pass on ALL scenarios using an expanded
  # warm-start matrix (warm_bank + best CV weights). This typically yields a stronger final solution
  # for deployment/publication while CV still provides unbiased model selection via fold validation.
  # -------- FINAL OUTPUT (correcto) --------
  if (isTRUE(final_retrain)) {
    cat(sprintf("\n[%s] Final retrain on ALL scenarios (warm-start)\n", dist_name))
    
    res_final <- run_fold(
      train_idx   = all_idx,
      val_idx     = all_idx,
      seed_fold   = seed + 999,
      warm_mat    = rbind(warm_bank, best_weights_cv),
      gens        = max(40L, as.integer(generations_per_fold)),
      pop_sz      = max(60L, as.integer(pop_size)),
      retrain_all = TRUE
    )
    
    return(list(
      fold_results = fold_results,
      final        = res_final,
      best_fold    = best_fold,
      val_scores   = val_scores
    ))
  } else {
    #  Alternative output path when final_retrain is disabled: returns the best fold’s
    # result directly as the final model. This is faster and preserves an unbiased validation estimate,
    # but may yield slightly weaker weights than training once on the full scenario set with warm-starts.
    # NO retrain: devolvemos el mejor fold
    return(list(
      fold_results = fold_results,
      final        = fold_results[[best_fold]],
      best_fold    = best_fold,
      val_scores   = val_scores
    ))
  }
  
  
}

# This is just a tag... if the module is successfully 
#loaded this tag will be printed out in the console
mark_module_done("05_ga_core.R")


  
