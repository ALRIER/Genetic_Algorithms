

# =========================== SETUP & PACKAGES ============================

required_packages <- c(
  "modeest",   # modos (hsm, parzen)
  "statmod",   # rinvgauss
  "gtools",    # rdirichlet
  "parallel",  # PSOCK cluster
  "MASS",      # rlm (Huber/Bisquare)
  "dplyr", "tibble", "readr"  # utilidades de datos
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
for (pkg in required_packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Small helper used later: x %||% y returns x if not NULL else y
`%||%` <- function(x, y) if (!is.null(x)) x else y

# =========================== WORKING DIRECTORY ===========================
setwd("C:/Users/AlvaroRivera-Eraso/Documents/simulating distributions and estimators/Estimator_Results")
options(stringsAsFactors = FALSE)

# ============================ RNG & REPRODUCIBILIDAD =====================
## Usamos L'Ecuyer-CMRG para reproducibilidad en paralelo y CRN
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
invisible(runif(1))

## Generador de índices CRN por (familia, n): mismos draws para TODAS las configs
make_crn_indices <- function(families, sample_sizes, B_boot = 100L, num_samples = 50L) {
  crn <- new.env(parent = emptyenv())
  ctr <- 1L
  for (fam in families) {
    for (n in sample_sizes) {
      ## Fijamos sub-semilla estable por par (fam, n)
      s_backup <- if (exists(".Random.seed", inherits = FALSE)) .Random.seed else NULL
      set.seed(1e6 + ctr)
      
      ## Índices/sortids que luego tu pipeline puede reutilizar
      sim_idx <- sample.int(1e9, num_samples, replace = FALSE)
      
      assign(paste(fam, n, sep = "__"),
             list(sim_idx = sim_idx),
             envir = crn)
      
      ## Restaurar estado RNG global (si existía)
      if (!is.null(s_backup)) .Random.seed <- s_backup
      ctr <- ctr + 1L
    }
  }
  crn
}

## ============================ SCENARIO BUILDERS ============================
## Full grid of scenarios (original behaviour, expanded)
build_scenarios_full <- function(
    rates  = c(0.00, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30),
    scales = c(1.5, 3, 4.5, 6, 9, 12, 15, 20, 30),
    types  = c("upper_tail", "lower_tail", "symmetric_t", "point_mass")
) {
  expand.grid(
    contamination_rate = rates,
    outlier_scale_mad  = scales,
    contamination_type = types,
    stringsAsFactors = FALSE
  )
}

## Light grid of scenarios (fast proxy)
## Cover clean, mild, and extreme; symmetric + asymmetric
build_scenarios_light <- function() {
  expand.grid(
    contamination_rate = c(0.00, 0.02, 0.20),   # clean, mild, extreme
    outlier_scale_mad  = c(3, 9, 20),           # small, med, huge
    contamination_type = c("upper_tail", "symmetric_t"),
    stringsAsFactors = FALSE
  )
}

## Cache para warm-start entre configuraciones vecinas (misma familia/seed)
.init_pop_cache <- new.env(parent = emptyenv())
get_warm_start <- function(fam, seed) {
  key <- paste(fam, seed, sep = "::")
  if (exists(key, envir = .init_pop_cache)) get(key, envir = .init_pop_cache) else NULL
}
set_warm_start <- function(fam, seed, population) {
  key <- paste(fam, seed, sep = "::")
  assign(key, population, envir = .init_pop_cache)
}

## Early stopping con umbral de mejora mínima (min_delta en validación)
should_stop <- function(val_hist, patience = 2L, min_delta = 0.005) {
  L <- length(val_hist)
  if (L < patience + 1L) return(FALSE)
  best_prev   <- min(val_hist[1:(L - patience)])
  best_recent <- min(val_hist[(L - patience + 1L):L])
  improvement <- (best_prev - best_recent) / (abs(best_prev) + 1e-12)
  improvement < min_delta
}

## Normalización defensiva al simplex (para vectores de pesos)
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





# ============================ UTILIDADES I/O =============================

RUN_TS <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
mkdirp  <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

write_manifest <- function(dir_out, meta, extra=NULL) {
  mkdirp(dir_out)
  scenario_mode_line <- meta$scenario_mode %||% "full"
  lines <- c(
    sprintf("# RUN_ID: %s", meta$run_id),
    sprintf("timestamp: %s", as.character(Sys.time())),
    sprintf("families: %s", paste(meta$families, collapse=", ")),
    sprintf("k_folds: %s", meta$k_folds),
    sprintf("objective: %s", meta$objective),
    sprintf("seeds: %s", meta$seeds %||% meta$seed),
    sprintf("pop_sizes: %s", paste(meta$pop_sizes, collapse=", ")),
    sprintf("generations_per_fold: %s", meta$generations_per_fold),
    sprintf("final_retrain: %s", meta$final_retrain),
    sprintf("bootstrap_B: %s", meta$bootstrap_B),
    sprintf("lambda_instab (default): %s", meta$lambda_instab_default),
    sprintf("scenario_mode: %s", scenario_mode_line)
  )
  if (!is.null(extra)) {
    lines <- c(lines, "", "## EXTRA", capture.output(print(extra)))
  }
  writeLines(lines, file.path(dir_out, "MANIFEST.txt"))
}

# (El bloque de Warm Bank persistente en RDS fue removido porque no se usa)

# ============================ RNG & FAMILIAS =============================

rexwald <- function(n, mu = 1, lambda = 1, rate = 1) {
  stopifnot(n >= 0, is.finite(mu), is.finite(lambda), is.finite(rate))
  if (n == 0) return(numeric(0))
  if (lambda <= 0 || rate <= 0) stop("rexwald: 'lambda' and 'rate' must be > 0")
  T0 <- rexp(n, rate = rate)
  W  <- statmod::rinvgauss(n, mean = mu, shape = lambda)
  T0 + W
}

rexgaussian <- function(n, mu = 0, sigma = 1, tau = 1) {
  stopifnot(n >= 0, is.finite(mu), is.finite(sigma), is.finite(tau))
  if (n == 0) return(numeric(0))
  if (sigma <= 0 || tau <= 0) stop("rexgaussian: 'sigma' and 'tau' must be > 0")
  rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1 / tau)
}

# Robust generator with defensive parameter handling and clear errors
generate_population <- function(distribution, n, params = list()) {
  if (!is.character(distribution) || length(distribution) != 1L)
    stop("generate_population: 'distribution' must be a single string.")
  if (!is.numeric(n) || length(n) != 1L || n < 0)
    stop("generate_population: 'n' must be a single non-negative number.")
  if (n == 0) return(numeric(0))
  params <- as.list(params)
  
  getp <- function(name, default = NULL, required = FALSE) {
    if (!is.null(params[[name]])) return(params[[name]])
    if (required) stop(sprintf("generate_population: missing required parameter '%s' for %s.",
                               name, distribution))
    default
  }
  
  out <- switch(tolower(distribution),
                normal = {
                  mean <- getp("mean", 0, required = FALSE)
                  sd   <- getp("sd",   1, required = FALSE)
                  if (!is.finite(sd) || sd <= 0) stop("normal: 'sd' must be > 0 and finite.")
                  rnorm(n, mean = mean, sd = sd)
                },
                lognormal = {
                  meanlog <- getp("meanlog", 0, required = FALSE)
                  sdlog   <- getp("sdlog",   1, required = FALSE)
                  if (!is.finite(sdlog) || sdlog <= 0) stop("lognormal: 'sdlog' must be > 0 and finite.")
                  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
                },
                weibull = {
                  shape <- getp("shape", NULL, required = TRUE)
                  scale <- getp("scale", NULL, required = TRUE)
                  if (!is.finite(shape) || shape <= 0) stop("weibull: 'shape' must be > 0 and finite.")
                  if (!is.finite(scale) || scale <= 0) stop("weibull: 'scale' must be > 0 and finite.")
                  rweibull(n, shape = shape, scale = scale)
                },
                invgauss = {
                  mean  <- getp("mean",  NULL, required = TRUE)
                  shape <- getp("shape", NULL, required = TRUE)
                  if (!is.finite(mean)  || mean  <= 0) stop("invgauss: 'mean' must be > 0 and finite.")
                  if (!is.finite(shape) || shape <= 0) stop("invgauss: 'shape' must be > 0 and finite.")
                  statmod::rinvgauss(n, mean = mean, shape = shape)
                },
                exgaussian = {
                  mu    <- getp("mu",    0, required = FALSE)
                  sigma <- getp("sigma", 1, required = FALSE)
                  tau   <- getp("tau",   1, required = FALSE)
                  if (!is.finite(sigma) || sigma <= 0) stop("exgaussian: 'sigma' must be > 0 and finite.")
                  if (!is.finite(tau)   || tau   <= 0) stop("exgaussian: 'tau' must be > 0 and finite.")
                  rexgaussian(n, mu = mu, sigma = sigma, tau = tau)
                },
                exwald = {
                  mu     <- getp("mu",     1, required = FALSE)
                  lambda <- getp("lambda", 1, required = FALSE)
                  rate   <- getp("rate",   1, required = FALSE)
                  if (!is.finite(lambda) || lambda <= 0) stop("exwald: 'lambda' must be > 0 and finite.")
                  if (!is.finite(rate)   || rate   <= 0) stop("exwald: 'rate' must be > 0 and finite.")
                  rexwald(n, mu = mu, lambda = lambda, rate = rate)
                },
                stop(sprintf("Unsupported distribution: %s", distribution))
  )
  
  # Final guard: drop non-finite draws if any numerical quirk sneaks in
  out[is.finite(out)]
}



# ============================ PARAMETER GRIDS ============================

param_grids <- list(
  normal     = expand.grid(mean = c(1, 5), sd = c(1, sqrt(2), 2, 3)),
  lognormal  = expand.grid(meanlog = c(0.5, 1), sdlog = c(0.25, 0.5, 1)),
  weibull    = expand.grid(shape = c(0.5, 1, 2, 3), scale = c(0.5, 1, 2)),
  invgauss   = expand.grid(mean = c(0.5, 1, 2), shape = c(0.5, 1, 2)),
  exgaussian = expand.grid(mu = c(0, 1, 2), sigma = c(0.5, 1, 1.5), tau = c(0.25, 0.5, 1)),
  exwald     = expand.grid(mu = c(1, 1.5, 2), lambda = c(0.5, 1, 2), rate = c(0.25, 0.5, 1))
)

# ======================== SET DE ESTIMADORES ============================

ESTIMATOR_NAMES <- c("mean","median","trimmed20","harmonic","geometric",
                     "mode_hsm","mode_parzen","trimean","huber","biweight")
N_EST <- length(ESTIMATOR_NAMES)

harmonic_mean_safe <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(NA_real_)
  xp <- x[is.finite(x) & x > 0]
  if (length(xp) == 0L) return(stats::median(x, na.rm = TRUE))
  length(xp) / sum(1 / xp)
}

geometric_mean_safe <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(NA_real_)
  xp <- x[is.finite(x) & x > 0]
  if (length(xp) == 0L) return(stats::median(x, na.rm = TRUE))
  exp(base::mean(log(xp)))
}


mode_hsm_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch(modeest::hsm(x), error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

mode_parzen_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch(as.numeric(modeest::mlv(x, method = "parzen")),
                  error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

trimean_safe <- function(x) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 8, na.rm = TRUE)
  (qs[1] + 2 * qs[2] + qs[3]) / 4
}

huber_mean_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch({
    fit <- MASS::rlm(x ~ 1, psi = MASS::psi.huber, scale.est = "MAD",
                     maxit = 50, na.action = na.omit)
    as.numeric(coef(fit)[1])
  }, error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

biweight_mean_safe <- function(x) {
  x <- as.numeric(x)
  out <- tryCatch({
    fit <- MASS::rlm(x ~ 1, psi = MASS::psi.bisquare, scale.est = "MAD",
                     maxit = 50, na.action = na.omit)
    as.numeric(coef(fit)[1])
  }, error = function(e) NA_real_)
  if (!is.finite(out)) stats::median(x, na.rm = TRUE) else out
}

# Acepta pesos en cualquier escala y los normaliza (usa .normalize_simplex global)
custom_estimator <- function(sample, weights) {
  stopifnot(length(weights) == N_EST)
  w <- .normalize_simplex(as.numeric(weights))
  x <- as.numeric(sample)
  comps <- c(
    base::mean(x, na.rm = TRUE),
    stats::median(x, na.rm = TRUE),
    base::mean(x, trim = 0.20, na.rm = TRUE),
    harmonic_mean_safe(x),
    geometric_mean_safe(x),
    mode_hsm_safe(x),
    mode_parzen_safe(x),
    trimean_safe(x),
    huber_mean_safe(x),
    biweight_mean_safe(x)
  )
  comps[!is.finite(comps)] <- stats::median(x, na.rm = TRUE)
  sum(w * comps)
}


weights_to_formula <- function(weights) {
  stopifnot(length(weights) == N_EST)
  w <- round(.normalize_simplex(as.numeric(weights)), 3)
  paste0(paste0(w, "*", ESTIMATOR_NAMES), collapse = " + ")
}







































# ======= POBLACIÓN, CROSSOVER y MUTACIÓN =================================

init_population <- function(size, n_estimators, alpha = 1, warm_start = NULL, jitter = 0.15) {
  stopifnot(size >= 1L, n_estimators >= 1L, is.finite(alpha), alpha > 0)
  # Generación base por Dirichlet
  pop <- matrix(gtools::rdirichlet(size, rep(alpha, n_estimators)),
                nrow = size, byrow = TRUE)
  # Warm-start opcional
  if (!is.null(warm_start)) {
    ws <- as.matrix(warm_start)
    # Ajuste defensivo por si vienen columnas de más/menos
    if (ncol(ws) != n_estimators) {
      if (ncol(ws) > n_estimators) ws <- ws[, seq_len(n_estimators), drop = FALSE]
      if (ncol(ws) < n_estimators) {
        add <- matrix(0, nrow = nrow(ws), ncol = n_estimators - ncol(ws))
        ws <- cbind(ws, add)
      }
    }
    # Normaliza cada fila del warm-start por seguridad
    ws <- t(apply(ws, 1, .normalize_simplex))
    k <- min(nrow(ws), size)
    # Inyecta k individuos del warm-start (con leve jitter para diversidad)
    if (k > 0) {
      eps <- matrix(gtools::rdirichlet(k, rep(alpha + jitter, n_estimators)),
                    nrow = k, byrow = TRUE)
      pop[seq_len(k), ] <- t(apply((ws[seq_len(k), , drop = FALSE] + eps) / 2, 1, .normalize_simplex))
    }
  }
  # Normaliza por si acaso
  pop <- t(apply(pop, 1, .normalize_simplex))
  pop
}

crossover <- function(parent1, parent2) {
  stopifnot(length(parent1) == length(parent2))
  alpha <- stats::runif(1)
  child <- alpha * as.numeric(parent1) + (1 - alpha) * as.numeric(parent2)
  .normalize_simplex(child)
}

mutate_weights <- function(weights, mutation_rate = 0.1, alpha_mut = 1) {
  w <- .normalize_simplex(weights)
  mutation_rate <- min(max(mutation_rate, 0), 1)   # clip a [0,1]
  if (stats::runif(1) < mutation_rate) {
    alpha_mut <- max(alpha_mut, 1e-6)  # robustez numérica
    perturb <- as.numeric(gtools::rdirichlet(1, rep(alpha_mut, length(w))))
    w <- (w + perturb) / 2
  }
  .normalize_simplex(w)
}

# ====================== CONTAMINACIÓN REALISTA ===========================

inject_outliers_realistic <- function(sample,
                                      contamination_rate = 0.05,
                                      outlier_scale_mad = 12,
                                      type = c("upper_tail", "symmetric_t", "lower_tail", "point_mass"),
                                      df_t = 3) {
  # Añadimos "lower_tail" para alinear con los escenarios "full"; mantenemos "point_mass" por compatibilidad
  type <- match.arg(type)
  n <- length(sample)
  k <- ceiling(n * contamination_rate)
  if (k <= 0) return(sample)
  
  med  <- stats::median(sample)
  madv <- stats::mad(sample, constant = 1)
  if (!is.finite(madv) || madv <= 0) madv <- stats::sd(sample)
  
  idx <- base::sample.int(n, size = k, replace = FALSE)
  
  if (type == "upper_tail") {
    bump <- abs(stats::rnorm(k))
    outlier_values <- med + outlier_scale_mad * madv * bump
  } else if (type == "symmetric_t") {
    bump <- stats::rt(k, df = df_t)
    outlier_values <- med + outlier_scale_mad * madv * bump
  } else if (type == "lower_tail") {
    bump <- abs(stats::rnorm(k))
    outlier_values <- med - outlier_scale_mad * madv * bump
  } else {  # point_mass
    outlier_values <- rep(med + outlier_scale_mad * madv, k)
  }
  
  sample[idx] <- outlier_values
  sample
}

# ============================ ESCENARIOS & DATA ==========================

## build_scenarios() is now just a dispatcher to light/full
build_scenarios <- function(mode = c("full","light")) {
  mode <- match.arg(mode)
  if (mode == "full") {
    build_scenarios_full()
  } else {
    build_scenarios_light()
  }
}

compute_components_vector <- function(x) {
  # defensivo: quita no finitos si aparecieran
  x <- x[is.finite(x)]
  comps <- c(
    base::mean(x, na.rm = TRUE),
    stats::median(x, na.rm = TRUE),
    base::mean(x, trim = 0.20, na.rm = TRUE),
    harmonic_mean_safe(x), geometric_mean_safe(x),
    mode_hsm_safe(x), mode_parzen_safe(x),
    trimean_safe(x), huber_mean_safe(x), biweight_mean_safe(x)
  )
  comps[!is.finite(comps)] <- stats::median(x, na.rm = TRUE)
  comps
}

# --- Helper: evaluar expr con semilla local sin afectar RNG global -------
.seed_scope <- function(seed, expr) {
  s_backup_exists <- exists(".Random.seed", inherits = FALSE)
  if (s_backup_exists) s_backup <- .Random.seed
  on.exit({
    if (s_backup_exists) .Random.seed <<- s_backup
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# --- PREP con CRN: usa sim-seeds por (familia, n) para muestrear reproducible
# 'crn_env' viene de make_crn_indices(); key = paste(family, n, sep="__")
prep_scenarios <- function(dist_name, dist_param_grid,
                           sample_sizes = c(50, 100, 300, 500, 1000, 2000),
                           num_samples = 80,
                           scenario_mode = c("full","light"),
                           seed = 123,
                           crn_env = NULL,
                           fam_key = NULL) {
  set.seed(seed)
  scenario_mode <- match.arg(scenario_mode)
  scenarios <- build_scenarios(scenario_mode)
  
  # 1) Población base combinada (grande) para aproximar "verdadero" mean
  combined_population <- unlist(lapply(seq_len(nrow(dist_param_grid)), function(i) {
    param <- as.list(dist_param_grid[i, , drop = FALSE])
    tryCatch(generate_population(dist_name, 1e5, param), error = function(e) NULL)
  }))
  combined_population <- combined_population[is.finite(combined_population)]
  if (length(combined_population) == 0) stop("Failed to generate combined population.")
  true_mean <- base::mean(combined_population)
  
  # 2) Para cada escenario y tamaño muestral: generar 'num_samples' replicaciones
  prepped <- lapply(seq_len(nrow(scenarios)), function(i) {
    sc <- scenarios[i, ]
    
    # asegurar que el tipo de contaminación sea character (no factor)
    sc$contamination_type <- as.character(sc$contamination_type)
    
    components_by_size <- lapply(sample_sizes, function(n) {
      # ¿Tenemos CRN para este (familia, n)?
      crn_key <- if (!is.null(fam_key)) paste(fam_key, n, sep = "__") else NULL
      sim_seeds <- NULL
      if (!is.null(crn_env) && !is.null(crn_key) && exists(crn_key, envir = crn_env)) {
        obj <- get(crn_key, envir = crn_env)
        sim_seeds <- obj$sim_idx
        # Asegurar longitud
        if (length(sim_seeds) < num_samples) {
          sim_seeds <- rep(sim_seeds, length.out = num_samples)  # reciclar determinísticamente
        } else {
          sim_seeds <- sim_seeds[seq_len(num_samples)]
        }
      }
      
      # Construir matriz de componentes (num_samples x N_EST)
      comps_mat <- matrix(NA_real_, nrow = num_samples, ncol = N_EST)
      for (r in seq_len(num_samples)) {
        # Muestreo reproducible con (familia, n, r) si hay CRN
        if (!is.null(sim_seeds)) {
          s <- .seed_scope(sim_seeds[r], base::sample(combined_population, n))
        } else {
          # si la población combinada es más pequeña que n, usa replace=TRUE por seguridad
          replace_flag <- length(combined_population) < n
          s <- base::sample(combined_population, n, replace = replace_flag)
        }
        
        # Inyectar outliers si aplica
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
    
    # --- construir el objeto de escenario y anexar el atributo scenario_mode ---
    S_i <- list(
      components_by_size = components_by_size,
      true_mean = true_mean,
      scenario  = sc
    )
    attr(S_i, "scenario_mode") <- scenario_mode  # <- atributo por-escenario
    S_i
  })
  
  # Atributos a nivel de lista completa (se mantienen también)
  attr(prepped, "sample_sizes")   <- sample_sizes
  attr(prepped, "scenario_mode")  <- scenario_mode
  prepped
}



# =============================== FITNESS =================================

# --- Huber loss para penalizar outliers de manera robusta ----------------
huber_loss <- function(r, delta = 1.0) {
  r <- as.numeric(r)
  a <- base::abs(r)
  ifelse(a <= delta, 0.5 * r^2, delta * (a - 0.5 * delta))
}

# --- Ponderación por dificultad del escenario ----------------------------
# Nota: toma en cuenta el tipo de contaminación (incluyendo lower_tail)
scenario_weight <- function(sc_row) {
  cr  <- sc_row$contamination_rate
  osm <- sc_row$outlier_scale_mad
  typ <- as.character(sc_row$contamination_type)
  base_w  <- 1 + 2 * cr + 0.03 * osm
  type_adj <- switch(typ,
                     "upper_tail"   = 1.05,
                     "lower_tail"   = 1.05,
                     "symmetric_t"  = 1.10,
                     1.00)
  base_w * type_adj
}

# --- Bootstrap q95 con opción determinística (seed o índices) -----------
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
    stats::runif(1)  # fuerza existencia de .Random.seed
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

fitness_universal <- function(w, prepped,
                              objective = c("mean","q95","max","mixed"),
                              idx = NULL, use_boot = TRUE, B = 300,
                              lambda_instab = 0.0,
                              lambda_entropy = 0.015,
                              max_w_dominance = 0.80,
                              lambda_dominance = 0.05,
                              mix_w_q95 = 0.7,
                              mix_w_max = 0.3,
                              lambda_bias = 0.05,
                              huber_delta = 1.0,
                              q95_seed = NULL,
                              q95_boot_indices = NULL) {
  w <- .normalize_simplex(as.numeric(w))
  objective <- match.arg(objective)
  
  # Subconjunto de escenarios si se pasa 'idx'
  Slist <- if (is.null(idx)) prepped else prepped[idx]
  if (length(Slist) == 0L) return(Inf)
  
  # Si es "light" y no se pasó B explícito, usa un B más chico
  if (!is.null(attr(prepped, "scenario_mode")) &&
      attr(prepped, "scenario_mode") == "light" &&
      missing(B)) {
    B <- 60L
  }
  
  # --- Pérdidas por escenario (Huber) ---
  losses_by_scen <- vapply(Slist, function(S) {
    ls_per_size <- vapply(S$components_by_size, function(C) {
      est <- as.vector(C %*% w)
      r   <- est - S$true_mean
      base::mean(huber_loss(r, delta = huber_delta), na.rm = TRUE)
    }, numeric(1))
    base::mean(ls_per_size) * scenario_weight(S$scenario)
  }, numeric(1))
  
  finite_mask <- is.finite(losses_by_scen)
  if (!all(finite_mask)) {
    losses_by_scen <- losses_by_scen[finite_mask]
    if (length(losses_by_scen) == 0L) return(Inf)
  }
  
  # --- Sesgo medio por escenario (evita stats::mean) ---
  bias_by_scen <- vapply(Slist, function(S) {
    base::mean(vapply(S$components_by_size, function(C) {
      est <- as.vector(C %*% w)
      base::mean(est, na.rm = TRUE) - S$true_mean
    }, numeric(1)))
  }, numeric(1))
  mean_abs_bias <- base::mean(base::abs(bias_by_scen), na.rm = TRUE)
  
  # --- Métrica base del objetivo ---
  q95_val  <- if (isTRUE(use_boot)) {
    q95_boot(losses_by_scen, B = B, rng_seed = q95_seed, boot_indices = q95_boot_indices)
  } else {
    as.numeric(stats::quantile(losses_by_scen, 0.95, names = FALSE, type = 8))
  }
  max_val  <- base::max(losses_by_scen)
  mean_val <- base::mean(losses_by_scen)
  
  base_obj <- switch(objective,
                     mean   = mean_val,
                     max    = max_val,
                     q95    = q95_val,
                     mixed  = mix_w_q95 * q95_val + mix_w_max * max_val)
  
  # --- Penalizaciones ---
  penalty_instab <- lambda_instab * stats::sd(losses_by_scen, na.rm = TRUE)
  
  eps <- 1e-12
  entropy <- -sum(w * base::log(w + eps))
  entropy_norm <- entropy / base::log(length(w))
  pen_entropy <- lambda_entropy * (1 - entropy_norm)
  
  dom_excess <- base::max(0, base::max(w) - max_w_dominance)
  pen_dominance <- lambda_dominance * dom_excess
  
  pen_bias <- lambda_bias * mean_abs_bias
  
  base_obj + penalty_instab + pen_entropy + pen_dominance + pen_bias
}









# ======================== SELECCIÓN / PARTICIONES ========================

tournament_select <- function(scores, pop_size, t_size = 2L, replace = TRUE) {
  stopifnot(is.numeric(scores), length(scores) >= 2L, pop_size >= 2L, t_size >= 2L)
  n <- length(scores)
  t_size  <- min(t_size, n)             # seguridad si t_size > n
  n_pairs <- max(1L, pop_size %/% 2L)   # número de ganadores a muestrear
  winners <- integer(n_pairs)
  
  for (i in seq_len(n_pairs)) {
    cand <- base::sample.int(n, size = t_size, replace = FALSE)
    # rompe empates de manera estable (mínimo índice si hay empate)
    local_best <- cand[which.min(scores[cand])]
    winners[i] <- local_best
  }
  
  if (!replace) winners <- unique(winners)
  if (!replace && length(winners) < n_pairs) {
    ord  <- base::order(scores)         # completa con los mejores restantes
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
      if (!is.null(s_bkp)) .Random.seed <<- s_bkp
      else rm(".Random.seed", inherits = FALSE)
    }, add = TRUE)
    set.seed(seed)
  }
  ids <- base::sample(seq_len(n))
  # Partición balanceada (diferencia de tamaños <= 1)
  split(ids, rep(seq_len(k), length.out = n))
}

# ================ WARM-STARTS DINÁMICOS (según N_EST) ====================
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

# =========================== MÉTRICAS / RESUMEN ==========================
summarize_scenario <- function(S, weights, distribution, sample_size,
                               q95_B = 300L, q95_seed = NULL, q95_boot_indices = NULL) {
  true_mu <- S$true_mean
  C <- S$components_by_size[[as.character(sample_size)]]
  
  # Si no hay datos para este tamaño, devolvemos filas con NA de forma segura
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
    na_metrics <- as.list(c(mse=NA_real_, mse_q95=NA_real_, mse_max=NA_real_,
                            bias=NA_real_, abs_bias=NA_real_,
                            variance=NA_real_, mad=NA_real_, iqr=NA_real_))
    return(dplyr::bind_rows(
      dplyr::bind_cols(base_row, estimator="robust",      na_metrics),
      dplyr::bind_cols(base_row, estimator="mean",        na_metrics),
      dplyr::bind_cols(base_row, estimator="median",      na_metrics),
      dplyr::bind_cols(base_row, estimator="trimmed20",   na_metrics),
      dplyr::bind_cols(base_row, estimator="harmonic",    na_metrics),
      dplyr::bind_cols(base_row, estimator="geometric",   na_metrics),
      dplyr::bind_cols(base_row, estimator="mode_hsm",    na_metrics),
      dplyr::bind_cols(base_row, estimator="mode_parzen", na_metrics),
      dplyr::bind_cols(base_row, estimator="trimean",     na_metrics),
      dplyr::bind_cols(base_row, estimator="huber",       na_metrics),
      dplyr::bind_cols(base_row, estimator="biweight",    na_metrics)
    ))
  }
  
  # Selector de columna robusto: si no existe, devuelve NA
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
  
  dplyr::bind_rows(
    dplyr::bind_cols(base_row, estimator = "robust",      as.list(m_robust)),
    dplyr::bind_cols(base_row, estimator = "mean",        as.list(m_mean)),
    dplyr::bind_cols(base_row, estimator = "median",      as.list(m_median)),
    dplyr::bind_cols(base_row, estimator = "trimmed20",   as.list(m_trim20)),
    dplyr::bind_cols(base_row, estimator = "harmonic",    as.list(m_harm)),
    dplyr::bind_cols(base_row, estimator = "geometric",   as.list(m_geom)),
    dplyr::bind_cols(base_row, estimator = "mode_hsm",    as.list(m_mode_h)),
    dplyr::bind_cols(base_row, estimator = "mode_parzen", as.list(m_mode_p)),
    dplyr::bind_cols(base_row, estimator = "trimean",     as.list(m_trimean)),
    dplyr::bind_cols(base_row, estimator = "huber",       as.list(m_huber)),
    dplyr::bind_cols(base_row, estimator = "biweight",    as.list(m_biwt))
  )
}









# ---- Helpers para reducir repetición ------------------------------------

.ga_default_ctrl <- list(
  objective         = "q95",
  lambda_instab     = 0.15,
  bootstrap_B       = 200L,
  t_size            = 2L,
  elitism           = 2L,
  immigrant_rate    = 0.10,
  mutation_rate_init= 0.18,
  init_alpha        = 1.0,
  alpha_mut         = 1.0,
  check_every       = 5L,
  patience          = 3L,
  min_delta         = 0.005,
  mix_w_q95         = 0.7,
  mix_w_max         = 0.3,
  lambda_entropy    = 0.015,
  max_w_dominance   = 0.80,
  lambda_dominance  = 0.05,
  lambda_bias       = 0.05,
  minibatch_frac    = NULL,
  minibatch_min     = 12L
)

.ga_pick_mode <- function(sample_sizes, num_samples, pop_size, generations) {
  n_scen_full <- nrow(build_scenarios_full())
  work_score  <- n_scen_full * length(sample_sizes) * num_samples * pop_size * generations
  if (work_score >= 5e7) "light" else "full"
}

.ga_setup_cluster <- function(use_parallel, seed) {
  if (!use_parallel) return(NULL)
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::clusterEvalQ(cl, { runif(1); NULL })
  cl
}

.ga_export_cluster <- function(cl, env) {
  if (is.null(cl)) return(invisible(NULL))
  parallel::clusterExport(
    cl,
    varlist = c(
      "fitness_universal","prepped","N_EST","ESTIMATOR_NAMES",
      "q95_boot","huber_loss","scenario_weight",".normalize_simplex",
      ".ga_eval_split","ctrl"   # <-- AÑADIDOS
    ),
    envir = env
  )
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(dplyr); library(tibble); library(MASS) })
    NULL
  })
}


.ga_decay_mut <- function(g, gmax, m0, mmin = 0.05) {
  mmin + (m0 - mmin) * (log(gmax + 1) - log(g + 1)) / log(gmax + 1)
}

.ga_eval_split <- function(w, prepped, idx, ctrl, q95_seed) {
  fitness_universal(
    w, prepped, objective = ctrl$objective,
    idx = idx, use_boot = TRUE, B = ctrl$bootstrap_B,
    lambda_instab   = ctrl$lambda_instab,
    lambda_entropy  = ctrl$lambda_entropy,
    max_w_dominance = ctrl$max_w_dominance,
    lambda_dominance= ctrl$lambda_dominance,
    mix_w_q95       = ctrl$mix_w_q95,
    mix_w_max       = ctrl$mix_w_max,
    lambda_bias     = ctrl$lambda_bias,
    q95_seed        = q95_seed
  )
}

# =========================== GA (single split) ===========================

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
                                                  # ---- minibatch only for LIGHT scenarios ----
                                                  minibatch_frac = NULL,
                                                  minibatch_min  = 12L,
                                                  # ---- CRN support ----
                                                  crn_env = NULL,
                                                  fam_key = NULL) {
  set.seed(seed)
  
  # ------ Unificar controles (evita repetir en 20 sitios) ---------------
  ctrl <- modifyList(.ga_default_ctrl, list(
    objective          = match.arg(objective),
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
    mix_w_q95          = mix_w_q95,
    mix_w_max          = mix_w_max,
    lambda_entropy     = lambda_entropy,
    max_w_dominance    = max_w_dominance,
    lambda_dominance   = lambda_dominance,
    lambda_bias        = lambda_bias,
    minibatch_frac     = minibatch_frac,
    minibatch_min      = minibatch_min
  ))
  
  # ---- Elegir light/full y preparar escenarios --------------------------
  scenario_mode <- .ga_pick_mode(sample_sizes, num_samples, pop_size, generations)
  prepped <- prep_scenarios(
    dist_name, dist_param_grid,
    sample_sizes = sample_sizes,
    num_samples  = num_samples,
    scenario_mode = scenario_mode,
    seed         = seed,
    crn_env      = crn_env,
    fam_key      = if (is.null(fam_key)) dist_name else fam_key
  )
  prepped_mode <- attr(prepped, "scenario_mode")
  
  # ---- Partición train/val si no viene dada -----------------------------
  all_idx <- seq_along(prepped)
  if (is.null(val_idx) || is.null(train_idx)) {
    folds <- make_folds(length(all_idx), k = 3, seed = seed)
    val_idx   <- sort(unlist(folds[[3]]))
    train_idx <- sort(setdiff(all_idx, val_idx))
  }
  
  # ---- Warm-starts ------------------------------------------------------
  warm_cache <- get_warm_start(dist_name, seed)
  warm <- build_warm_starts()
  if (!is.null(warm_start_matrix)) warm <- unique(rbind(warm_start_matrix, warm))
  if (!is.null(warm_cache))        warm <- unique(rbind(warm_cache, warm))
  rest <- max(0, pop_size - nrow(warm))
  ga_population <- if (rest > 0) rbind(warm, init_population(rest, N_EST, alpha = ctrl$init_alpha)) else warm
  ga_population <- t(apply(ga_population, 1, .normalize_simplex))
  
  # ---- Paralelización (encapsulada) ------------------------------------
  cl <- .ga_setup_cluster(use_parallel, seed)
  on.exit({ if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  .ga_export_cluster(cl, environment())
  
  # ---- Tracking ---------------------------------------------------------
  conv <- if (record_convergence)
    data.frame(gen=integer(0), best_train=numeric(0), med_train=numeric(0),
               best_val=numeric(0),   med_val=numeric(0)) else NULL
  val_hist <- numeric(0); best_val_so_far <- Inf; no_improve <- 0L
  
  # ======================== Bucle evolutivo ==============================
  for (gen in seq_len(generations)) {
    mut_rate     <- .ga_decay_mut(gen, generations, ctrl$mutation_rate_init)
    q95_seed_gen <- 1e6 + seed*1000 + gen
    
    # ---- Minibatch LIGHT-only ------------------------------------------
    train_sub <- train_idx
    if (identical(prepped_mode, "light") && !is.null(ctrl$minibatch_frac)) {
      k <- max(ctrl$minibatch_min, ceiling(length(train_idx) * ctrl$minibatch_frac))
      k <- min(k, length(train_idx))
      train_sub <- .seed_scope(9e6 + seed*1000 + gen, base::sample(train_idx, k, replace = FALSE))
    }
    
    # ---- Fitness para toda la población --------------------------------
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
    
    
    # ---- Elitismo + selección + variación -------------------------------
    elites   <- ga_population[order(fitness_scores)[1:min(ctrl$elitism, nrow(ga_population))], , drop = FALSE]
    sel_idx  <- tournament_select(fitness_scores, nrow(ga_population), t_size = ctrl$t_size, replace = TRUE)
    selected <- ga_population[sel_idx, , drop = FALSE]
    
    n_off <- nrow(ga_population) - nrow(elites)
    offspring <- if (n_off > 0) t(replicate(n_off, {
      p1 <- selected[base::sample(nrow(selected), 1), ]
      p2 <- selected[base::sample(nrow(selected), 1), ]
      mutate_weights(crossover(p1, p2), mutation_rate = mut_rate, alpha_mut = ctrl$alpha_mut)
    })) else matrix(numeric(0), 0, ncol(ga_population))
    
    ga_population <- rbind(elites, offspring)
    
    # ---- Inmigración periódica -----------------------------------------
    if (gen %% 15L == 0 && ctrl$immigrant_rate > 0) {
      m <- max(1L, floor(nrow(ga_population) * ctrl$immigrant_rate))
      worst <- order(fitness_scores, decreasing = TRUE)[1:m]
      ga_population[worst, ] <- init_population(length(worst), N_EST, alpha = ctrl$init_alpha)
    }
    
    # ---- Validación y early-stopping -----------------------------------
    if (gen %% ctrl$check_every == 0 || gen == generations) {
      best_idx <- which.min(fitness_scores)
      w_best   <- ga_population[best_idx, ]
      
      eval_on <- function(w, idx)
        .ga_eval_split(w, prepped, idx, ctrl, q95_seed = q95_seed_gen)
      
      best_train <- eval_on(w_best, train_idx)   # reporte (train completo)
      best_val   <- eval_on(w_best,   val_idx)   # val completo (nunca minibatch)
      
      if (record_convergence) {
        med_train <- stats::median(fitness_scores)
        probe_val <- sapply(base::sample(seq_len(nrow(ga_population)),
                                         min(20, nrow(ga_population))),
                            function(i) eval_on(ga_population[i, ], val_idx))
        conv <- rbind(conv, data.frame(gen=gen, best_train=best_train, med_train=med_train,
                                       best_val=best_val,   med_val=stats::median(probe_val)))
      }
      
      val_hist <- c(val_hist, best_val)
      if (best_val < best_val_so_far - 1e-12) { best_val_so_far <- best_val; no_improve <- 0L }
      else no_improve <- no_improve + 1L
      if (no_improve > 0L) mut_rate <- min(0.40, mut_rate * (1 + 0.10 * no_improve))
      
      cat(sprintf("[%-10s] Gen %3d | Train %.6g | Val %.6g | mut=%.3f | mode=%s%s\n",
                  dist_name, gen, best_train, best_val, mut_rate, prepped_mode,
                  if (!is.null(ctrl$minibatch_frac) && identical(prepped_mode, "light"))
                    sprintf(" (mb=%d/%d)", length(train_sub), length(train_idx)) else ""))
      
      if (should_stop(val_hist, patience = ctrl$patience, min_delta = ctrl$min_delta)) {
        message(sprintf("[%s] Early stopping at gen %d.", dist_name, gen))
        break
      }
    }
  }
  
  # ======================== Selección final ==============================
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
  
  # warm-start bank
  keep_k <- min(nrow(ga_population), max(10L, as.integer(ctrl$elitism)))
  set_warm_start(dist_name, seed, ga_population[order(final_scores_val)[1:keep_k], , drop = FALSE])
  
  # ======================== Resultados por escenario =====================
  scenario_rows <- lapply(seq_along(prepped), function(i) {
    S <- prepped[[i]]
    dplyr::bind_rows(lapply(sample_sizes, function(n)
      summarize_scenario(S, as.numeric(best_weights), dist_name, n,
                         q95_B = ctrl$bootstrap_B, q95_seed = q95_seed_final)))
  })
  scenario_df <- dplyr::bind_rows(scenario_rows)
  
  agg_by_est <- scenario_df %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(mean_mse = base::mean(mse, na.rm = TRUE),
                     q95_mse  = base::mean(mse_q95, na.rm = TRUE),
                     max_mse  = base::max(mse_max, na.rm = TRUE),
                     mean_bias= base::mean(bias, na.rm = TRUE),
                     mean_variance = base::mean(variance, na.rm = TRUE),
                     .groups = "drop")
  robust_row <- dplyr::filter(agg_by_est, estimator == "robust")
  
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
  
  list(weights        = as.numeric(best_weights),
       estimator_str  = weights_to_formula(as.numeric(best_weights)),
       topk_weights   = ga_population[topk_idx, , drop = FALSE],
       scenario_table = scenario_df,
       overall        = overall,
       convergence    = conv,
       best_val_score = if (return_val_best) best_val_score else NA_real_)
}

  

# =========================== GA con K-fold CV (K=3) ======================

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
                                                     # --- loss mix and penalties ---
                                                     mix_w_q95 = 0.7,
                                                     mix_w_max = 0.3,
                                                     lambda_entropy = 0.015,
                                                     max_w_dominance = 0.80,
                                                     lambda_dominance = 0.05,
                                                     lambda_bias = 0.05,
                                                     # --- LIGHT-only minibatch knobs ---
                                                     minibatch_frac = NULL,
                                                     minibatch_min  = 12L,
                                                     # --- CRN / seeds ---
                                                     crn_env = NULL,
                                                     fam_key = NULL) {
  set.seed(seed)
  
  # --- Elegir scenario_mode con la MISMA heurística que usa el GA interno ---
  n_scen_full <- nrow(build_scenarios_full())
  work_score  <- n_scen_full * length(sample_sizes) * num_samples * pop_size * generations_per_fold
  scenario_mode <- if (work_score >= 5e7) "light" else "full"
  
  # --- Preparar una sola vez para conocer cuántos escenarios hay y crear folds ---
  prepped_ref <- prep_scenarios(
    dist_name, dist_param_grid,
    sample_sizes  = sample_sizes,
    num_samples   = num_samples,
    scenario_mode = scenario_mode,
    seed          = seed,
    crn_env       = crn_env,
    fam_key       = if (is.null(fam_key)) dist_name else fam_key
  )
  all_idx <- seq_along(prepped_ref)
  
  # --- Folds sobre el mismo # de escenarios que usará el GA interno ---
  folds <- make_folds(length(all_idx), k = k_folds, seed = seed)
  fold_results <- vector("list", k_folds)
  warm_bank <- NULL
  
  # Helper para evitar repetir el largo llamado
  run_fold <- function(train_idx, val_idx, seed_fold, warm_mat, gens, pop_sz,
                       retrain_all = FALSE, retrain_warm = NULL) {
    evolve_universal_estimator_per_family(
      dist_name       = dist_name,
      dist_param_grid = dist_param_grid,
      sample_sizes    = sample_sizes,
      num_samples     = num_samples,
      pop_size        = pop_sz,
      generations     = gens,
      seed            = seed_fold,
      objective       = objective,
      use_parallel    = use_parallel,
      train_idx       = train_idx,
      val_idx         = val_idx,
      lambda_instab   = lambda_instab,
      bootstrap_B     = bootstrap_B,
      t_size          = t_size,
      elitism         = elitism,
      immigrant_rate  = immigrant_rate,
      mutation_rate_init = mutation_rate_init,
      init_alpha      = init_alpha,
      alpha_mut       = alpha_mut,
      check_every     = check_every,
      patience        = patience,
      min_delta       = min_delta,
      warm_start_matrix = warm_mat,
      return_val_best = !retrain_all,
      mix_w_q95       = mix_w_q95,
      mix_w_max       = mix_w_max,
      lambda_entropy  = lambda_entropy,
      max_w_dominance = max_w_dominance,
      lambda_dominance= lambda_dominance,
      lambda_bias     = lambda_bias,
      minibatch_frac  = minibatch_frac,
      minibatch_min   = minibatch_min,
      crn_env         = crn_env,
      fam_key         = if (is.null(fam_key)) dist_name else fam_key
    )
  }
  
  # --- CV loop ---
  for (k in seq_len(k_folds)) {
    val_idx   <- sort(unlist(folds[[k]]))
    train_idx <- sort(setdiff(all_idx, val_idx))
    cat(sprintf("\n[%s] Fold %d/%d | Train=%d | Val=%d\n",
                dist_name, k, k_folds, length(train_idx), length(val_idx)))
    
    res_k <- run_fold(train_idx, val_idx, seed + k, warm_bank,
                      gens = generations_per_fold, pop_sz = pop_size)
    fold_results[[k]] <- res_k
    warm_bank <- rbind(warm_bank, res_k$weights)
  }
  
  # --- Seleccionar mejor fold ---
  val_scores <- sapply(fold_results, function(x) x$best_val_score)
  best_fold <- which.min(val_scores)
  best_weights_cv <- fold_results[[best_fold]]$weights
  
  # --- Retrain final opcional en TODOS los escenarios ---
  if (final_retrain) {
    cat(sprintf("\n[%s] Final retrain on all scenarios with warm-start\n", dist_name))
    res_final <- run_fold(
      train_idx = all_idx,
      val_idx   = all_idx,
      seed_fold = seed + 999,
      warm_mat  = rbind(warm_bank, best_weights_cv),
      gens      = max(40, generations_per_fold),
      pop_sz    = max(60, pop_size)
    )
    return(list(fold_results = fold_results, final = res_final))
  } else {
    return(list(fold_results = fold_results,
                final = fold_results[[best_fold]]))
  }
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
  fam_dir  <- file.path(root_out, family_tag)
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
  for (k in seq_along(fr_list)) {
    fr <- fr_list[[k]]
    .write_csv_safe(fr$overall,        file.path(folds_dir, sprintf("fold%02d_OVERALL.csv",    k)))
    .write_csv_safe(fr$scenario_table, file.path(folds_dir, sprintf("fold%02d_SCENARIOS.csv",  k)))
    .write_csv_safe(fr$convergence,    file.path(folds_dir, sprintf("fold%02d_CONVERGENCE.csv", k)))
  }
}


# ====================== LAUNCHER: TODO EN UNA TIRADA ======================

run_all_one_shot <- function(
    families_to_run = c("lognormal","weibull"),
    # grids de sensibilidad y tuning
    pop_sizes       = c(80, 100),
    mutation_rates  = c(0.12, 0.18, 0.24),
    init_alphas     = c(0.5, 1.0),
    alpha_mut_set   = c(0.5, 1.0),
    immigrant_rates = c(0.05, 0.10),
    t_sizes         = c(2L, 3L),
    elitism_set     = c(1L, 2L),
    seeds           = c(101, 202),
    # setup GA / escenarios
    sample_sizes = c(300, 500, 1000, 2000, 5000),
    num_samples  = 80,
    generations_per_fold = 60,
    k_folds = 3,
    final_retrain = TRUE,
    # métrica y robustez
    objective    = "q95",
    lambda_instab_default = 0.15,
    bootstrap_B  = 200L,
    check_every  = 5L,
    patience     = 3L,
    min_delta    = 0.005,
    # infra
    use_parallel = TRUE,
    mix_w_q95 = 0.7, mix_w_max = 0.3,
    out_root = ".",
    # --- LIGHT-only minibatching (random subset of scenarios per gen) ---
    minibatch_frac = 0.5,   # e.g., 0.25 to enable; NULL disables
    minibatch_min  = 8L
) {
  run_id  <- paste0("GA_UNIFIED_", RUN_TS())
  root_out <- file.path(out_root, run_id); mkdirp(root_out)
  
  meta <- list(
    run_id = run_id, families = families_to_run, k_folds = k_folds, objective = objective,
    seeds = paste(seeds, collapse = ","),
    pop_sizes = pop_sizes, generations_per_fold = generations_per_fold,
    final_retrain = final_retrain, bootstrap_B = bootstrap_B,
    lambda_instab_default = lambda_instab_default,
    mutation_rates = mutation_rates, init_alphas = init_alphas,
    alpha_mut_set = alpha_mut_set, immigrant_rates = immigrant_rates,
    t_sizes = t_sizes, elitism_set = elitism_set,
    minibatch_frac = minibatch_frac, minibatch_min = minibatch_min,
    # Dejar explícito en MANIFEST que el modo de escenarios es automático
    scenario_mode = "auto"   # (elige full/light según carga; minibatch solo en light)
  )
  write_manifest(root_out, meta)
  
  # ---- CRN global para toda la corrida ----
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
  
  master_overall <- list()
  master_index   <- list()
  
  # Argumentos comunes pasados al CV (menos repetición al llamar)
  common_args <- list(
    sample_sizes         = sample_sizes,
    num_samples          = num_samples,
    generations_per_fold = generations_per_fold,
    objective            = objective,
    use_parallel         = use_parallel,
    k_folds              = k_folds,
    bootstrap_B          = bootstrap_B,
    check_every          = check_every,
    patience             = patience,
    min_delta            = min_delta,
    final_retrain        = final_retrain,
    mix_w_q95            = mix_w_q95,
    mix_w_max            = mix_w_max,
    minibatch_frac       = minibatch_frac,
    minibatch_min        = minibatch_min,
    crn_env              = crn_env
  )
  
  for (dist in families_to_run) {
    message(sprintf("\n========== FAMILY: %s ==========\n", dist))
    grid <- param_grids[[dist]]
    
    for (i in seq_len(nrow(cfg_grid))) {
      cfg <- cfg_grid[i, ]
      config_tag <- sprintf("ps%d_mr%.2f_a0%.2f_am%.2f_im%.2f_t%d_e%d_seed%d",
                            cfg$pop_size, cfg$mutation_rate_init, cfg$init_alpha, cfg$alpha_mut,
                            cfg$immigrant_rate, cfg$t_size, cfg$elitism, cfg$seed)
      cat(sprintf("[RUN] %s | %s\n", dist, config_tag))
      
      res_cv <- do.call(
        evolve_universal_estimator_per_family_cv,
        c(list(
          dist_name          = dist,
          dist_param_grid    = grid,
          pop_size           = cfg$pop_size,
          seed               = cfg$seed,
          lambda_instab      = cfg$lambda_instab,
          t_size             = cfg$t_size,
          elitism            = cfg$elitism,
          immigrant_rate     = cfg$immigrant_rate,
          mutation_rate_init = cfg$mutation_rate_init,
          init_alpha         = cfg$init_alpha,
          alpha_mut          = cfg$alpha_mut
        ),
        common_args)
      )
      
      # Guardar artefactos
      save_family_config_outputs(root_out, dist, config_tag, res_cv)
      
      # Rollups
      tmp_overall <- res_cv$final$overall
      tmp_overall$config_tag <- config_tag
      tmp_overall$family     <- dist
      master_overall[[length(master_overall) + 1L]] <- tmp_overall
      
      master_index[[length(master_index) + 1L]] <- data.frame(
        family = dist, config_tag = config_tag,
        pop_size = cfg$pop_size, mutation_rate_init = cfg$mutation_rate_init,
        init_alpha = cfg$init_alpha, alpha_mut = cfg$alpha_mut,
        immigrant_rate = cfg$immigrant_rate, t_size = cfg$t_size, elitism = cfg$elitism,
        seed = cfg$seed, lambda_instab = cfg$lambda_instab,
        stringsAsFactors = FALSE
      )
    }
  }
  
  all_overall <- dplyr::bind_rows(master_overall)
  all_index   <- dplyr::bind_rows(master_index)
  readr::write_csv(all_overall, file.path(root_out, "ROLLUP__OVERALL_ALL.csv"))
  readr::write_csv(all_index,   file.path(root_out, "ROLLUP__CONFIG_INDEX.csv"))
  
  write_manifest(root_out, meta, extra = utils::head(all_overall, 10))
  message(sprintf("\nLISTO. Artefactos en: %s", root_out))
  invisible(list(out_dir = root_out, overall = all_overall, index = all_index))
}


# # =============================== FINAL RUN ===============================
# # Objetivo: ejecutar la corrida completa sobre TODAS las familias principales.
# Resultados se guardan en "./Final/GA_UNIFIED_<timestamp>"
out <- run_all_one_shot(
  families_to_run = c("normal", "lognormal", "weibull",
                      "invgauss", "exgaussian", "exwald"),
  pop_sizes       = c(80, 100),
  mutation_rates  = c(0.12, 0.18),
  init_alphas     = c(0.5, 1.0),
  alpha_mut_set   = c(0.5, 1.0),
  immigrant_rates = c(0.05, 0.10),
  t_sizes         = c(2, 3),
  elitism_set     = c(1, 2),
  seeds           = c(101, 202),
  sample_sizes    = c(300, 500, 1000, 2000, 5000),
  num_samples     = 100,
  generations_per_fold = 60,
  k_folds         = 3,
  final_retrain   = TRUE,
  objective       = "mixed",
  lambda_instab_default = 0.15,
  bootstrap_B     = 300L,
  check_every     = 5L,
  patience        = 3L,
  use_parallel    = TRUE,
  out_root        = "./Final"
)


# # ============================ QUICK TEST (Smoke test) ============================
# # Objetivo: validar que TODO el pipeline corre sin errores, con settings mínimos.
# out <- run_all_one_shot(
#   families_to_run = c("lognormal"),   # solo una familia
#   pop_sizes       = c(20),            # población pequeña
#   mutation_rates  = c(0.15),
#   init_alphas     = c(1.0),
#   alpha_mut_set   = c(1.0),
#   immigrant_rates = c(0.10),
#   t_sizes         = c(2),
#   elitism_set     = c(1),
#   seeds           = c(101),
#   sample_sizes    = c(200, 500),      # tamaños pequeños
#   num_samples     = 10,               # muy pocas réplicas
#   generations_per_fold = 3,           # solo 3 generaciones
#   k_folds         = 2,                # 2 folds para rapidez
#   final_retrain   = FALSE,
#   objective       = "q95",
#   lambda_instab_default = 0.0,
#   bootstrap_B     = 50L,
#   check_every     = 1L,
#   patience        = 2L,
#   use_parallel    = TRUE,
#   out_root        = "./TEST_LIGHT"
# )


# ====================== SECOND LEVEL (Confidence test) ==========================
# Objetivo: probar estabilidad con más réplicas y 2 familias.
# out <- run_all_one_shot(
#   families_to_run = c("lognormal", "weibull"),
#   pop_sizes       = c(100),
#   mutation_rates  = c(0.15),
#   init_alphas     = c(1.0),
#   alpha_mut_set   = c(1.0),
#   immigrant_rates = c(0.10),
#   t_sizes         = c(2),
#   elitism_set     = c(2),
#   seeds           = c(101, 202),
#   sample_sizes    = c(200, 500, 1000),
#   num_samples     = 40,                # más réplicas (menor varianza)
#   generations_per_fold = 15,
#   k_folds         = 3,
#   final_retrain   = FALSE,             # aún sin retrain largo
#   objective       = "mixed",
#   lambda_instab_default = 0.30,        # ligera penalización extra
#   bootstrap_B     = 100L,
#   check_every     = 2L,
#   patience        = 3L,
#   use_parallel    = TRUE,
#   out_root        = "./CONF_TEST_3"
# )

