# =============================================================================
# 00_utils_debug_io.R
#
# Purpose
# - Provide small, reusable utilities used across the whole experiment.
# - Centralize lightweight run tracking, debugging toggles, and safe I/O helpers.
# - Keep side-effects explicit (file writes, progress bars, console messages).
#
# Design note
# This module is sourced first so that later modules can rely on these helpers.
# =============================================================================
.MOD_STATUS <- if (exists(".MOD_STATUS", inherits = TRUE) && is.environment(.MOD_STATUS)) {
  # Shared module-status environment: reused across files to confirm load order and detect missing modules.
  # This avoids silent failures where downstream code runs without required helpers, improving auditability.
  .MOD_STATUS
} else {
  # Create a fresh status environment when running the project in a new R session (clean state, no parents).
  new.env(parent = emptyenv())
}

mark_module_done <- function(module_id, extra = NULL) {
  # Persist a “module completed” flag with timestamp and optional metadata, and print a standardized log line.
  # Used for reproducibility and debugging across multi-file pipelines and long runs (clear sourcing trace).
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  .MOD_STATUS[[module_id]] <- list(done = TRUE, time = ts, extra = extra)
  
  cat(sprintf("[MODULE DONE] %s | %s%s\n",
              ts, module_id,
              if (!is.null(extra)) paste0(" | ", extra) else ""))
  flush.console()
  
  invisible(TRUE)
}

is_module_done <- function(module_id) {
  # Query helper to confirm whether a specific module_id was already marked done; returns TRUE/FALSE safely.
  # This supports defensive “preflight checks” before expensive GA stages that require earlier definitions.
  x <- try(.MOD_STATUS[[module_id]], silent = TRUE)
  is.list(x) && isTRUE(x$done)
}


# =============================================================================
# Package availability
#
# We rely on a small set of CRAN packages for: robust statistics, random draws,
# and data I/O. This block ensures they are installed and available at runtime.
# =============================================================================

required_packages <- c(
  # Declares the external dependencies used across modules; keeping them centralized reduces hidden imports.
  # Installing missing packages here ensures a new machine/node can run the pipeline without manual setup.
  "modeest",   # modos (hsm, parzen)
  "statmod",   # rinvgauss
  "gtools",    # rdirichlet
  "parallel",  # PSOCK cluster
  "MASS",      # rlm (Huber/Bisquare)
  "dplyr", "tibble", "readr"  # utilidades de datos
)

for (pkg in required_packages) {
  # Ensures each required package is available in the library. If missing, installs it with dependencies.
  # This is a deliberate side-effect: it trades strict reproducibility for practical “it runs anywhere” setup.
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}


# =============================================================================
# Debug logging helpers (opt-in)
#
# These helpers provide trace messages when diagnosing unexpected behavior.
# They are disabled by default and can be enabled via environment variables.
# =============================================================================
# Enable debugging with an environment variable, e.g.:
#   Sys.setenv(DEBUG_TOPK="1")
# On Windows Command Prompt, the equivalent is:
#   set DEBUG_TOPK=1
.debug_on <- function() identical(Sys.getenv("DEBUG_TOPK", "0"), "1")

.dbg <- function(...) {
  # Lightweight debug logger: prints timestamped messages only when DEBUG_TOPK=1.
  # Keeps normal runs quiet, but makes it easy to trace structural issues in result objects during diagnosis.
  if (!.debug_on()) return(invisible(NULL))
  msg <- paste0(..., collapse = "")
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[DEBUG %s] %s", ts, msg))
  invisible(NULL)
}

#
# Top-K extraction diagnostics
#
# Some runs may return results in slightly different shapes depending on
# the execution path (early-stopping, error handling, etc.). This helper
# prints the structure of a result object to explain why a Top-K matrix
# could not be recovered.
.debug_topk_diagnose <- function(res, label = "") {
  # Structural introspection for GA/CV outputs: prints key fields and shapes to explain Top-K extraction.
  # Useful when wrappers differ (CV vs single-run) or when early-stopping/errors change object structure.
  if (!.debug_on()) return(invisible(NULL))
  if (is.null(res)) {
    .dbg(label, " res is NULL (likely error caught by tryCatch).")
    return(invisible(NULL))
  }
  if (!is.list(res)) {
    .dbg(label, " res is not a list; class=", paste(class(res), collapse=","), ".")
    return(invisible(NULL))
  }
  nms <- names(res); if (is.null(nms)) nms <- character()
  .dbg(label, " names(res)=", paste(nms, collapse=", "))
  has_final <- !is.null(res$final) && is.list(res$final)
  .dbg(label, " has_final=", has_final)
  if (has_final) {
    fnms <- names(res$final); if (is.null(fnms)) fnms <- character()
    .dbg(label, " names(res$final)=", paste(fnms, collapse=", "))
    has_topk <- !is.null(res$final$topk_weights)
    has_w    <- !is.null(res$final$weights)
    .dbg(label, " final$topk_weights present=", has_topk, " | final$weights present=", has_w)
    if (!is.null(res$final$overall)) {
      ov <- res$final$overall
      if (is.data.frame(ov)) {
        wcols <- grep("^w_", colnames(ov), value = TRUE)
        .dbg(label, " overall rows=", nrow(ov), " cols=", ncol(ov),
             " | w_ cols=", paste(wcols, collapse=", "))
      } else {
        .dbg(label, " final$overall is not a data.frame; class=", paste(class(ov), collapse=","))
      }
    } else {
      .dbg(label, " final$overall is NULL")
    }
  }
  invisible(NULL)
}

# Capture a compact call stack string (useful inside error handlers).
.debug_stack <- function() {
  # Captures the current call chain as a single string; helpful inside tryCatch error handlers when you want
  # context without printing full tracebacks. This is diagnostics-only and should not affect RNG or results.
  calls <- sys.calls()
  if (length(calls) == 0) return("")
  paste(vapply(calls, function(x) paste(deparse(x), collapse=" "), character(1)), collapse="  <-  ")
}

for (pkg in required_packages) {
  # Attaches the required packages quietly (suppresses startup messages) so functions are available by name.
  # Centralizing library() calls here avoids repeated imports in every module and keeps logs cleaner.
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Small helper used later: x %||% y returns x if not NULL else y
# Convenience operator used widely for defaulting optional values (metadata fields, nested list elements).
# This reduces verbose if/else checks and makes result-handling code easier to read and less error-prone.
`%||%` <- function(x, y) if (!is.null(x)) x else y


# --- Extract Top-K candidate weights from a result object (robust across wrappers) ---
# Returns: matrix [k x p] or NULL
.extract_topk_from_res <- function(res_cv, K = 5L) {
  # Robust Top-K extractor: supports both CV wrapper outputs (list(final=...)) and direct single-run outputs.
  # It prioritizes explicit topk_weights, falls back to a single best weights vector, then tries w_* columns in
  # an 'overall' data.frame. Returns NULL if nothing consistent is found (caller may trigger diagnostics).
  if (is.null(res_cv)) return(NULL)
  K <- as.integer(K)
  if (!is.finite(K) || K < 1L) K <- 1L
  
  # res_cv can be a CV wrapper or a direct result. Prefer $final when present.
  final <- res_cv$final %||% res_cv
  
  # 1) Explicit Top-K
  if (!is.null(final$topk_weights)) {
    w <- as.matrix(final$topk_weights)
    if (nrow(w) > K) w <- w[seq_len(K), , drop = FALSE]
    return(w)
  }
  
  # 2) Single best vector -> Top-1
  if (!is.null(final$weights)) {
    wv <- as.numeric(final$weights)
    return(matrix(wv, nrow = 1L))
  }
  
  # 3) Fallback: if overall has w_* columns
  if (!is.null(final$overall) && is.data.frame(final$overall)) {
    wcols <- grep("^w_", names(final$overall), value = TRUE)
    if (length(wcols) > 0) {
      wv <- as.numeric(final$overall[1, wcols, drop = TRUE])
      return(matrix(wv, nrow = 1L))
    }
  }
  
  return(NULL)
}




# =============================================================================
# File-system and logging utilities
#
# These helpers standardize how we create directories, write CSV files safely,
# and print consistent timestamped messages for experiment logs.
# =============================================================================

RUN_TS <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

mkdirp <- function(p) {
  # “mkdir -p” equivalent: creates directories recursively if missing, then returns a normalized path.
  # Keeping directory creation centralized reduces scattered side-effects and prevents write failures mid-run.
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(normalizePath(p, winslash = "/", mustWork = FALSE))
}

# Simple timestamped logger used throughout the pipeline for reproducible logs.
catf <- function(fmt, ...) {
  # Standard console logger: prepends a human-readable timestamp and formats messages with sprintf().
  # Used across modules so logs remain consistent and searchable (especially when debugging long experiments).
  msg <- sprintf(fmt, ...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# --- Safe scalar extractor (top-level) -----------------------------------
# Used when building small tibbles/data.frames from possibly length>1 vectors.
.safe_scalar1 <- function(x, default = NA_real_) {
  # Extracts a single numeric value defensively from potentially empty/non-numeric inputs.
  # Prevents accidental recycling or list-columns when building summaries; returns default on invalid values.
  if (is.null(x) || length(x) < 1L) return(default)
  v <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(v)) return(default)
  v
}

# Write text lines to disk (creating parent folders if needed). Returns the written path.
safe_write_lines <- function(lines, path) {
  # Writes a character vector to a UTF-8 text file, ensuring the parent directory exists first.
  # This is used for MANIFEST/SESSION logs and other run artifacts where exact text output matters.
  mkdirp(dirname(path))
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  invisible(path)
}

# Robust saving helpers for RDS and CSV outputs.
safe_save_rds <- function(obj, path) {
  # Saves arbitrary R objects as RDS after creating the parent directory.
  # Keeping RDS writes centralized helps ensure consistent file layout across runs and nodes.
  mkdirp(dirname(path))
  saveRDS(obj, file = path)
  invisible(path)
}

safe_write_csv <- function(df, path) {
  # Writes a data.frame to CSV defensively: skips empty/invalid frames, creates parent directories, and uses
  # readr::write_csv when available for consistent formatting. Falls back to base write.csv for portability.
  if (!is.data.frame(df) || nrow(df) == 0L) return(invisible(NULL))
  mkdirp(dirname(path))
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(df, path)
  } else {
    utils::write.csv(df, path, row.names = FALSE)
  }
  invisible(path)
}



# =============================================================================
# Experiment index (run-level ledger)
#
# The experiment index is a single CSV that records each configuration evaluated
# (family, sample size, GA settings, stage, seed, etc.) and where its artifacts were
# written. This makes it easy to audit a run and to reproduce a specific row later.
# =============================================================================

.index_rows <- new.env(parent = emptyenv())
index_init <- function() .index_rows$rows <- list()

index_add <- function(family, stage, output_csv, n_candidates = NA_integer_,
                      n_scenarios = NA_integer_, seed = NA_integer_) {
  # Appends one row to an in-memory ledger of outputs produced during this run. Each row links a logical
  # “stage” (e.g., HPF1/HALVING/FINAL) to its artifact path and key counts, enabling later audit/reproduction.
  .index_rows$rows[[length(.index_rows$rows)+1]] <<- data.frame(
    family = family,
    stage = stage,
    output_csv = output_csv,
    n_candidates = n_candidates,
    n_scenarios = n_scenarios,
    seed = seed,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
}

index_flush <- function(run_dir) {
  # Materializes the in-memory index ledger as experiment_index.csv under run_dir. This file acts as the
  # run’s “table of contents” so collaborators can locate artifacts and reproduce a row without guessing paths.
  df <- dplyr::bind_rows(.index_rows$rows)
  safe_write_csv(df, file.path(run_dir, "experiment_index.csv"))
  df
}

# Attempt to read git metadata (short commit hash and branch) when available.
.git_info <- function() {
  # Best-effort capture of git state for reproducibility: commit hash and branch name. If the code is not
  # in a git repo (or git is unavailable), returns NA rather than failing the run (non-critical metadata).
  tryCatch({
    hash <- suppressWarnings(system("git rev-parse --short HEAD", intern = TRUE))
    br   <- suppressWarnings(system("git rev-parse --abbrev-ref HEAD", intern = TRUE))
    if (length(hash) && length(br)) sprintf("%s (%s)", hash[1], br[1]) else NA_character_
  }, error = function(e) NA_character_)
}

# Write MANIFEST.txt and a session snapshot (SESSION.txt).
# These files capture run metadata and the R environment for reproducibility.
write_manifest <- function(dir_out, meta, extra = NULL, write_session = TRUE) {
  # Creates human-readable run metadata files:
  # - MANIFEST.txt captures configuration (families, seeds, GA settings, scenario mode, minibatching, etc.)
  # - SESSION.txt captures sessionInfo() and an installed.packages snapshot. Together they support auditing
  #   and reproduction of a run on a different machine or at a later date with minimal ambiguity.
  mkdirp(dir_out)
  fmt_vec <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.list(x))  return(paste(capture.output(str(x, give.attr = FALSE)), collapse = " "))
    if (length(x) > 1) paste0("[", paste(x, collapse = ", "), "]") else as.character(x)
  }
  
  scenario_mode_line <- meta$scenario_mode %||% "full"
  git_line <- .git_info()
  
  lines <- c(
    sprintf("# RUN_ID: %s", meta$run_id %||% RUN_TS()),
    sprintf("timestamp: %s", as.character(Sys.time())),
    sprintf("git: %s", ifelse(is.na(git_line), "NA", git_line)),
    sprintf("platform: %s | R: %s | locale: %s",
            paste(R.version$platform, R.version$arch, sep = "/"),
            getRversion(), Sys.getlocale()),
    "",
    "## CONFIG",
    sprintf("families: %s", fmt_vec(meta$families)),
    sprintf("k_folds: %s", meta$k_folds %||% NA),
    sprintf("objective: %s", meta$objective %||% NA),
    sprintf("seeds: %s", fmt_vec(meta$seeds %||% meta$seed)),
    sprintf("pop_sizes: %s", fmt_vec(meta$pop_sizes)),
    sprintf("generations_per_fold: %s", meta$generations_per_fold %||% NA),
    sprintf("final_retrain: %s", meta$final_retrain %||% NA),
    sprintf("bootstrap_B: %s", meta$bootstrap_B %||% NA),
    sprintf("lambda_instab (default): %s", meta$lambda_instab_default %||% NA),
    sprintf("scenario_mode: %s", scenario_mode_line),
    sprintf("minibatch_frac: %s", meta$minibatch_frac %||% NA),
    sprintf("minibatch_min: %s", meta$minibatch_min %||% NA)
  )
  
  if (!is.null(extra)) {
    # Optional preview block for additional metadata objects (e.g., scenario subset tables). We only print
    # the head() to keep MANIFEST readable while still exposing enough context for quick sanity checks.
    lines <- c(lines, "", "## EXTRA (preview)", capture.output(print(utils::head(extra, 10))))
  }
  
  mf_path <- file.path(dir_out, "MANIFEST.txt")
  safe_write_lines(lines, mf_path)
  catf("MANIFEST escrito en: %s", mf_path)
  
  if (isTRUE(write_session)) {
    # SESSION.txt is produced via sink() so we capture rich printed output exactly as R would display it.
    # This helps reproducing environments across nodes, especially when package versions differ subtly.
    sess_path <- file.path(dir_out, "SESSION.txt")
    sink(sess_path)
    on.exit({ if (sink.number() > 0) sink() }, add = TRUE)
    cat("# sessionInfo()\n\n"); print(utils::sessionInfo()); cat("\n")
    pkgs <- tryCatch(installed.packages()[, c("Package","Version")], error = function(e) NULL)
    if (!is.null(pkgs)) {
      cat("# installed.packages() snapshot (Package / Version)\n\n")
      print(utils::head(as.data.frame(pkgs), 50))
    }
    sink()
    catf("SESSION escrito en: %s", sess_path)
  }
  
  invisible(mf_path)
}

# =============================================================================
# Scenario universe export
#
# We export the full scenario grids ("full" and "light") to CSV so collaborators can
# see exactly which contamination settings were considered. Stage-specific subsets
# used by HPF1/HPF2/HALVING are exported elsewhere.
# =============================================================================

add_scenario_ids <- function(df) {
  # Adds a stable sequential scenario_id to a scenario grid so subsets can be referenced unambiguously.
  # This enables exporting subsets and re-loading them later without relying on row order or implicit indices.
  df <- as.data.frame(df)
  df$scenario_id <- sprintf("SC%05d", seq_len(nrow(df)))
  df
}

write_scenarios_universe <- function(dir_out, scenario_mode = c("full","light")) {
  # Exports the complete scenario grid (full or light) as a CSV artifact under the run directory. This is a
  # transparency/reproducibility step: collaborators can inspect the contamination design space used by runs.
  scenario_mode <- match.arg(scenario_mode)
  sc <- build_scenarios(scenario_mode)
  sc <- add_scenario_ids(sc)
  
  safe_write_csv(sc, file.path(dir_out, sprintf("scenarios_universe_%s.csv", scenario_mode)))
  invisible(sc)
}




# =============================================================================
# Progress reporting
#
# A lightweight progress bar for long runs. This is purely cosmetic and does not
# affect the experiment's stochastic behavior.
# =============================================================================
.progress <- new.env(parent = emptyenv())

progress_init <- function(total_steps) {
  # Initializes a simple text progress bar where progress is tracked as percentage of total_steps.
  # Stored in an environment to avoid global variables and to keep progress state encapsulated and mutable.
  .progress$done  <- 0L
  .progress$total <- max(1L, as.integer(total_steps))
  .progress$pb    <- utils::txtProgressBar(min = 0, max = 100, style = 3)
  utils::setTxtProgressBar(.progress$pb, 0)
}

progress_step <- function(label = NULL, inc = 1L) {
  # Advances progress by 'inc' steps, updates the progress bar, and optionally logs a labeled checkpoint.
  # This is intentionally separated from computation logic so it never changes RNG state or experiment output.
  .progress$done <- min(.progress$done + as.integer(inc), .progress$total)
  pct <- round(100 * .progress$done / .progress$total)
  if (!is.null(.progress$pb)) utils::setTxtProgressBar(.progress$pb, pct)
  if (!is.null(label)) catf("Progreso %3d%% | %s", pct, label)
}

progress_finalize <- function() {
  # Closes the progress bar cleanly and prints a final “DONE” log line. Ensures the console is left in a clean
  # state (no partially drawn progress bar) which is important when logs are being captured to files on nodes.
  if (!is.null(.progress$pb)) {
    utils::setTxtProgressBar(.progress$pb, 100)
    close(.progress$pb); .progress$pb <- NULL
  }
  catf("Progreso 100%% | DONE")
}

# Mark this module as successfully sourced (used for sanity checks in the console).
mark_module_done("00_utils_debug_io.R")
