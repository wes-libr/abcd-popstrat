##############################################################
##  00_utils.R — Shared utilities for the music-ed pipeline ##
##  ABCD Music Education & Academic Achievement              ##
##  Wes Thompson (wes.stat@gmail.com)                        ##
##############################################################

# ── Required packages ──────────────────────────────────────

required_pkgs <- c(
 "NBDCtools",  # ABCD data loading (also loads NBDCtoolsData)
 "dplyr",      # data wrangling
 "tidyr",      # reshaping
 "stringr",    # string helpers
 "grf",        # causal forests
 "nlme",       # mixed-effects models
 "mgcv",       # GAMMs
 "lavaan",     # SEM / CLPM / RI-CLPM / LGC
 "fixest",     # fast fixed-effects estimation
 "ggplot2",    # plotting
 "patchwork"   # combining plots
)

# Additional packages used by specific modules (loaded conditionally)
optional_pkgs <- c(
 "OpenMx"      # twin models / DoC / MRDoC2 (12, 14)
)

# ── SAIQ ordinal → midpoint mappings ──────────────────────
# Coding verified against ABCD 7.0 data dictionary (ph_p_saiq.json,
# ph_y_saiq.json). The non-trivial ones:
#
# Parent instrument _003 (minutes per session, ORDINAL):
#   0=0, 1=<30 (mid 15), 2=30, 3=45, 4=60, 5=90, 6=120,
#   7=150, 8=180, 9=>180 (mid 210)
saiq_minutes_map <- c(
 `0` = 0,   `1` = 15,  `2` = 30,  `3` = 45,
 `4` = 60,  `5` = 90,  `6` = 120, `7` = 150,
 `8` = 180, `9` = 210
)
map_saiq_minutes <- function(x) saiq_minutes_map[as.character(x)]

# Parent instrument _002 (days per week, ORDINAL):
#   0-7 = direct; 8 = "Once every 2 weeks" (0.5 days/wk);
#   9 = "One day every month" (1/4.345 days/wk);
#   10 = "Less than one day per month" (conservative 0.1 days/wk)
# IMPORTANT: the previous pipeline treated 8/9/10 as literal 8/9/10 days/
# week, which grossly overestimates activity hours for those subjects.
saiq_dayswk_map <- c(
 `0` = 0, `1` = 1, `2` = 2, `3` = 3, `4` = 4,
 `5` = 5, `6` = 6, `7` = 7,
 `8` = 0.5, `9` = 1/4.345, `10` = 0.1
)
map_saiq_dayswk <- function(x) saiq_dayswk_map[as.character(x)]

# Parent instrument _001__l (months since last visit):
# baseline _001 has max 12. Longitudinal _001__l adds code 13 = ">12 months".
# Cap 13 at 12 (conservative; typical ABCD visits are 12-14 months apart).
map_saiq_months_long <- function(x) {
 x <- saiq_num(x)
 ifelse(is.na(x), NA_real_, pmin(x, 12))
}

# Bulletproof numeric casting for SAIQ raw inputs.
# NBDCtools returns baseline SAIQ _001/_002/_003 as R factors with levels
# "0","1",...,"999". Factors cannot be used directly in arithmetic (R silently
# returns NA for the entire product via factor->integer->NA coercion), which
# historically caused all baseline activity hours to collapse to NA. This
# helper casts to numeric and sets the "decline to answer" code (999) to NA.
# Placed below the map_* functions so they can call it after definition.
saiq_num <- function(x) {
 v <- suppressWarnings(as.numeric(as.character(x)))
 v[!is.na(v) & v == 999] <- NA_real_
 v
}

# Youth yr3/4 lmusic _002 (hours per DAY listening to music, ORDINAL):
#   0 = <1 hr (mid 0.5); 1-7 = [k, k+1) (mid k+0.5); 8 = >8 hr (mid 10)
saiq_y_lmusic_hrday_map <- c(
 `0` = 0.5, `1` = 1.5, `2` = 2.5, `3` = 3.5, `4` = 4.5,
 `5` = 5.5, `6` = 6.5, `7` = 7.5, `8` = 10
)
map_y_lmusic_hrday <- function(x) saiq_y_lmusic_hrday_map[as.character(x)]

# Youth yr3/4 read _002 (hours per WEEK reading for pleasure, ORDINAL
# with uneven bins):
#   0=<1 (0.5), 1=[1,2) (1.5), 2=[2,3) (2.5), 3=[3,4) (3.5), 4=[4,5) (4.5),
#   5=[5,6) (5.5), 6=[6,8) (7), 7=[8,12) (10), 8=[12,18) (15),
#   9=[18,24) (21), 10=[24,36) (30), 11=[36,72) (54), 12=>72 (90)
saiq_y_read_hrwk_map <- c(
 `0` = 0.5,  `1` = 1.5,  `2` = 2.5,  `3` = 3.5,  `4` = 4.5,
 `5` = 5.5,  `6` = 7,    `7` = 10,   `8` = 15,   `9` = 21,
 `10` = 30,  `11` = 54,  `12` = 90
)
map_y_read_hrwk <- function(x) saiq_y_read_hrwk_map[as.character(x)]

load_packages <- function(pkgs = required_pkgs) {
 for (pkg in pkgs) {
   if (!requireNamespace(pkg, quietly = TRUE)) {
     message("Package '", pkg, "' not installed — skipping.")
   } else {
     suppressPackageStartupMessages(library(pkg, character.only = TRUE))
   }
 }
}


# ── Paths ──────────────────────────────────────────────────
#
# DESIGN NOTE (blinded analyses):
# This repo (abcd-analyses/) is shared with Claude / Cowork. Claude
# is blinded to participant-level data: it writes code, the user runs it
# locally, and only AGGREGATE outputs return to the shared tree. Therefore
# every individual-level intermediate (the analytic panel and all *_results
# model objects, whose $model frames / residuals / random-effect posteriors
# carry per-participant rows) is written to `dir_derived` — a host-only path
# OUTSIDE this repo. The shared tree gets figures + suppressed CSV tables only.
#
# `output` resolves relative to the pipeline working directory (analyses/R/);
# run_all.R sets that working directory before sourcing.

## ---- dir_derived (host-only, outside the shared repo) -----------------------
## Hard-coded absolute literal so a single source of truth can be parsed by
## .Rprofile (R_HISTFILE redirect) and scripts/preflight.sh with the same
## regex: ^dir_derived\s*<-\s*"...". It is a sibling of the raw release, on
## the user's machine, and NOT inside abcd-analyses/. Adjust to your
## host if running elsewhere — but keep it outside the shared repo.
## Set ABCD_DATA_ROOT to the folder holding your NBDC download (never inside this repo).
abcd_root <- Sys.getenv("ABCD_DATA_ROOT", "~/nbdc-datasets")
dir_derived <- file.path(abcd_root, "abcd_7-0_phenotype/7_0/derived")

## ---- Project root (the shared repo) -----------------------------------------
## Resolved from the pipeline working directory (analyses/R/). run_all.R and
## the project setup instruct the user to setwd() there before sourcing.
.pipeline_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
dir_project   <- normalizePath(file.path(.pipeline_dir, "..", ".."),
                               winslash = "/", mustWork = FALSE)
## Sanity fallback: if that doesn't look like the repo root (no analyses/R/
## beneath it), walk upward from the working directory to find it.
if (!dir.exists(file.path(dir_project, "analyses", "R"))) {
  .cand <- .pipeline_dir
  for (.i in 1:4) {
    if (dir.exists(file.path(.cand, "analyses", "R"))) { dir_project <- .cand; break }
    .cand <- normalizePath(file.path(.cand, ".."), winslash = "/", mustWork = FALSE)
  }
}

paths <- list(
 # Local ABCD 7.0 tabulated data directory (raw release).
 abcd_data = file.path(abcd_root, "abcd_7-0_phenotype/7_0/tabulated/release_abcd"),

 # Genetic ancestry admixture file (ancestry.tsv).
 ancestry  = file.path(abcd_root, "abcd_7-0_phenotype/ancestry.tsv"),

 # Polygenic scores, one file per score. Named vector maps the analytic column
 # name -> file. Two formats are read (read_plink_profile): PLINK 1 --score
 # ".profile" (FID IID PHENO CNT CNT2 SCORE) and PLINK 2 ".sscore"
 # (#FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG).
 #   pgs_dep = DEPRESSION PGS (common_pgs_*)          negative-control instrument (13, 14)
 #   pgs_soc = SOCIAL-CONNECTEDNESS PGS (soc_pgs_*)   negative-control instrument (13)
 #   pgs_ea  = EDUCATIONAL ATTAINMENT PGS, SSGAC 2022 GWAS (Okbay et al. 2022),
 #             PRS-CS weights, PLINK 2 --score on the ABCD MAF>0.01 genotypes
 #             (11,389 samples). Added 2026-10-19. Used only in the MR-style
 #             comparison of 13_pgs_iv.R (Method 11b).
 # PLINK files key participants as NDAR_INV...; ABCD 7.0 uses sub-XXXXXXXX.
 # Joins therefore go through canon_id() (below), never on raw id strings.
 # Set any entry to NA, or pgs to NA, to skip -> PGS-IV (13) / MRDoC2 (14) skip cleanly.
 pgs       = c(
   pgs_dep = file.path(abcd_root, "abcd_7-0_phenotype/common_pgs_abcd.profile"),
   pgs_soc = file.path(abcd_root, "abcd_7-0_phenotype/soc_pgs_abcd.profile"),
   pgs_ea  = file.path(abcd_root, "abcd-prs/plink_SSGAC_Educational_Attainment_2022.txt.sscore")
 ),

 # Pipeline intermediate data (RDS cache) — HOST-ONLY, outside the shared repo.
 # Holds the analytic panel + every model-fit object (individual-level).
 cache     = dir_derived,

 # Output figures and SUPPRESSED tables — inside the shared repo (aggregate only)
 output    = "output"
)

## ---- Blinding safety check --------------------------------------------------
## Refuse to start if the RDS cache resolves to anything inside the shared
## repo, because writing individual-level intermediates there would leak data
## into the Cowork-connected folder. Compares canonicalized paths.
.assert_outside_project <- function(path, project) {
  np_path    <- sub("/+$", "", normalizePath(path,    winslash = "/", mustWork = FALSE))
  np_project <- sub("/+$", "", normalizePath(project, winslash = "/", mustWork = FALSE))
  if (identical(np_path, np_project) ||
      startsWith(paste0(np_path, "/"), paste0(np_project, "/"))) {
    stop(
      "[00_utils] BLINDING VIOLATION: cache (dir_derived) resolves inside the shared repo.\n",
      "  cache       = ", np_path, "\n",
      "  dir_project = ", np_project, "\n",
      "Set dir_derived to a path OUTSIDE abcd-analyses/ so individual-level ",
      "RDS intermediates do not become readable by Claude/Cowork.",
      call. = FALSE)
  }
  invisible(TRUE)
}
.assert_outside_project(paths$cache, dir_project)

## ---- Blinding safety check: stray individual-level artifacts in the repo ----
## Broadened 2026-09-10. Previously this only looked at the legacy in-repo cache
## analyses/R/data/. Two escapes were possible and one of them actually happened:
##   (a) a save_rds()/saveRDS() to any OTHER in-repo path, not just data/;
##   (b) an R workspace image (.RData), written at quit, holding the whole
##       session — dat_long / dat_cumul / every fitted model object. A 267 MB
##       one was found in analyses/R/ on 2026-09-10.
## Scope: the analyses/ subtree (where the pipeline writes) plus the repo root's
## top level. Deliberately NOT the whole repo — grant_proposals/ is a large
## historical archive and scanning it on every startup is slow; `preflight.sh`
## sweeps the full tree instead. The synthetic fixtures are allowlisted: they
## are parametric draws, not real rows (see tests/fixtures/generate_fixtures.R).
.analyses_dir  <- normalizePath(file.path(.pipeline_dir, ".."), winslash = "/", mustWork = FALSE)
.fixtures_dir  <- normalizePath(file.path(.pipeline_dir, "tests", "fixtures", "data"),
                                winslash = "/", mustWork = FALSE)

.scan_repo <- function(pattern) {
  hits <- character(0)
  if (dir.exists(.analyses_dir)) {
    hits <- c(hits, list.files(.analyses_dir, pattern = pattern, recursive = TRUE,
                               all.files = TRUE, full.names = TRUE, no.. = TRUE))
  }
  if (dir.exists(dir_project)) {
    hits <- c(hits, list.files(dir_project, pattern = pattern, recursive = FALSE,
                               all.files = TRUE, full.names = TRUE, no.. = TRUE))
  }
  unique(normalizePath(hits, winslash = "/", mustWork = FALSE))
}

## (a) stray .rds anywhere under analyses/ except the synthetic fixtures
.stray_rds <- Filter(function(f) !startsWith(f, paste0(.fixtures_dir, "/")),
                     .scan_repo("\\.rds$"))
if (length(.stray_rds)) {
  stop("[00_utils] BLINDING VIOLATION: RDS file(s) inside the shared repo.\n",
       paste0("  ", .stray_rds, collapse = "\n"), "\n",
       "Every individual-level intermediate belongs in dir_derived (outside the ",
       "repo):\n  ", dir_derived, "\n",
       "Move or delete these before running. Only tests/fixtures/data/*.rds ",
       "(fully synthetic) may live in the repo.",
       call. = FALSE)
}

## (b) R workspace image. run_all.R sweeps these before we get here, so reaching
## this means 00_utils.R was sourced directly — still a violation, so stop.
.stray_rdata <- .scan_repo("^\\.RData$")
if (length(.stray_rdata)) {
  stop("[00_utils] BLINDING VIOLATION: R workspace image(s) in the shared repo.\n",
       paste0("  ", .stray_rdata, collapse = "\n"), "\n",
       "A .RData saved after a pipeline run contains the entire session, including ",
       "per-participant rows. Delete it (it regenerates by re-running), then quit R ",
       "with q(save = \"no\") — analyses/R/.Rprofile now defaults q() that way and ",
       "purges a stray .RData at startup.",
       call. = FALSE)
}

# Create the host-side cache and the in-repo output dir. Other `paths` entries
# are input files, not directories.
ensure_dirs <- function() {
 for (d in c(paths$cache, paths$output)) {
   if (!dir.exists(d)) dir.create(d, recursive = TRUE)
 }
}

# ── Disclosure control: small-cell suppression ─────────────
# Project-wide minimum cell size for any aggregate written to paths$output.
# Mirrors the HBCD pipeline's safe_write_csv(): rows whose count column(s)
# fall below the threshold have their value + count columns blanked, and a
# sidecar <name>.suppressed.csv logs how many rows were redacted.
min_cell_size      <- 10L   # ABCD data use terms: no cells below 10
suppression_marker <- "<10"

#' Write an aggregate data frame to CSV with small-cell suppression.
#'
#' @param df          data frame / tibble to write
#' @param path        output CSV path (under paths$output)
#' @param count_cols  REQUIRED. Columns holding raw counts; a row is suppressed
#'                    if ANY non-NA count column is below `min_cell`. Pass
#'                    character(0) if the table genuinely has no count columns.
#' @param value_cols  columns to redact on suppression (default: all non-key,
#'                    non-count columns)
#' @param key_cols    columns always kept visible (labels: method, exposure,
#'                    outcome, variable, level, wave, term)
#' Writes <path> (redacted) and <path>.suppressed.csv (one-row summary).
safe_write_csv <- function(df, path, count_cols,
                           value_cols = NULL,
                           key_cols   = c("method", "exposure", "outcome",
                                          "variable", "level", "wave", "term"),
                           min_cell   = min_cell_size,
                           marker     = suppression_marker) {
  if (!inherits(df, "data.frame")) stop("safe_write_csv() expects a data frame")
  if (missing(count_cols)) {
    stop("safe_write_csv(): `count_cols` is required. Pass an explicit ",
         "character vector (or character(0) if the table has no count columns).",
         call. = FALSE)
  }
  count_cols <- intersect(count_cols, names(df))
  key_cols   <- intersect(key_cols,   names(df))
  if (is.null(value_cols)) {
    value_cols <- setdiff(names(df), c(count_cols, key_cols))
  } else {
    value_cols <- intersect(value_cols, names(df))
  }

  if (length(count_cols) == 0L) {
    is_small <- rep(FALSE, nrow(df))
  } else {
    counts_mat <- vapply(count_cols, function(cc) {
      x <- suppressWarnings(as.numeric(df[[cc]]))
      ifelse(is.na(x), Inf, x)          # NA counts -> not small
    }, numeric(nrow(df)))
    if (is.null(dim(counts_mat))) counts_mat <- matrix(counts_mat, ncol = 1)
    is_small <- apply(counts_mat, 1L, function(r) any(r < min_cell))
  }
  is_small[is.na(is_small)] <- FALSE

  out <- df
  if (any(is_small)) {
    for (cc in c(count_cols, value_cols)) {
      if (is.factor(out[[cc]])) out[[cc]] <- as.character(out[[cc]])
      out[[cc]][is_small] <- if (is.numeric(out[[cc]])) NA else marker
    }
  }

  sup_path <- sub("\\.csv$", ".suppressed.csv", path)
  sup_log  <- data.frame(
    path            = basename(path),
    n_rows          = nrow(df),
    n_rows_redacted = sum(is_small),
    pct_redacted    = if (nrow(df)) round(100 * sum(is_small) / nrow(df), 2) else 0,
    min_cell        = min_cell,
    count_cols      = paste(count_cols, collapse = "|"),
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(out,     path,     row.names = FALSE)
  utils::write.csv(sup_log, sup_path, row.names = FALSE)
  message("Saved (suppressed k<", min_cell, "): ", path,
          if (sum(is_small)) sprintf("  [%d/%d rows redacted]", sum(is_small), nrow(df)) else "")
  invisible(out)
}


# ── Activity definitions ───────────────────────────────────

# ABCD 7.0 activity short names (SAIQ instrument)
activities_6 <- c(
 "baseb", "basketb", "fhock", "footb", "ihock", "lax", "rugby",
 "soccer", "swim", "horse", "gymnast", "tennis", "run", "mma",
 "volleyb", "dance", "climb", "martarts", "iskate", "skatebd",
 "surf", "yoga", "art", "drama", "craft", "games", "collect",
 "musicinstr"
)

# Mapping: analysis short name → 7.0 SAIQ prefix
act_map <- c(
 base    = "baseb",    basket  = "basketb",
 fhock   = "fhock",    fball   = "footb",
 ihock   = "ihock",    lax     = "lax",
 rugby   = "rugby",    soc     = "soccer",
 wpolo   = "swim",     polo    = "horse",
 gym     = "gymnast",  tennis  = "tennis",
 run     = "run",      mma     = "mma",
 vball   = "volleyb",  dance   = "dance",
 climb   = "climb",    m_arts  = "martarts",
 iskate  = "iskate",   sboard  = "skatebd",
 surf    = "surf",     yoga    = "yoga",
 art     = "art",      drama   = "drama",
 crafts  = "craft",    chess   = "games",
 collect = "collect",  music   = "musicinstr"
)

# Composite activity groups
Tsports_acts      <- c("base", "basket", "fball", "soc", "run")
richTsports_acts  <- c("fhock", "ihock", "lax", "rugby", "soc",
                        "wpolo", "tennis", "vball", "gym", "polo")
Isports_acts      <- c("dance", "climb", "m_arts", "iskate", "sboard",
                        "surf", "yoga", "mma")
sports_acts       <- c("fhock", "lax", "soc", "wpolo", "tennis", "run",
                        "vball", "gym", "polo", "dance", "climb", "m_arts",
                        "iskate", "sboard", "surf", "yoga", "mma")
Bsports_acts      <- c("base", "basket", "fball", "ihock")

# Standard covariate set
covars_base <- c("site", "age", "sex", "educ_prnt", "married", "income",
                 "ethnrace", "sch_achiev", "sch_slope", "sch_frl",
                 "pmq1", "pmq2", "pmq3", "pmq4", "pmq5", "fam")

# Extended covariate set. 2026-11-06: only FES conflict was measured at baseline,
# so it is the only FES variable in a baseline confounder set; the year-2
# subscales (fes_*_y2) are pre-exposure only for designs with time zero >= year 2.
covars_extended <- c(covars_base, "fes_confl")

# Additional outcome variables
outcomes_all <- c("reading", "reading_agecorr", "math", "grades")

# Variables for genetically informed methods
twin_vars <- c("zygosity_exp", "zygosity_gen", "pihat")
ancestry_pc_vars <- paste0("anc_pc", sprintf("%02d", 1:10))

# Genetic ancestry admixture proportions (from ancestry.tsv)
ancestry_prop_vars <- c("anc_EUR", "anc_EAS", "anc_AMR", "anc_AFR")

# Polygenic score variables (computed externally, loaded from PLINK .profile)
# pgs_dep = depression PGS (common_pgs_*); pgs_soc = social-connectedness PGS
# (soc_pgs_*). Used as NEGATIVE-CONTROL / falsification instruments, not as
# educational-attainment or intelligence instruments.
pgs_vars <- c("pgs_dep", "pgs_soc", "pgs_ea")   # pgs_ea: educational attainment (SSGAC 2022)

# Standard exposure composites to create
exposure_names <- c("music_sc", "arts", "other", "chess_sc", "sports",
                    "Tsports", "Isports", "richTsports", "Bsports")


# ── Helper functions ───────────────────────────────────────

#' Latest-wave anchor: one row per subject at their LATEST wave with a
#' non-missing value of `var`. Used to anchor cross-sectional / fixed-outcome
#' analyses on each subject's most recent observed outcome. ABCD 7.0 extends
#' the annual waves through 7yr, so the outcome wave is no longer a fixed
#' "4yr" — this picks up 7yr where present and falls back per subject.
#' @param d   long data frame with columns `id`, `event` (ordered factor), `var`
#' @param var outcome column to anchor on (default "reading")
#' @return one row per id; subjects with no non-NA `var` are dropped
anchor_latest <- function(d, var = "reading") {
 stopifnot("id" %in% names(d), "event" %in% names(d), var %in% names(d))
 d <- d[!is.na(d[[var]]), , drop = FALSE]
 if (!nrow(d)) return(d)
 ord <- if (is.factor(d$event)) as.integer(d$event)
        else match(as.character(d$event), event_levels_global())
 d <- d[order(d$id, -ord), , drop = FALSE]
 d[!duplicated(d$id), , drop = FALSE]
}

#' Canonical wave ordering, used by anchor_latest() when `event` is not a
#' factor. Mirrors the event_levels in 01_data.R.
event_levels_global <- function() {
 c("baseline", "1yr", "2yr", "3yr", "4yr", "5yr", "6yr", "7yr")
}

#' Activity hours accumulated strictly BEFORE each subject's outcome wave.
#' Pairs with anchor_latest(): the outcome wave is the subject's latest wave with
#' a non-missing `outcome`, and the exposure sums `var` over earlier waves only,
#' so music recorded at or after the test cannot enter (audit finding F3).
#' Subjects whose outcome wave is their first wave get NA. Added 2026-10-19.
#' @return data.frame(id, outcome_wave, <var>_pre, <var>_pre_mean, n_pre_waves)
#' Join each subject's baseline reading (reading_bl) onto a cross-section by id.
#' Used by the "+ baseline reading" ladder rows (1b, 3b, 4b, 8c), which refit a
#' module's main model with and without reading_bl on the same children so the
#' change is the adjustment, not the sample. Added 2026-10-19.
add_baseline_reading <- function(d, dl) {
 bl <- dl[as.character(dl$event) == "baseline" & !is.na(dl$reading), c("id", "reading")]
 bl <- bl[!duplicated(bl$id), ]
 d$reading_bl <- suppressWarnings(as.numeric(bl$reading[match(as.character(d$id), as.character(bl$id))]))
 d
}

pre_outcome_exposure <- function(d, var = "music", outcome = "reading") {
 ord <- function(e) if (is.factor(e)) as.integer(e) else match(as.character(e), event_levels_global())
 anc <- anchor_latest(d, outcome)
 ow  <- stats::setNames(ord(anc$event), as.character(anc$id))
 x <- d[as.character(d$id) %in% names(ow), c("id", "event", var)]
 x$w  <- ord(x$event)
 x$ow <- ow[as.character(x$id)]
 x <- x[!is.na(x$w) & x$w < x$ow & !is.na(x[[var]]), , drop = FALSE]
 agg <- x |>
   dplyr::group_by(id = as.character(id)) |>
   dplyr::summarise(pre = sum(.data[[var]]), pre_mean = mean(.data[[var]]),
                    n_pre_waves = dplyr::n(), .groups = "drop")
 out <- dplyr::left_join(data.frame(id = names(ow), outcome_wave = unname(ow),
                                    stringsAsFactors = FALSE), agg, by = "id")
 out$n_pre_waves[is.na(out$n_pre_waves)] <- 0L
 names(out)[names(out) == "pre"]      <- paste0(var, "_pre")
 names(out)[names(out) == "pre_mean"] <- paste0(var, "_pre_mean")
 out
}

#' Create scaled composite exposure variables from activity columns
#' @param df data frame with activity columns
#' @return df with additional composite columns
add_composites <- function(df) {
 # Wrap scale() in as.numeric() to return vectors (not 1-column matrices),
 # which prevents downstream issues with model.matrix() and as.matrix().
 df |>
   mutate(
     Tsports     = as.numeric(scale(rowSums(across(all_of(Tsports_acts)), na.rm = TRUE))),
     richTsports = as.numeric(scale(rowSums(across(all_of(richTsports_acts)), na.rm = TRUE))),
     Isports     = as.numeric(scale(rowSums(across(all_of(Isports_acts)), na.rm = TRUE))),
     sports      = as.numeric(scale(rowSums(across(all_of(sports_acts)), na.rm = TRUE))),
     Bsports     = as.numeric(scale(rowSums(across(all_of(Bsports_acts)), na.rm = TRUE))),
     arts        = as.numeric(scale(art + drama)),
     other       = as.numeric(scale(collect + chess + crafts)),
     chess_sc    = as.numeric(scale(chess)),
     music_sc    = as.numeric(scale(music))
   )
}

#' Add baseline-level composite columns (suffixed _bl)
#' @param df data frame with activity_bl columns
#' @return df with additional baseline composite columns
add_composites_bl <- function(df) {
 Tsports_bl_acts <- paste0(Tsports_acts, "_bl")
 Isports_bl_acts <- paste0(Isports_acts, "_bl")

 df |>
   mutate(
     Tsports_bl  = as.numeric(scale(rowSums(across(all_of(Tsports_bl_acts)), na.rm = TRUE))),
     Isports_bl  = as.numeric(scale(rowSums(across(all_of(Isports_bl_acts)), na.rm = TRUE))),
     sports_bl   = as.numeric(scale(rowSums(across(all_of(c(Tsports_bl_acts, Isports_bl_acts))),
                                             na.rm = TRUE))),
     arts_bl     = as.numeric(scale(art_bl + drama_bl)),
     other_bl    = as.numeric(scale(collect_bl + chess_bl + crafts_bl)),
     chess_sc_bl = as.numeric(scale(chess_bl)),
     music_sc_bl = as.numeric(scale(music_bl))
   )
}

#' Run GRF causal forest for a set of exposures and one outcome
#' @param dat data frame (model-matrix ready)
#' @param exps character vector of exposure column names
#' @param outcome character: "reading" or "math"
#' @param fam_col character: family clustering column
#' @param covar_cols integer or column indices for X matrix
#' @return data frame of ATE results
run_grf_battery <- function(dat, exps, outcome, fam_col = "fam", covar_cols) {
 results <- data.frame(Variable = character(), ATE_Estimate = double(),
                       ATE_Std_Error = double(), stringsAsFactors = FALSE)
 # dat_cumul / dat_blctrl can now contain rows with NA outcomes (reading or
 # math) — 01_data.R relaxed complete.cases so math (3yr-only) and reading
 # (baseline/2yr/4yr) can be populated independently. Drop rows with NA
 # outcome here so causal_forest never sees NA Y.
 keep   <- !is.na(dat[[outcome]])
 dat_ok <- dat[keep, , drop = FALSE]
 for (exp in exps) {
   W   <- as.vector(dat_ok[, exp])
   Y   <- as.vector(dat_ok[, outcome])
   X   <- as.matrix(dat_ok[, covar_cols])
   fam <- dat_ok[[fam_col]]

   model    <- causal_forest(X = X, Y = Y, W = W, clusters = fam,
                             compute.oob.predictions = TRUE)
   ate      <- average_treatment_effect(model)

   results <- rbind(results,
                    data.frame(Variable = exp,
                               ATE_Estimate  = ate[1],
                               ATE_Std_Error = ate[2]))
 }
 results$Lower_CI <- results$ATE_Estimate - qnorm(0.975) * results$ATE_Std_Error
 results$Upper_CI <- results$ATE_Estimate + qnorm(0.975) * results$ATE_Std_Error
 results
}

#' Save an RDS file to the cache directory with a message
save_rds <- function(obj, filename) {
 fp <- file.path(paths$cache, filename)
 saveRDS(obj, fp)
 message("Saved: ", fp)
}

#' Load an RDS file from the cache directory
load_rds <- function(filename) {
 fp <- file.path(paths$cache, filename)
 if (!file.exists(fp)) stop("Cache file not found: ", fp)
 readRDS(fp)
}

#' Canonical participant key for joins across ID conventions.
#' ABCD 7.0 tabulated data use "sub-XXXXXXXX" (8 characters, no NDARINV);
#' PLINK genotype-derived files (PGS .profile/.sscore, possibly ancestry.tsv)
#' use "NDAR_INVXXXXXXXX". Strips a leading "sub-", removes underscores,
#' upper-cases, then strips a leading "NDARINV", so both reduce to the same
#' 8-character key (confirmed by the masked-pattern diagnostic in
#' scripts/check_pgs_load.R, 2026-10-19). A raw-string join between the two
#' conventions matches nobody: the reason the 2026-07-24 PGS-IV ran in
#' ancestry-PC fallback mode (technical report finding F1).
canon_id <- function(x)
  sub("^NDARINV", "", toupper(gsub("_", "", sub("^sub-", "", as.character(x)))))

#' Report (aggregate counts only) how many participants a genetic file matched.
#' Warns below 50%: that pattern means an ID-format mismatch, not missing data.
report_id_match <- function(label, dat_keys, file_keys) {
  u <- unique(dat_keys); m <- sum(u %in% file_keys)
  message(sprintf("  ID match %-30s %d of %d participants (%.1f%%); file has %d IDs",
                  label, m, length(u), 100 * m / max(1, length(u)), length(unique(file_keys))))
  if (m < 0.5 * length(u))
    warning(sprintf("[01_data] only %.1f%% of participants matched %s: check the ID format",
                    100 * m / max(1, length(u)), label), call. = FALSE)
  invisible(m)
}

#' Load genetic ancestry admixture proportions from ancestry.tsv
#' @return data.frame with columns: id, anc_EUR, anc_EAS, anc_AMR, anc_AFR
load_ancestry <- function() {
 fp <- paths$ancestry
 if (!file.exists(fp)) {
   message("Ancestry file not found: ", fp)
   return(NULL)
 }
 anc <- read.delim(fp, stringsAsFactors = FALSE)
 names(anc)[names(anc) == "ID"]  <- "id"
 names(anc)[names(anc) == "EUR"] <- "anc_EUR"
 names(anc)[names(anc) == "EAS"] <- "anc_EAS"
 names(anc)[names(anc) == "AMR"] <- "anc_AMR"
 names(anc)[names(anc) == "AFR"] <- "anc_AFR"
 anc$id <- canon_id(anc$id)                      # join key; see canon_id()
 message("Loaded ancestry data: ", nrow(anc), " subjects")
 anc
}

#' Read a single PLINK score file: PLINK 1 ".profile" (FID IID PHENO CNT CNT2
#' SCORE) or PLINK 2 ".sscore" (#FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM
#' SCORE1_AVG). Keys on IID, returned as canon_id(IID); takes SCORE, else
#' SCORE1_AVG, else the first SCORE* column. The PLINK 2 header starts with
#' "#", so comment.char must be "" or read.table drops the header line.
#' @param fp path to a .profile file
#' @param score_name name to give the returned (standardized) score column
#' @return data.frame(id, <score_name>) or NULL if unreadable
read_plink_profile <- function(fp, score_name) {
 if (is.na(fp) || !file.exists(fp)) {
   message("  PGS file not found for ", score_name, ": ", fp)
   return(NULL)
 }
 prof <- read.table(fp, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, fill = TRUE, comment.char = "")
 names(prof) <- sub("^#", "", names(prof))
 id_col    <- if ("IID" %in% names(prof)) "IID" else names(prof)[1]
 score_col <- if ("SCORE" %in% names(prof)) "SCORE" else
              if ("SCORE1_AVG" %in% names(prof)) "SCORE1_AVG" else {
                sc <- grep("^SCORE", names(prof), value = TRUE)
                if (length(sc)) sc[1] else names(prof)[ncol(prof)] }
 df <- data.frame(id = canon_id(prof[[id_col]]),
                  score = as.numeric(scale(as.numeric(prof[[score_col]]))),
                  stringsAsFactors = FALSE)
 names(df)[2] <- score_name
 message(sprintf("  Loaded %s from %s: %d subjects (id=%s, score=%s)",
                 score_name, basename(fp), nrow(df), id_col, score_col))
 df
}

#' Load polygenic scores from the PLINK .profile files named in paths$pgs.
#' paths$pgs is a named vector: analytic column name -> .profile path.
#' For this project: pgs_dep = depression PGS (common_pgs_*),
#'                   pgs_soc = social-connectedness PGS (soc_pgs_*).
#' Scores are z-standardized within file. NOTE: these are confound-/falsification
#' scores, NOT educational-attainment or intelligence instruments.
#' @return data.frame with columns: id, pgs_dep, pgs_soc (or NULL if none load)
load_pgs <- function() {
 fps <- paths$pgs
 if (length(fps) == 0 || all(is.na(fps)) || is.null(names(fps))) {
   message("PGS files not configured (paths$pgs is NA) — PGS columns set to NA.")
   return(NULL)
 }
 parts <- Filter(Negate(is.null),
                 lapply(names(fps), function(nm) read_plink_profile(fps[[nm]], nm)))
 if (length(parts) == 0) {
   message("No PGS files could be loaded.")
   return(NULL)
 }
 pgs <- Reduce(function(a, b) dplyr::full_join(a, b, by = "id"), parts)
 message("Loaded PGS data: ", nrow(pgs), " subjects, ",
         length(parts), " score(s): ",
         paste(setdiff(names(pgs), "id"), collapse = ", "))
 pgs
}

cat("✓ 00_utils.R loaded\n")
