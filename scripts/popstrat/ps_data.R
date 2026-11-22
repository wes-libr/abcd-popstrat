########################################################################
##  ps_data.R — extend the Marek analytic object with what the PopStrat
##  analyses need, cached in dir_derived/popstrat_data.rds (host-only).
##
##  Adds to md$base (baseline, ses-00A):
##    race     : legacy 5-level ethno-racial identity (ab_g_stc__cohort_ethnrace__leg),
##               used ONLY to set the diversity proportion of resampled
##               analytic samples, as in the 2024 draft
##    birth_id, zyg  : twin/sibling structure for the within-family analyses
##    anc_EUR/AFR/AMR/EAS : ADMIXTURE proportions (ancestry.tsv)
##    motion_dmri, motion_nback, frames_rs, motion_rs_raw : motion QC
##    IDP sets (feat): smri  = Desikan thickness + area + aseg volume
##                     dti   = AtlasTrack FA (40 tracts)
##                     rs    = Gordon network-pair RSFC (from mr_data.R)
##                     nback = Desikan 2-back vs 0-back betas
########################################################################

ps_build <- function(md, cache_file = file.path(paths$cache, "popstrat_data.rds"), rebuild = FALSE) {
  if (file.exists(cache_file) && !rebuild) { message("Using cached ", cache_file); return(readRDS(cache_file)) }
  suppressPackageStartupMessages(library(NBDCtools))
  nm <- as.data.frame(NBDCtools::get_dd("abcd", release = "7.0"))$name
  pick <- function(re) grep(re, nm, value = TRUE)
  smri <- c(pick("^mr_y_smri__thk__dsk__[a-z0-9]+__[lr]h_mean$"),
            pick("^mr_y_smri__area__dsk__[a-z0-9]+__[lr]h_sum$"),
            pick("^mr_y_smri__vol__aseg__(ab|ag|cd|hc|pl|pt|tp|vdc)__[lr]h_sum$"))
  dti   <- pick("^mr_y_dti__fs__fa__at__[a-z0-9]+__[lr]h_wmean$")
  nback <- pick("^mr_y_tfmri__nback__2bv0b__dsk__[a-z0-9]+__[lr]h_beta$")
  extra <- c("ab_g_stc__cohort_ethnrace__leg", "ab_g_stc__design_id__birth", "gn_y_genrel_zyg__01",
             "mr_y_qc__incl__dmri_indicator", "mr_y_qc__incl__tfmri__nback_indicator",
             "mr_y_qc__mot__dmri__mot_mean", "mr_y_qc__mot__tfmri__nback__mot_mean",
             "mr_y_qc__mot__rsfmri__mot_mean", "mr_y_qc__mot__rsfmri__vol__censor_count")
  want <- intersect(c(extra, smri, dti, nback), nm)
  message("Pulling ", length(want), " additional variables (first run only)...")
  raw <- as.data.frame(NBDCtools::create_dataset_abcd(dir_data = paths$abcd_data, vars = want,
                                                      tables = "ab_g_stc", value_to_na = TRUE))
  to_num <- function(x) {
    if (is.numeric(x)) return(x)
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (sum(!is.na(x)) && mean(is.na(y[!is.na(x)])) > 0.5) y <- as.numeric(factor(x))
    y
  }
  raw$id <- as.character(raw$participant_id)
  static <- intersect(c("ab_g_stc__cohort_ethnrace__leg", "ab_g_stc__design_id__birth", "gn_y_genrel_zyg__01"), names(raw))
  st <- raw[!duplicated(raw$id), "id", drop = FALSE]
  for (v in static) st[[v]] <- tapply(as.character(raw[[v]]), raw$id, function(x) x[!is.na(x)][1])[st$id]
  b <- raw[raw$session_id %in% "ses-00A", , drop = FALSE]

  d <- md$base
  i <- match(d$id, b$id); j <- match(d$id, st$id)
  ## race: legacy codes 1 Hispanic, 2 White, 3 Black, 4 Asian, 13 Other (labels if NBDCtools returns a factor)
  rc <- as.character(st$ab_g_stc__cohort_ethnrace__leg[j])
  rc <- c("1" = "Hispanic", "2" = "White", "3" = "Black", "4" = "Asian", "13" = "Other")[rc] %||% rc
  d$race <- ifelse(is.na(rc), as.character(st$ab_g_stc__cohort_ethnrace__leg[j]), rc)
  d$birth_id <- st$ab_g_stc__design_id__birth[j]
  d$zyg <- to_num(st$gn_y_genrel_zyg__01[j])
  d$motion_dmri  <- to_num(b$mr_y_qc__mot__dmri__mot_mean[i])
  d$motion_nback <- to_num(b$mr_y_qc__mot__tfmri__nback__mot_mean[i])
  d$motion_rs_raw <- to_num(b$mr_y_qc__mot__rsfmri__mot_mean[i])
  d$frames_rs <- to_num(b$mr_y_qc__mot__rsfmri__vol__censor_count[i])
  for (v in c(smri, dti, nback)) if (v %in% names(b)) d[[v]] <- to_num(b[[v]][i])
  smri <- intersect(smri, names(d)); dti <- intersect(dti, names(d)); nback <- intersect(nback, names(d))
  d$ok_smri  <- d$ok_ct & stats::complete.cases(d[, smri])
  d$ok_dti   <- to_num(b$mr_y_qc__incl__dmri_indicator[i]) %in% 1 & stats::complete.cases(d[, dti])
  d$ok_nback <- to_num(b$mr_y_qc__incl__tfmri__nback_indicator[i]) %in% 1 & stats::complete.cases(d[, nback])

  ## ADMIXTURE proportions (participant-level; joined through canon_id)
  anc <- if (exists("load_ancestry")) load_ancestry() else NULL
  if (!is.null(anc)) {
    k <- match(canon_id(d$id), anc$id)
    for (v in c("anc_EUR", "anc_AFR", "anc_AMR", "anc_EAS")) d[[v]] <- if (v %in% names(anc)) anc[[v]][k] else NA_real_
  } else for (v in c("anc_EUR", "anc_AFR", "anc_AMR", "anc_EAS")) d[[v]] <- NA_real_

  pd <- md; pd$base <- d
  pd$feat <- c(md$feat, list(smri = smri, dti = dti, nback = nback))
  pd$y2 <- NULL
  saveRDS(pd, cache_file)
  message("Cached PopStrat analytic object (participant-level, host-only): ", cache_file)
  pd
}

`%||%` <- function(a, b) { a[is.na(a)] <- b[is.na(a)]; a }
