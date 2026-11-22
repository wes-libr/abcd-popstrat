########################################################################
##  sc_data.R — longitudinal panel for the causal analysis of Marek et
##  al.'s exposure (neighborhood SES: COI social and economic domain) on
##  the brain and on IQ, ABCD 7.0 imaging waves 0, 2, 4, 6.
##  Cached host-only in dir_derived/ses_causal_panel.rds.
##
##  Outcomes
##    sig_rs : RSFC "SES signature" = each child's network-pair RSFC
##             projected on the baseline SES map (cov-adjusted r per pair)
##             estimated in the DISCOVERY half only. Causal analyses of the
##             signature use the REPLICATION half only, so the outcome
##             definition is independent of the estimation sample.
##    sig_ct : the same for Destrieux cortical thickness.
##    ms_rs  : mean RSFC of motor/sensory network pairs (smh, smm, vis,
##             aud): an a-priori, not data-derived, version of Marek's map.
##    iq     : NIH Toolbox total composite (uncorrected), at the waves it
##             was given; fiq = fluid composite.
##  Exposure: coi = le_l_coi__addr1__se__total__national_zscore at each
##  wave. Changes within child arise from moves (COI is fixed per tract).
########################################################################

sc_build <- function(md, cache_file = file.path(paths$cache, "ses_causal_panel.rds"), rebuild = FALSE) {
  if (file.exists(cache_file) && !rebuild) { message("Using cached ", cache_file); return(readRDS(cache_file)) }
  suppressPackageStartupMessages(library(NBDCtools))
  nm <- as.data.frame(NBDCtools::get_dd("abcd", release = "7.0"))$name
  rs <- md$feat$rs; ct <- md$feat$ct
  vars <- c("le_l_coi__addr1__se__total__national_zscore", "le_l_adi__addr1__national_prcnt",
            "nc_y_nihtb__comp__tot__uncor_score", "nc_y_nihtb__comp__fluid__uncor_score",
            "ab_g_dyn__cohort_income__hhold__6lvl", "ab_g_dyn__cohort_edu__cgs",
            "ab_g_dyn__design_site", "ab_g_dyn__visit_age", "mr_y_adm__info__dev_manufact",
            "mr_y_qc__incl__rsfmri_indicator", "mr_y_qc__incl__smri__t1_indicator",
            "mr_y_qc__mot__rsfmri__vol__censor_count", "mr_y_qc__mot__rsfmri__mot__censor_mean",
            "ph_y_anthr__height_mean", rs, ct)
  vars <- intersect(vars, nm)
  message("Pulling ", length(vars), " variables across imaging waves (first run only)...")
  raw <- as.data.frame(NBDCtools::create_dataset_abcd(dir_data = paths$abcd_data, vars = vars, value_to_na = TRUE))
  to_num <- function(x) {
    if (is.numeric(x)) return(x)
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (sum(!is.na(x)) && mean(is.na(y[!is.na(x)])) > 0.5) y <- as.numeric(factor(x))
    y
  }
  waves <- c("ses-00A" = 0, "ses-02A" = 2, "ses-04A" = 4, "ses-06A" = 6)
  x <- raw[as.character(raw$session_id) %in% names(waves), , drop = FALSE]
  b <- md$base
  p <- data.frame(id = as.character(x$participant_id), wave = unname(waves[as.character(x$session_id)]))
  k <- match(p$id, b$id)
  p$fam <- b$fam[k]; p$sex <- b$sex[k]; p$split <- b$split[k]
  for (g in grep("^gpc", names(b), value = TRUE)) p[[g]] <- b[[g]][k]
  p$site  <- factor(x$ab_g_dyn__design_site); p$age <- to_num(x$ab_g_dyn__visit_age)
  if (isTRUE(stats::median(p$age, na.rm = TRUE) > 50)) p$age <- p$age / 12
  p$manuf <- factor(x$mr_y_adm__info__dev_manufact)
  p$motion <- to_num(x$mr_y_qc__mot__rsfmri__mot__censor_mean)
  p$coi <- to_num(x$le_l_coi__addr1__se__total__national_zscore)
  p$adi <- to_num(x$le_l_adi__addr1__national_prcnt)
  p$iq  <- to_num(x$nc_y_nihtb__comp__tot__uncor_score)
  p$fiq <- to_num(x$nc_y_nihtb__comp__fluid__uncor_score)
  p$height <- to_num(x$ph_y_anthr__height_mean)
  ## family SES at baseline (time-fixed confounders)
  p$income <- b[["ab_g_dyn__cohort_income__hhold__6lvl"]][k]
  p$edu    <- b[["ab_g_dyn__cohort_edu__cgs"]][k]
  ok_rs <- to_num(x$mr_y_qc__incl__rsfmri_indicator) %in% 1 &
           (to_num(x$mr_y_qc__mot__rsfmri__vol__censor_count) >= 375) %in% TRUE
  ok_ct <- to_num(x$mr_y_qc__incl__smri__t1_indicator) %in% 1
  R <- sapply(rs, function(v) to_num(x[[v]])); C <- sapply(ct, function(v) to_num(x[[v]]))
  ok_rs <- ok_rs & stats::complete.cases(R); ok_ct <- ok_ct & stats::complete.cases(C)

  ## SES signature weights: baseline DISCOVERY half, cov-adjusted r per feature
  disc0 <- b[b$split %in% "disc", ]
  w_rs <- assoc_map(disc0[disc0$ok_rs %in% TRUE, ], rs, md$key[["ses"]], "rs", "cov")$r
  w_ct <- assoc_map(disc0[disc0$ok_ct %in% TRUE, ], ct, md$key[["ses"]], "ct", "cov")$r
  ## standardize features on the baseline discovery distribution, then project
  sig <- function(M, ok, ref, w) {
    mu <- colMeans(ref); s <- apply(ref, 2, stats::sd)
    z <- sweep(sweep(M, 2, mu), 2, s, "/"); out <- as.numeric(z %*% w) / sum(abs(w))
    out[!ok] <- NA; out
  }
  ref_rs <- as.matrix(disc0[disc0$ok_rs %in% TRUE, rs]); ref_ct <- as.matrix(disc0[disc0$ok_ct %in% TRUE, ct])
  p$sig_rs <- sig(R, ok_rs, ref_rs, w_rs); p$sig_ct <- sig(C, ok_ct, ref_ct, w_ct)
  msp <- md$rs_net$feature[md$rs_net$net1 %in% mr_opts$ms_nets & md$rs_net$net2 %in% mr_opts$ms_nets]
  p$ms_rs <- rowMeans(R[, msp, drop = FALSE]); p$ms_rs[!ok_rs] <- NA
  p$ok_rs <- ok_rs; p$ok_ct <- ok_ct
  p <- p[order(p$id, p$wave), ]
  saveRDS(list(panel = p, w_rs = w_rs, w_ct = w_ct, ms_pairs = msp), cache_file)
  message("Cached SES causal panel (participant-level, host-only): ", cache_file)
  list(panel = p, w_rs = w_rs, w_ct = w_ct, ms_pairs = msp)
}
