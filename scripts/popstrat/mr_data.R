########################################################################
##  mr_data.R — build the analytic object for the Marek replication
##  from the ABCD 7.0 tabulated release (NBDCtools), cached OUTSIDE the
##  shared repo in dir_derived/marek_rep_data.rds.
##
##  Brain features (what 7.0 tabulated data allow; see README):
##    rs : Gordon network-pair resting-state correlations, 13 networks,
##         91 unique pairs (mr_y_rsfmri__corr__gpnet). Marek used 55,278
##         parcel edges (333 Gordon parcels) from the DCAN pipeline, which
##         are not in the tabulated release.
##    ct : Destrieux cortical thickness, 148 regions (mr_y_smri__thk__dst).
##         Marek used 59,412 vertices.
##  QC: rsfMRI recommended for inclusion and >= 375 retained frames
##  (5 min at TR 0.8 s; ABCD censors at FD 0.2 mm, Marek at 0.08 mm);
##  T1 recommended for inclusion.
##  Variables of interest: an automated exposome/phenome built from the
##  table families Marek drew on (summary scores, not items), filtered in
##  mr_analysis.R by Marek's n >= 2000 rule.
########################################################################

mr_key <- c(
  ses    = "le_l_coi__addr1__se__total__national_zscore",  # COI social & economic domain (Marek's "SES")
  iq     = "nc_y_nihtb__comp__tot__uncor_score",           # NIH Toolbox total composite
  pfac   = "mh_p_cbcl_sum",                                 # CBCL total problems ("p-factor")
  sleep  = "ph_p_sds__dims_001",                            # hours of sleep on most nights (higher = less)
  screen = "screen_total",                                  # youth weekday + weekend screen time (derived)
  adi    = "le_l_adi__addr1__national_prcnt",
  income = "ab_g_dyn__cohort_income__hhold__6lvl",
  edu    = "ab_g_dyn__cohort_edu__cgs",
  height = "ph_y_anthr__height_mean")

## table families -> category (Marek's 12 categories, collapsed where 7.0 tables differ)
mr_tables <- list(
  socioeconomic = c("le_l_coi", "le_l_adi", "le_l_svi", "le_l_ssvi", "le_l_eci"),
  culture_environment = c("le_l_crime", "le_l_nbhsoc", "le_l_pm25", "le_l_no2", "le_l_leadrisk",
                          "le_l_urban", "le_l_walk", "fc_p_nsc", "fc_y_nsc", "fc_y_srpf", "fc_p_meim", "fc_y_meim"),
  screen_time = c("nt_y_stq", "nt_p_stq"),
  cognition = c("nc_y_nihtb", "nc_y_lmt", "nc_y_ravlt", "nc_y_wisc", "nc_y_cct"),
  mental_health = c("mh_p_cbcl", "mh_p_gbi", "mh_y_upps", "mh_y_bisbas", "mh_y_pps", "mh_p_ple", "mh_y_ple"),
  physical_health = c("ph_p_sds", "ph_y_anthr", "ph_y_pds", "ph_p_pds", "ph_p_dhx", "ph_p_mhx", "ph_y_pa"),
  family_parenting = c("fc_p_fes", "fc_y_fes", "fc_y_pm", "fc_y_crpbi", "fc_p_psb", "fc_y_psb"),
  substance = c("su_y_sui", "su_p_crpf", "su_y_crpf"))

mr_build <- function(cache_file = file.path(paths$cache, "marek_rep_data.rds"), rebuild = FALSE) {
  if (file.exists(cache_file) && !rebuild) { message("Using cached ", cache_file); return(readRDS(cache_file)) }
  suppressPackageStartupMessages(library(NBDCtools))
  dd <- as.data.frame(NBDCtools::get_dd("abcd", release = "7.0"))
  nm <- dd$name
  lab_col <- intersect(c("label", "label_en", "description"), names(dd))[1]
  lab <- if (is.na(lab_col)) rep("", nrow(dd)) else as.character(dd[[lab_col]])

  rs_vars <- grep("^mr_y_rsfmri__corr__gpnet__[a-z]{3}__[a-z]{3}_mean$", nm, value = TRUE)
  pr <- do.call(rbind, strsplit(sub("^mr_y_rsfmri__corr__gpnet__([a-z]{3})__([a-z]{3})_mean$", "\\1 \\2", rs_vars), " "))
  rs_vars <- rs_vars[pr[, 1] <= pr[, 2]]                          # 169 ordered -> 91 unique pairs
  rs_net <- data.frame(feature = rs_vars, net1 = pr[pr[, 1] <= pr[, 2], 1], net2 = pr[pr[, 1] <= pr[, 2], 2])
  ct_vars <- grep("^mr_y_smri__thk__dst__[a-z0-9]+__[lr]h_mean$", nm, value = TRUE)

  admin <- "(_dtt|_age|_lang|_dev|_adm.*|_nm|_ver|_dt|__v0[0-9]+)$"
  vi_rows <- list()
  for (cat_ in names(mr_tables)) for (tb in mr_tables[[cat_]]) {
    v <- grep(paste0("^", tb, "(_|$)"), nm, value = TRUE)
    v <- v[!grepl(admin, v)]
    if (startsWith(tb, "le_l_")) v <- v[grepl("addr1", v) & !grepl("metro|state|_c5$", v)]
    else {
      ## summary scores only: derived suffixes, no item codes (_001, ___3, __dk, __l)
      derived <- "(_sum|_mean|_score|_tscore|_zscore|_count|_rt|_prcnt|_prcntile|_rate|_total)$"
      item    <- "_[0-9]{3}([a-z]|__|_|$)|___[0-9]+|__dk$|__l$|_dk$"
      v <- v[(grepl(derived, v) & !grepl(item, v)) | v %in% mr_key]
    }
    if (tb == "nc_y_nihtb") v <- v[grepl("(uncor|agecor|agecorr)_score$|computed_score$", v)]
    if (length(v)) vi_rows[[length(vi_rows) + 1]] <- data.frame(variable = v, table = tb, category = cat_)
  }
  vi <- do.call(rbind, vi_rows)
  ## household income and caregiver education sit in ab_g_dyn
  vi <- rbind(vi, data.frame(variable = c(mr_key[["income"]], mr_key[["edu"]], mr_key[["sleep"]], "screen_total"),
                             table = c("ab_g_dyn", "ab_g_dyn", "ph_p_sds", "nt_y_stq"),
                             category = c("socioeconomic", "socioeconomic", "physical_health", "screen_time")))
  vi <- vi[!duplicated(vi$variable), ]
  vi$label <- substr(lab[match(vi$variable, nm)], 1, 120)

  design <- c("ab_g_dyn__design_site", "ab_g_dyn__visit_age", "mr_y_adm__info__dev_manufact",
              "mr_y_qc__incl__rsfmri_indicator", "mr_y_qc__incl__smri__t1_indicator",
              "mr_y_qc__mot__rsfmri__vol__censor_count", "mr_y_qc__mot__rsfmri__mot__censor_mean",
              "nt_y_stq__screen__wkdy_sum", "nt_y_stq__screen__wknd_sum")
  static <- c("ab_g_stc__cohort_sex", "ab_g_stc__design_id__fam", sprintf("gn_y_popstruct_pc__%02d", 1:10))
  want <- unique(c(design, static, rs_vars, ct_vars, setdiff(vi$variable, "screen_total")))
  miss <- setdiff(want, nm)
  if (length(miss)) message("Dropping ", length(miss), " names absent from the 7.0 dictionary: ",
                            paste(utils::head(miss, 20), collapse = ", "))
  want <- intersect(want, nm); vi <- vi[vi$variable %in% c(want, "screen_total"), ]

  message("Pulling ", length(want), " variables from ABCD 7.0 (first run only; slow)...")
  raw <- as.data.frame(NBDCtools::create_dataset_abcd(dir_data = paths$abcd_data, vars = want,
                                                      tables = "ab_g_stc", value_to_na = TRUE))
  to_num <- function(x) {
    if (is.numeric(x)) return(x)
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (sum(!is.na(x)) && mean(is.na(y[!is.na(x)])) > 0.5) y <- as.numeric(factor(x))
    y
  }
  raw$id <- as.character(raw$participant_id); raw$ses <- as.character(raw$session_id)
  ## static fields arrive on the Static (NA) session: carry them to every row
  for (v in intersect(static, names(raw))) {
    s <- tapply(as.character(raw[[v]]), raw$id, function(x) x[!is.na(x)][1])
    raw[[v]] <- unname(s[raw$id])
  }

  build_wave <- function(ses_id) {
    x <- raw[raw$ses %in% ses_id, , drop = FALSE]
    out <- data.frame(id = x$id, fam = x$ab_g_stc__design_id__fam,
                      site = factor(x$ab_g_dyn__design_site), age = to_num(x$ab_g_dyn__visit_age),
                      sex = factor(x$ab_g_stc__cohort_sex), manuf = factor(x$mr_y_adm__info__dev_manufact),
                      motion = to_num(x$mr_y_qc__mot__rsfmri__mot__censor_mean))
    for (k in 1:10) out[[sprintf("gpc%02d", k)]] <- to_num(x[[sprintf("gn_y_popstruct_pc__%02d", k)]])
    out$ok_rs <- to_num(x$mr_y_qc__incl__rsfmri_indicator) %in% 1 &
                 (to_num(x$mr_y_qc__mot__rsfmri__vol__censor_count) >= 375) %in% TRUE
    out$ok_ct <- to_num(x$mr_y_qc__incl__smri__t1_indicator) %in% 1
    for (v in c(rs_vars, ct_vars)) out[[v]] <- to_num(x[[v]])
    for (v in setdiff(vi$variable, "screen_total")) out[[v]] <- to_num(x[[v]])
    out$screen_total <- to_num(x$nt_y_stq__screen__wkdy_sum) + to_num(x$nt_y_stq__screen__wknd_sum)
    out$ok_rs <- out$ok_rs & stats::complete.cases(out[, rs_vars])
    out$ok_ct <- out$ok_ct & stats::complete.cases(out[, ct_vars])
    out
  }
  base <- build_wave("ses-00A"); y2 <- build_wave("ses-02A")

  ## linked external (le_l_*) and static variables: give year 2 the baseline
  ## value where year 2 has none (Marek used baseline linked data in year 2)
  carry <- vi$variable[startsWith(vi$table, "le_l_")]
  for (v in carry) {
    b <- base[[v]][match(y2$id, base$id)]; y2[[v]] <- ifelse(is.na(y2[[v]]), b, y2[[v]])
  }

  ## discovery / replication: random half of FAMILIES within site (Marek used
  ## the matched ARMS split of Feczko et al., not in the 7.0 release)
  set.seed(20261116)
  base$fam[is.na(base$fam)] <- paste0("solo_", base$id[is.na(base$fam)])
  fam_site <- base[!duplicated(base$fam), c("fam", "site")]
  fam_site$split <- NA_character_
  for (s in unique(as.character(fam_site$site))) {
    i <- which(as.character(fam_site$site) %in% s)
    fam_site$split[i] <- sample(rep(c("disc", "rep"), length.out = length(i)))
  }
  base$split <- fam_site$split[match(base$fam, fam_site$fam)]
  y2$split <- base$split[match(y2$id, base$id)]

  md <- list(base = base, y2 = y2, vi = vi, key = mr_key,
             feat = list(rs = rs_vars, ct = ct_vars), rs_net = rs_net,
             built = format(Sys.time()), release = "7.0")
  saveRDS(md, cache_file)
  message("Cached analytic object (participant-level, host-only): ", cache_file)
  md
}
