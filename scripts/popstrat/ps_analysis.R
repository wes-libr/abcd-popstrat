########################################################################
##  ps_analysis.R — the PopStrat paper's analyses on ABCD 7.0, plus the
##  extensions prompted by Marek et al. (2026). Uses helpers from
##  mr_analysis.R (make_Z, resid_on, ridge_fit_predict, prep_resid,
##  abs_q, mr_safe). Aggregate output only.
##
##  A  Diversity-proportion resampling (draft Figs 2-4): variance in each
##     IDP explained by genetic PCs 1-10; % reduction in IDP-VI association
##     strength when the PCs are added; similarity of each VI map to the
##     genetic-PC1 map — by diversity proportion 0, .25, .5, .75, 1.
##  B  Adjustment ladder in the full sample: cov / + family SES (income,
##     caregiver education) / + genetic PCs / + both; observed attenuation
##     vs the attenuation expected if the PCs acted only through the VI;
##     structure beyond every SES measure; ridge r_oos under each set.
##  C  Within-family (sibling fixed effects) brain-phenotype maps vs
##     between-family and PC-adjusted maps.
##  D  Genetically similar subsamples (ADMIXTURE EUR >= .9; AFR >= .5),
##     compared with size-matched draws from the full sample.
##  E  Head motion: SES-motion correlation, selection into the RSFC
##     sample, and the SES map under stricter motion control.
########################################################################

ps_opts <- modifyList(list(
  seed = 20261121, props = c(0, .25, .5, .75, 1), n_div = 2000, n_div_boot = 100,
  modalities = c("smri", "dti", "rs", "nback"), n_match_boot = 100,
  eur_cut = 0.9, afr_cut = 0.5, min_sub = 300
), if (exists("ps_opts_user")) ps_opts_user else list())

ps_vis <- c(ses = "le_l_coi__addr1__se__total__national_zscore", adi = "le_l_adi__addr1__national_prcnt",
            fiq = "nc_y_nihtb__comp__fluid__uncor_score", iq = "nc_y_nihtb__comp__tot__uncor_score",
            internalizing = "mh_p_cbcl__synd__int_tscore", height = "ph_y_anthr__height_mean",
            weight = "ph_y_anthr__weight_mean")
fam_ses <- c("ab_g_dyn__cohort_income__hhold__6lvl", "ab_g_dyn__cohort_edu__cgs")

ps_motion <- c(rs = "motion", dti = "motion_dmri", nback = "motion_nback", smri = NA, ct = NA)
ps_terms <- function(adj, mod, extra = NULL) {
  base <- c("age", "sex", "site", "manuf", if (!is.na(ps_motion[[mod]])) ps_motion[[mod]])
  t <- switch(adj, raw = character(0), cov = base, fam = c(base, fam_ses),
              gpc = c(base, gpc_names()), famgpc = c(base, fam_ses, gpc_names()))
  c(t, extra)
}

## partial correlation of VI with each IDP after residualizing both on
## `terms` (Frisch-Waugh). Partial r rather than standardized b: b inflates
## with estimation noise when an adjustment shrinks the VI's residual
## variance, which made correlated adjustors look like they *strengthen*
## associations in testing.
assoc_t <- function(d, feats, y, terms) {
  X <- as.matrix(d[, feats, drop = FALSE]); yy <- as.numeric(d[[y]])
  if (length(terms)) { Z <- make_Z(d, terms); X <- resid_on(X, Z); yy <- as.numeric(resid_on(yy, Z)) }
  r <- suppressWarnings(as.numeric(stats::cor(X, yy)))
  list(r = stats::setNames(r, feats), n = nrow(d) - if (length(terms)) ncol(make_Z(d, terms)) else 1)
}

## association strength above the noise floor: sum |r| minus its expectation
## under no association, p * sqrt(2 / (pi * n_eff))
strength <- function(a) sum(abs(a$r), na.rm = TRUE) - sum(!is.na(a$r)) * sqrt(2 / (pi * max(a$n, 1)))

## rows complete for the IDPs, y and every term in `term_sets`
rows_ok <- function(d, feats, y, term_sets) {
  v <- unique(c(feats, y, unlist(term_sets)))
  stats::complete.cases(d[, intersect(v, names(d)), drop = FALSE])
}

## variance in each IDP explained by `add` beyond `base` terms
dR2 <- function(d, feats, base, add) {
  X <- as.matrix(d[, feats, drop = FALSE]); tss <- colSums(sweep(X, 2, colMeans(X))^2)
  r0 <- resid_on(X, make_Z(d, base)); r1 <- resid_on(X, make_Z(d, c(base, add)))
  (colSums(r0^2) - colSums(r1^2)) / tss
}

ridge_oos_t <- function(dtr, dte, feats, y, terms) {
  P <- prep_resid(dtr, dte, feats, y, terms)
  if (stats::sd(P$ytr) == 0) return(NA_real_)
  suppressWarnings(stats::cor(ridge_fit_predict(P$Xtr, P$ytr, P$Xte), P$yte))
}

nz <- function(x) if (is.null(x) || !length(x)) NA_real_ else x
## null SD of strength() for p independent features: sqrt(p * (1 - 2/pi) / n)
strength_null_sd <- function(a) sqrt(sum(!is.na(a$r)) * (1 - 2 / pi) / max(a$n, 1))
## % reduction in strength; NA unless the reference strength is > 3 null SDs
## (otherwise the ratio is noise over noise)
pct_red <- function(a0, a1) {
  s0 <- strength(a0)
  if (!is.finite(s0) || s0 <= 3 * strength_null_sd(a0)) return(NA_real_)
  100 * (1 - strength(a1) / s0)
}
qsum <- function(x, nm) stats::setNames(c(mean(x, na.rm = TRUE), stats::quantile(x, c(.025, .975), na.rm = TRUE)),
                                        paste0(nm, c("_mean", "_lo", "_hi")))

########################################################################
ps_run <- function(pd, out_dir, write = function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE) }) {
  set.seed(ps_opts$seed); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- pd$base; feats <- pd$feat
  mods <- intersect(ps_opts$modalities, names(feats)); mods <- mods[vapply(mods, function(m) length(feats[[m]]) > 0, TRUE)]
  V <- ps_vis[ps_vis %in% names(d)]
  sub_mod <- function(x, m) x[x[[paste0("ok_", m)]] %in% TRUE, , drop = FALSE]
  d$nonwhite <- ifelse(is.na(d$race), NA, d$race != "White")

  ## ------------------------------------------------------------ sample
  cat("\n[PS0] Sample\n")
  smp <- do.call(rbind, lapply(mods, function(m) {
    x <- sub_mod(d, m)
    data.frame(modality = m, n_features = length(feats[[m]]), n = nrow(x),
               n_white = sum(x$race %in% "White"), n_nonwhite = sum(x$nonwhite %in% TRUE),
               n_gpc = sum(stats::complete.cases(x[, gpc_names()])),
               n_in_multichild_family = sum(duplicated(x$fam) | duplicated(x$fam, fromLast = TRUE)),
               n_eur90 = sum(x$anc_EUR >= ps_opts$eur_cut, na.rm = TRUE),
               n_afr50 = sum(x$anc_AFR >= ps_opts$afr_cut, na.rm = TRUE))
  }))
  write(smp, "ps_sample.csv")

  ## -------------------------------------- A. diversity-proportion resampling
  cat("[PS-A] Diversity-proportion resampling\n")
  A1 <- list(); A2 <- list(); A3 <- list()
  for (m in mods) {
    pool <- sub_mod(d, m)
    pool <- pool[!is.na(pool$nonwhite) & stats::complete.cases(pool[, c(gpc_names(), ps_terms("cov", m))]), ]
    W <- which(!pool$nonwhite); NW <- which(pool$nonwhite)
    for (p in ps_opts$props) {
      n_use <- min(ps_opts$n_div, length(W), length(NW))   # same n at every diversity level
      n_nw <- round(p * n_use); n_w <- n_use - n_nw
      res <- lapply(seq_len(ps_opts$n_div_boot), function(b) {
        s <- pool[c(sample(W, n_w), sample(NW, n_nw)), ]
        g <- dR2(s, feats[[m]], ps_terms("cov", m), gpc_names())
        gmap <- assoc_t(s, feats[[m]], "gpc01", ps_terms("cov", m))$r
        vi <- lapply(names(V), function(k) {
          ok <- !is.na(s[[V[[k]]]]); ss <- s[ok, ]
          if (nrow(ss) < 200 || stats::sd(ss[[V[[k]]]]) == 0) return(c(red = NA, sim = NA))
          a0 <- assoc_t(ss, feats[[m]], V[[k]], ps_terms("cov", m))
          a1 <- assoc_t(ss, feats[[m]], V[[k]], ps_terms("gpc", m))
          c(red = pct_red(a0, a1), sim = stats::cor(a0$r, gmap[names(a0$r)]))
        })
        list(dr2 = c(mean = mean(g), p50 = stats::median(g), p75 = abs_q(g, .75), p99 = abs_q(g, .99)),
             vi = do.call(rbind, vi))
      })
      dr <- do.call(rbind, lapply(res, `[[`, "dr2"))
      A1[[length(A1) + 1]] <- data.frame(modality = m, diversity = p, n = n_w + n_nw, n_nonwhite = n_nw,
                                         t(apply(dr, 2, mean)), dR2_p99_lo = stats::quantile(dr[, "p99"], .025),
                                         dR2_p99_hi = stats::quantile(dr[, "p99"], .975), n_boot = length(res))
      for (j in seq_along(V)) {
        red <- vapply(res, function(z) z$vi[j, "red"], 0); sim <- vapply(res, function(z) z$vi[j, "sim"], 0)
        A2[[length(A2) + 1]] <- data.frame(modality = m, diversity = p, variable = names(V)[j], n = n_w + n_nw,
                                           t(qsum(red, "pct_reduction")))
        A3[[length(A3) + 1]] <- data.frame(modality = m, diversity = p, variable = names(V)[j], n = n_w + n_nw,
                                           t(qsum(sim, "r_with_gpc1_map")))
      }
      cat(sprintf("    %s diversity %.2f done\n", m, p))
    }
  }
  write(do.call(rbind, A1), "ps_A_gpc_dR2_by_diversity.csv")
  write(do.call(rbind, A2), "ps_A_effect_reduction_by_diversity.csv")
  write(do.call(rbind, A3), "ps_A_gpc1_similarity_by_diversity.csv")

  ## ------------------------------------------------ B. adjustment ladder
  cat("[PS-B] Adjustment ladder: cov / family SES / genetic PCs / both\n")
  B1 <- list(); B2 <- list(); B3 <- list()
  adjs <- c("cov", "fam", "gpc", "famgpc")
  for (m in mods) {
    dm <- sub_mod(d, m)
    for (k in names(V)) {
      ok <- rows_ok(dm, feats[[m]], V[[k]], lapply(adjs, ps_terms, mod = m)); s <- dm[ok, ]
      if (nrow(s) < 500) next
      A <- lapply(adjs, function(a) assoc_t(s, feats[[m]], V[[k]], ps_terms(a, m))); names(A) <- adjs
      B <- lapply(A, `[[`, "r")
      gmap <- assoc_t(s, feats[[m]], "gpc01", ps_terms("cov", m))$r
      ## R2 of the VI on the PCs beyond cov (and beyond cov + family SES)
      r2v <- function(base) {
        Z0 <- make_Z(s, base); Z1 <- make_Z(s, c(base, gpc_names())); y <- s[[V[[k]]]]
        e0 <- sum(resid_on(y, Z0)^2); e1 <- sum(resid_on(y, Z1)^2); (e0 - e1) / e0
      }
      R2c <- r2v(ps_terms("cov", m)); R2f <- r2v(ps_terms("fam", m))
      B1[[length(B1) + 1]] <- data.frame(modality = m, variable = k, n = nrow(s),
        strength_cov = strength(A$cov), strength_cov_null_sd = strength_null_sd(A$cov), q90_abs_r_cov = abs_q(B$cov, .9),
        pct_red_fam = pct_red(A$cov, A$fam), pct_red_gpc = pct_red(A$cov, A$gpc),
        pct_red_famgpc = pct_red(A$cov, A$famgpc), pct_red_gpc_given_fam = pct_red(A$fam, A$famgpc),
        partialR2_vi_on_gpc = R2c, expected_pct_red_gpc = 100 * (1 - sqrt(1 - R2c)),
        partialR2_vi_on_gpc_given_fam = R2f, expected_pct_red_gpc_given_fam = 100 * (1 - sqrt(1 - R2f)),
        r_map_gpc1_cov = stats::cor(B$cov, gmap), r_map_gpc1_fam = stats::cor(B$fam, gmap),
        r_map_gpc1_gpc = stats::cor(B$gpc, gmap), r_map_cov_vs_gpc = stats::cor(B$cov, B$gpc),
        r_map_cov_vs_fam = stats::cor(B$cov, B$fam))
    }
    ## population structure beyond every SES measure: dR2 of PCs on IDPs
    ## beyond cov + family SES + COI + ADI (null expectation ~ 10 / n)
    allses <- c(fam_ses, unname(V[intersect(c("ses", "adi"), names(V))]))
    ok <- rows_ok(dm, feats[[m]], NULL, list(ps_terms("famgpc", m), allses)); s <- dm[ok, ]
    g0 <- dR2(s, feats[[m]], ps_terms("cov", m), gpc_names())
    g1 <- dR2(s, feats[[m]], c(ps_terms("cov", m), allses), gpc_names())
    sesR <- dR2(s, feats[[m]], ps_terms("cov", m), allses)
    sesR_g <- dR2(s, feats[[m]], ps_terms("gpc", m), allses)
    B2[[length(B2) + 1]] <- data.frame(modality = m, n = nrow(s), null_expected = 10 / nrow(s),
      gpc_dR2_mean = mean(g0), gpc_dR2_p99 = abs_q(g0, .99),
      gpc_dR2_beyond_allSES_mean = mean(g1), gpc_dR2_beyond_allSES_p99 = abs_q(g1, .99),
      allSES_dR2_mean = mean(sesR), allSES_dR2_beyond_gpc_mean = mean(sesR_g),
      pct_features_gpc_beyond_allSES_gt_null = 100 * mean(g1 > 3 * 10 / nrow(s)))
    ## ridge, discovery -> replication, for SES and IQ under each adjustment
    for (k in intersect(c("ses", "iq", "fiq"), names(V))) {
      ok <- rows_ok(dm, feats[[m]], V[[k]], lapply(adjs, ps_terms, mod = m)); s <- dm[ok, ]
      tr <- s$split %in% "disc"
      for (a in adjs) B3[[length(B3) + 1]] <- data.frame(modality = m, variable = k, adjustment = a,
        r_oos = nz(mr_safe("ridge", ridge_oos_t(s[tr, ], s[!tr, ], feats[[m]], V[[k]], ps_terms(a, m)))),
        n_train = sum(tr), n_test = sum(!tr))
    }
    cat(sprintf("    %s done\n", m))
  }
  write(do.call(rbind, B1), "ps_B_adjustment_ladder.csv")
  write(do.call(rbind, B2), "ps_B_structure_beyond_SES.csv")
  write(do.call(rbind, B3), "ps_B_ridge_by_adjustment.csv")

  ## ------------------------------------------- C. within-family (sibling FE)
  cat("[PS-C] Within-family (sibling fixed-effects) maps\n")
  C1 <- list()
  for (m in mods) for (excl_mz in c(FALSE, TRUE)) {
    dm <- sub_mod(d, m)
    if (excl_mz) dm <- dm[!(dm$fam %in% dm$fam[dm$zyg %in% 1]), ]
    ## child-level variables only: COI and ADI are address-level, ~constant within family
    for (k in intersect(c("iq", "fiq", "internalizing", "height", "weight"), names(V))) {
      tc <- setdiff(ps_terms("cov", m), c("site", "manuf"))
      ok <- rows_ok(dm, feats[[m]], V[[k]], list(ps_terms("gpc", m))); s <- dm[ok, ]
      s <- s[s$fam %in% s$fam[duplicated(s$fam)], ]
      if (nrow(s) < 200) next
      y <- s[[V[[k]]]]; wsd <- stats::sd(y - stats::ave(y, s$fam))
      if (!is.finite(wsd) || wsd < 1e-8) next
      a_cov <- assoc_t(s, feats[[m]], V[[k]], ps_terms("cov", m)); b_cov <- a_cov$r
      a_gpc <- assoc_t(s, feats[[m]], V[[k]], ps_terms("gpc", m)); b_gpc <- a_gpc$r
      gmap <- assoc_t(s, feats[[m]], "gpc01", ps_terms("cov", m))$r
      ## demean IDPs, VI and the person-varying covariates within family
      w <- s; num <- c(feats[[m]], V[[k]], intersect(tc, c("age", ps_motion[[m]])))
      for (v in num) w[[v]] <- w[[v]] - stats::ave(w[[v]], w$fam)
      w$sex_w <- as.numeric(w$sex) - stats::ave(as.numeric(w$sex), w$fam)
      X <- as.matrix(w[, feats[[m]]]); yy <- w[[V[[k]]]]
      Z <- cbind(1, as.matrix(w[, c(intersect(tc, c("age", ps_motion[[m]])), "sex_w")]))
      X <- resid_on(X, Z); yy <- as.numeric(resid_on(yy, Z))
      b_w <- stats::setNames(suppressWarnings(as.numeric(stats::cor(X, yy))), feats[[m]])
      a_w <- list(r = b_w, n = nrow(s) - length(unique(s$fam)) - ncol(Z))
      C1[[length(C1) + 1]] <- data.frame(modality = m, sample = if (excl_mz) "siblings+DZ (MZ families excluded)" else "all multi-child families",
        variable = k, n = nrow(s), n_families = length(unique(s$fam)),
        within_family_sd_ratio = wsd / stats::sd(y),
        pct_red_within_vs_cov = pct_red(a_cov, a_w), pct_red_gpc_vs_cov = pct_red(a_cov, a_gpc),
        r_map_within_vs_cov = stats::cor(b_w, b_cov), r_map_within_vs_gpc = stats::cor(b_w, b_gpc),
        r_map_within_vs_gpc1 = stats::cor(b_w, gmap), r_map_cov_vs_gpc1 = stats::cor(b_cov, gmap))
    }
  }
  write(do.call(rbind, C1), "ps_C_within_family.csv")

  ## ---------------------------- D. genetically similar subsamples
  cat("[PS-D] Genetically similar subsamples vs size-matched full-sample draws\n")
  D1 <- list()
  groups <- list(EUR = function(x) x$anc_EUR >= ps_opts$eur_cut, AFR = function(x) x$anc_AFR >= ps_opts$afr_cut)
  for (m in mods) {
    dm <- sub_mod(d, m)
    for (k in intersect(c("ses", "adi", "iq", "fiq", "height"), names(V))) {
      ok <- rows_ok(dm, feats[[m]], V[[k]], list(ps_terms("gpc", m))); s <- dm[ok, ]
      full_cov <- assoc_t(s, feats[[m]], V[[k]], ps_terms("cov", m))$r
      full_gpc <- assoc_t(s, feats[[m]], V[[k]], ps_terms("gpc", m))$r
      for (g in names(groups)) {
        sel <- groups[[g]](s) %in% TRUE; if (sum(sel) < ps_opts$min_sub) next
        sg <- s[sel, ]
        ac <- assoc_t(sg, feats[[m]], V[[k]], ps_terms("cov", m)); bc <- ac$r
        ag <- assoc_t(sg, feats[[m]], V[[k]], ps_terms("gpc", m))
        mt <- replicate(ps_opts$n_match_boot, {
          sm <- s[sample.int(nrow(s), nrow(sg)), ]; strength(assoc_t(sm, feats[[m]], V[[k]], ps_terms("cov", m))) })
        D1[[length(D1) + 1]] <- data.frame(modality = m, group = g, variable = k, n = nrow(sg),
          sd_ratio_vs_full = stats::sd(sg[[V[[k]]]]) / stats::sd(s[[V[[k]]]]),
          strength_group_cov = strength(ac), strength_fullsample_sizematched = mean(mt),
          strength_ratio = strength(ac) / mean(mt), pct_red_gpc_in_group = pct_red(ac, ag),
          r_map_group_vs_full_cov = stats::cor(bc, full_cov), r_map_group_vs_full_gpc = stats::cor(bc, full_gpc))
      }
    }
  }
  if (length(D1)) write(do.call(rbind, D1), "ps_D_genetically_similar_subsamples.csv")

  ## ------------------------------------------------------ E. head motion
  if ("rs" %in% mods && "ses" %in% names(V)) {
    cat("[PS-E] Head motion and selection into the RSFC sample\n")
    sv <- V[["ses"]]; all <- d[!is.na(d[[sv]]), ]
    sdiff <- function(v) { a <- all[[v]][all$ok_rs %in% TRUE]; b <- all[[v]][!(all$ok_rs %in% TRUE)]
      (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / stats::sd(all[[v]], na.rm = TRUE) }
    E1 <- data.frame(term = c("r(SES, mean FD all frames)", "r(SES, frames retained)", "r(gPC1, mean FD all frames)",
                              "std. diff SES: RSFC-included minus excluded", "std. diff gPC1: included minus excluded",
                              "std. diff IQ: included minus excluded"),
                     value = c(stats::cor(all[[sv]], all$motion_rs_raw, use = "pairwise"),
                               stats::cor(all[[sv]], all$frames_rs, use = "pairwise"),
                               stats::cor(all$gpc01, all$motion_rs_raw, use = "pairwise"),
                               sdiff(sv), sdiff("gpc01"), if ("iq" %in% names(V)) sdiff(V[["iq"]]) else NA),
                     n = c(sum(stats::complete.cases(all[, c(sv, "motion_rs_raw")])),
                           sum(stats::complete.cases(all[, c(sv, "frames_rs")])),
                           sum(stats::complete.cases(all[, c("gpc01", "motion_rs_raw")])),
                           nrow(all), sum(!is.na(all$gpc01)), nrow(all)))
    write(E1, "ps_E_motion_selection.csv")
    dm <- sub_mod(d, "rs"); extra_m <- c("motion_rs_raw", "frames_rs")
    ok <- rows_ok(dm, feats$rs, sv, list(c(ps_terms("gpc", "rs"), extra_m))); s <- dm[ok, ]
    A_std <- assoc_t(s, feats$rs, sv, ps_terms("cov", "rs")); std <- A_std$r
    A_more <- assoc_t(s, feats$rs, sv, c(ps_terms("cov", "rs"), extra_m)); more <- A_more$r
    strict <- s[s$motion <= stats::median(s$motion), ]
    A_sb <- assoc_t(strict, feats$rs, sv, ps_terms("cov", "rs")); st_b <- A_sb$r
    A_sg <- assoc_t(strict, feats$rs, sv, ps_terms("gpc", "rs")); st_g <- A_sg$r
    mt <- replicate(ps_opts$n_match_boot, strength(assoc_t(s[sample.int(nrow(s), nrow(strict)), ], feats$rs, sv,
                                                          ps_terms("cov", "rs"))))
    E2 <- data.frame(analysis = c("standard", "+ raw mean FD + frames retained", "low-motion half", "low-motion half + genetic PCs"),
                     n = c(nrow(s), nrow(s), nrow(strict), nrow(strict)),
                     strength = c(strength(A_std), strength(A_more), strength(A_sb), strength(A_sg)),
                     strength_sizematched_full = c(NA, NA, mean(mt), NA),
                     r_map_with_standard = c(1, stats::cor(std, more), stats::cor(std, st_b), stats::cor(std, st_g)),
                     ms_share = vapply(list(std, more, st_b, st_g), function(b) net_profile(b, pd$rs_net)$ms_share, 0))
    write(E2, "ps_E_motion_ses_map.csv")
  }
  cat("\nDone. Aggregate outputs in ", out_dir, "\n")
  invisible(TRUE)
}
