########################################################################
##  px_analysis.R — extension analyses for the PopStrat manuscript
##  (ABCD 7.0), answering the issues raised in the November 2026 review
##  of the draft. Uses helpers from mr_analysis.R and ps_analysis.R.
##  Aggregate output only.
##
##  X1  Map similarity, done properly (review items 2 and 6).
##      The draft's test (Fig. 5) compares a variable's brain map with the
##      genetic-PC1 map. That similarity is expected to be high even with
##      no confounding whenever the variable correlates with PC1: if PC1
##      reached the brain only through the variable, the PC1 map would be
##      the variable's map scaled by their correlation. Two corrections:
##        sim_unadj  r(VI map, PC1 map)                    (the draft's metric)
##        sim_expect the same similarity simulated under "PC1 reaches the
##                   brain only through the VI" (PC1* = rho*VI + noise)
##        sim_adj    r(VI map, map of PC1 residualized on the VI); about 0
##                   if PC1 reaches the brain only through the VI
##        z_gadj     strength of the residualized-PC1 map in null SDs
##      Permutation null for sim_adj (residualized PC1 permuted across
##      children), in every diversity resample and in the full sample;
##      family-clustered bootstrap CIs in the full sample.
##  X2  Head and body size (review item 4). Variance in IDPs explained by
##      genetic PCs before and after adjusting for intracranial volume,
##      height and weight, by sMRI measure type (thickness, area, volume)
##      and modality; whether the PC attenuation of COI and IQ associations
##      survives that adjustment; genetically similar subsamples by sMRI
##      measure type.
##  X3  Robustness of the population-structure adjustment (review item 5):
##      PCs 1-5, 1-10, all available; ADMIXTURE proportions; PCs + ADMIXTURE;
##      site x scanner (and site x serial number, if available) in place of
##      separate site and manufacturer terms. % reduction in association
##      strength relative to covariates, and given family SES.
##      PCs computed within ancestry groups need genotypes and are not run.
##  Figures (review item 1): px_figures() draws the 7.0 versions of draft
##      Figs 3-5 from ps_A_*.csv and X1, as PDFs from aggregate tables.
########################################################################

px_opts <- modifyList(list(
  seed = 20261122, props = c(0, .25, .5, .75, 1), n_div = 2000, n_div_boot = 100,
  n_perm_div = 100, n_perm_full = 1000, n_cboot = 200, n_expect = 20, n_match_boot = 100,
  vis = c("ses", "adi", "iq", "fiq", "internalizing", "height", "weight"),
  robust_vis = c("ses", "adi", "iq", "fiq", "height")
), if (exists("px_opts_user")) px_opts_user else list())

## ---------------------------------------------------------------- helpers
## standardized residual columns (mean 0, unit SD); zero-variance columns -> NA
px_std <- function(M) { M <- as.matrix(M); s <- apply(M, 2, stats::sd); M <- sweep(M, 2, colMeans(M)); sweep(M, 2, ifelse(s > 0, s, NA), "/") }

## correlation of each column of standardized Xs with a vector
px_r <- function(Xs, v) { v <- (v - mean(v)) / stats::sd(v); as.numeric(crossprod(Xs, v)) / (length(v) - 1) }

## X1 core: map-similarity statistics for one sample `s`
px_map_stats <- function(s, feats, vi, terms, n_perm, n_expect) {
  Z <- make_Z(s, terms)
  Xs <- px_std(resid_on(as.matrix(s[, feats, drop = FALSE]), Z))
  y <- as.numeric(resid_on(s[[vi]], Z)); g <- as.numeric(resid_on(s$gpc01, Z))
  if (stats::sd(y) < 1e-10 || stats::sd(g) < 1e-10) return(NULL)
  n_eff <- nrow(s) - ncol(Z); p <- length(feats)
  r_vi <- px_r(Xs, y); r_g <- px_r(Xs, g)
  g_perp <- as.numeric(stats::lm.fit(cbind(1, y), g)$residuals)
  r_ga <- px_r(Xs, g_perp)
  rho <- stats::cor(y, g)
  ## expected unadjusted similarity if PC1 reached the brain only through the VI
  ys <- (y - mean(y)) / stats::sd(y)
  sim_exp <- mean(replicate(n_expect, {
    gs <- rho * ys + sqrt(1 - rho^2) * stats::rnorm(length(ys)); stats::cor(r_vi, px_r(Xs, gs)) }))
  ## permutation null for sim_adj
  sim_adj <- stats::cor(r_vi, r_ga)
  null <- if (n_perm > 0) replicate(n_perm, stats::cor(r_vi, px_r(Xs, sample(g_perp)))) else NA_real_
  st <- function(r) sum(abs(r)) - p * sqrt(2 / (pi * n_eff))
  c(rho_vi_gpc1 = rho, sim_unadj = stats::cor(r_vi, r_g), sim_expect = sim_exp, sim_adj = sim_adj,
    sim_adj_null_q975 = if (n_perm > 0) as.numeric(stats::quantile(abs(null), .975)) else NA,
    sim_adj_perm_p = if (n_perm > 0) (1 + sum(abs(null) >= abs(sim_adj))) / (1 + n_perm) else NA,
    z_gadj = st(r_ga) / sqrt(p * (1 - 2 / pi) / n_eff), z_vi = st(r_vi) / sqrt(p * (1 - 2 / pi) / n_eff))
}

## % reduction in strength for an arbitrary pair of term sets on the same rows
px_red <- function(s, feats, vi, t0, t1) pct_red(assoc_t(s, feats, vi, t0), assoc_t(s, feats, vi, t1))

## proportion of variance in the VI explained by `add` beyond `base`
px_R2 <- function(s, vi, base, add) {
  y <- s[[vi]]; e0 <- sum(resid_on(y, make_Z(s, base))^2); e1 <- sum(resid_on(y, make_Z(s, c(base, add)))^2)
  (e0 - e1) / e0
}

########################################################################
px_run <- function(pd, out_dir, write = function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE) }) {
  set.seed(px_opts$seed); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- pd$base; feats <- pd$feat
  mods <- intersect(ps_opts$modalities, names(feats)); mods <- mods[vapply(mods, function(m) length(feats[[m]]) > 0, TRUE)]
  V <- ps_vis[intersect(px_opts$vis, names(ps_vis))]; V <- V[V %in% names(d)]
  sub_mod <- function(x, m) x[x[[paste0("ok_", m)]] %in% TRUE, , drop = FALSE]
  d$nonwhite <- ifelse(is.na(d$race), NA, d$race != "White")
  gall <- px_gpc_all(d)

  ## ---------------------------------------------- X1a. diversity resampling
  cat("\n[PX-1a] Map similarity by diversity proportion (", px_opts$n_div_boot, " resamples each)\n", sep = "")
  X1 <- list()
  for (m in mods) {
    pool <- sub_mod(d, m)
    pool <- pool[!is.na(pool$nonwhite) & stats::complete.cases(pool[, c("gpc01", ps_terms("cov", m))]), ]
    W <- which(!pool$nonwhite); NW <- which(pool$nonwhite)
    n_use <- min(px_opts$n_div, length(W), length(NW))
    for (p in px_opts$props) {
      n_nw <- round(p * n_use); n_w <- n_use - n_nw
      res <- lapply(seq_len(px_opts$n_div_boot), function(b) {
        s <- pool[c(sample(W, n_w), sample(NW, n_nw)), ]
        lapply(V, function(v) { ss <- s[!is.na(s[[v]]), ]
          if (nrow(ss) < 200) return(NULL)
          mr_safe("map stats", px_map_stats(ss, feats[[m]], v, ps_terms("cov", m), px_opts$n_perm_div, px_opts$n_expect)) })
      })
      for (k in names(V)) {
        M <- do.call(rbind, lapply(res, `[[`, k)); if (is.null(M) || !nrow(M)) next
        q <- function(v) c(mean(M[, v], na.rm = TRUE), stats::quantile(M[, v], c(.025, .975), na.rm = TRUE))
        X1[[length(X1) + 1]] <- data.frame(modality = m, diversity = p, variable = k, n = n_use,
          rho_vi_gpc1 = mean(M[, "rho_vi_gpc1"]),
          sim_unadj = q("sim_unadj")[1], sim_unadj_lo = q("sim_unadj")[2], sim_unadj_hi = q("sim_unadj")[3],
          sim_expect = mean(M[, "sim_expect"]),
          sim_adj = q("sim_adj")[1], sim_adj_lo = q("sim_adj")[2], sim_adj_hi = q("sim_adj")[3],
          sim_adj_null_q975 = mean(M[, "sim_adj_null_q975"]),
          pct_resamples_perm_p05 = 100 * mean(M[, "sim_adj_perm_p"] < .05),
          z_gadj = mean(M[, "z_gadj"]), z_vi = mean(M[, "z_vi"]), n_boot = nrow(M), row.names = NULL)
      }
      cat(sprintf("    %s diversity %.2f done\n", m, p))
    }
  }
  write(do.call(rbind, X1), "px_X1_map_similarity_by_diversity.csv")

  ## ------------------------------------------------- X1b. full sample
  cat("[PX-1b] Map similarity in the full sample (permutation null, family bootstrap)\n")
  X1b <- list()
  for (m in mods) {
    dm <- sub_mod(d, m)
    for (k in names(V)) {
      ok <- rows_ok(dm, feats[[m]], V[[k]], list(c(ps_terms("cov", m), "gpc01"))); s <- dm[ok, ]
      if (nrow(s) < 500) next
      st <- px_map_stats(s, feats[[m]], V[[k]], ps_terms("cov", m), px_opts$n_perm_full, px_opts$n_expect)
      byf <- split(seq_len(nrow(s)), s$fam); fams <- names(byf)
      bs <- t(replicate(px_opts$n_cboot, {
        ix <- unlist(byf[sample(fams, replace = TRUE)], use.names = FALSE)
        x <- mr_safe("boot", px_map_stats(s[ix, ], feats[[m]], V[[k]], ps_terms("cov", m), 0, 1))
        if (is.null(x)) c(NA, NA) else x[c("sim_unadj", "sim_adj")] }))
      ## normal interval from the bootstrap SD: resampled map correlations are
      ## biased towards 0 (more estimation noise), so percentile and basic
      ## intervals can exclude the estimate
      bi <- function(est, v) est + c(-1, 1) * 1.96 * stats::sd(v, na.rm = TRUE)
      X1b[[length(X1b) + 1]] <- data.frame(modality = m, variable = k, n = nrow(s), n_families = length(fams),
        t(st), sim_unadj_lo = bi(st[["sim_unadj"]], bs[, 1])[1], sim_unadj_hi = bi(st[["sim_unadj"]], bs[, 1])[2],
        sim_adj_lo = bi(st[["sim_adj"]], bs[, 2])[1], sim_adj_hi = bi(st[["sim_adj"]], bs[, 2])[2],
        n_perm = px_opts$n_perm_full, n_boot = px_opts$n_cboot, row.names = NULL)
    }
    cat(sprintf("    %s done\n", m))
  }
  write(do.call(rbind, X1b), "px_X1_map_similarity_full.csv")

  ## ------------------------------------------------- X2. head and body size
  cat("[PX-2] Head and body size\n")
  anthro <- c("icv", V[["height"]], V[["weight"]]); anthro <- anthro[!is.na(anthro)]
  have_icv <- any(!is.na(d$icv))
  if (!have_icv) cat("    intracranial volume not available: ICV terms dropped\n")
  if (!have_icv) anthro <- setdiff(anthro, "icv")
  X2a <- list(); X2b <- list(); X2c <- list()
  for (m in mods) {
    dm <- sub_mod(d, m); tc <- ps_terms("cov", m)
    ok <- rows_ok(dm, feats[[m]], NULL, list(c(tc, gpc_names(), anthro))); s <- dm[ok, ]
    if (nrow(s) < 500) next
    sets <- if (m == "smri") split(feats[[m]], pd$feat_type[feats[[m]]]) else stats::setNames(list(feats[[m]]), m)
    for (ty in names(sets)) {
      f <- sets[[ty]]
      g0 <- dR2(s, f, tc, gpc_names())
      gi <- if (have_icv) dR2(s, f, c(tc, "icv"), gpc_names()) else rep(NA, length(f))
      ga <- dR2(s, f, c(tc, anthro), gpc_names())
      an <- dR2(s, f, tc, anthro)
      X2a[[length(X2a) + 1]] <- data.frame(modality = m, measure_type = ty, n = nrow(s), n_features = length(f),
        null_expected = length(gpc_names()) / nrow(s),
        gpc_dR2_mean = mean(g0), gpc_dR2_p99 = abs_q(g0, .99),
        gpc_dR2_beyond_icv_mean = mean(gi), gpc_dR2_beyond_icv_height_weight_mean = mean(ga),
        gpc_dR2_beyond_icv_height_weight_p99 = abs_q(ga, .99),
        pct_gpc_dR2_remaining = 100 * mean(ga) / mean(g0),
        anthro_dR2_mean = mean(an))
      for (k in intersect(c("ses", "adi", "iq", "fiq"), names(V))) {
        ok2 <- !is.na(s[[V[[k]]]]); s2 <- s[ok2, ]
        X2b[[length(X2b) + 1]] <- data.frame(modality = m, measure_type = ty, variable = k, n = nrow(s2),
          pct_red_gpc = px_red(s2, f, V[[k]], tc, c(tc, gpc_names())),
          pct_red_anthro = px_red(s2, f, V[[k]], tc, c(tc, anthro)),
          pct_red_gpc_given_anthro = px_red(s2, f, V[[k]], c(tc, anthro), c(tc, anthro, gpc_names())),
          pct_red_anthro_given_gpc = px_red(s2, f, V[[k]], c(tc, gpc_names()), c(tc, gpc_names(), anthro)))
      }
      ## genetically similar subsamples by measure type (COI and IQ)
      for (k in intersect(c("ses", "iq"), names(V))) {
        s2 <- s[!is.na(s[[V[[k]]]]), ]
        for (gname in c("EUR", "AFR")) {
          sel <- if (gname == "EUR") s2$anc_EUR >= ps_opts$eur_cut else s2$anc_AFR >= ps_opts$afr_cut
          sel <- sel %in% TRUE; if (sum(sel) < ps_opts$min_sub) next
          ag <- strength(assoc_t(s2[sel, ], f, V[[k]], tc))
          am <- lapply(seq_len(px_opts$n_match_boot), function(b) assoc_t(s2[sample.int(nrow(s2), sum(sel)), ], f, V[[k]], tc))
          mt <- vapply(am, strength, 0)
          ## a ratio over a reference that is itself noise is meaningless: require the
          ## size-matched reference to exceed 3 null SDs (as pct_red does)
          if (mean(mt) <= 3 * strength_null_sd(am[[1]])) mt[] <- NA
          X2c[[length(X2c) + 1]] <- data.frame(modality = m, measure_type = ty, variable = k, group = gname,
            n = sum(sel), strength_ratio = ag / mean(mt),
            strength_ratio_lo = ag / stats::quantile(mt, .975, na.rm = TRUE), strength_ratio_hi = ag / stats::quantile(mt, .025, na.rm = TRUE),
            sd_ratio_vs_full = stats::sd(s2[[V[[k]]]][sel]) / stats::sd(s2[[V[[k]]]]), n_boot = px_opts$n_match_boot,
            row.names = NULL)
        }
      }
    }
    cat(sprintf("    %s done\n", m))
  }
  write(do.call(rbind, X2a), "px_X2_gpc_variance_vs_head_body_size.csv")
  write(do.call(rbind, X2b), "px_X2_attenuation_given_head_body_size.csv")
  if (length(X2c)) write(do.call(rbind, X2c), "px_X2_subsamples_by_measure_type.csv")

  ## ------------------------------------------------- X3. robustness
  cat("[PX-3] Robustness of the population-structure adjustment\n")
  anc <- c("anc_AFR", "anc_AMR", "anc_EAS")   # EUR is the reference (proportions sum to 1)
  have_anc <- all(anc %in% names(d)) && any(!is.na(d$anc_AFR))
  X3 <- list()
  for (m in mods) {
    dm <- sub_mod(d, m)
    tc <- ps_terms("cov", m)
    tc_ss <- c(setdiff(tc, c("site", "manuf")), "scanner_serial")
    variants <- list(
      "PCs 1-5" = c(tc, gpc_names()[1:5]),
      "PCs 1-10" = c(tc, gpc_names()),
      "all available PCs" = c(tc, gall),
      "ADMIXTURE proportions" = if (have_anc) c(tc, anc),
      "PCs 1-10 + ADMIXTURE" = if (have_anc) c(tc, gpc_names(), anc),
      "site x scanner, no PCs" = tc_ss,
      "site x scanner + PCs 1-10" = c(tc_ss, gpc_names()))
    variants <- variants[!vapply(variants, is.null, TRUE)]
    for (k in intersect(px_opts$robust_vis, names(V))) {
      need <- unique(c(unlist(variants), fam_ses))
      ok <- rows_ok(dm, feats[[m]], V[[k]], list(need)); s <- dm[ok, ]
      if (nrow(s) < 500) next
      a0 <- assoc_t(s, feats[[m]], V[[k]], tc); af <- assoc_t(s, feats[[m]], V[[k]], c(tc, fam_ses))
      for (vn in names(variants)) {
        tv <- variants[[vn]]; add <- setdiff(tv, tc)
        ref <- if (grepl("site x scanner", vn)) tc_ss else tc
        a1 <- assoc_t(s, feats[[m]], V[[k]], tv)
        a1f <- assoc_t(s, feats[[m]], V[[k]], c(tv, fam_ses))
        R2 <- if (length(setdiff(tv, ref))) px_R2(s, V[[k]], ref, setdiff(tv, ref)) else 0
        aref <- if (identical(ref, tc)) a0 else assoc_t(s, feats[[m]], V[[k]], ref)
        X3[[length(X3) + 1]] <- data.frame(modality = m, variable = k, adjustment = vn, reference = if (identical(ref, tc)) "covariates" else "covariates with site x scanner",
          n = nrow(s), n_terms_added = length(setdiff(tv, ref)),
          pct_red_vs_cov = pct_red(a0, a1), pct_red_vs_reference = pct_red(aref, a1), pct_red_given_fam = pct_red(af, a1f),
          partialR2_vi_on_added = R2, expected_pct_red = 100 * (1 - sqrt(1 - R2)),
          r_map_vs_cov = stats::cor(a0$r, a1$r))
      }
    }
    cat(sprintf("    %s done\n", m))
  }
  X3 <- do.call(rbind, X3)
  X3$n_pcs_available <- length(gall)
  write(X3, "px_X3_robustness.csv")
  cat("\nDone. Aggregate outputs in ", out_dir, "\n")
  invisible(TRUE)
}

########################################################################
##  Figures: 7.0 versions of draft Figs 3-5 from aggregate tables.
##  main_dir: output/popstrat_main (ps_A_*.csv); ext_dir: popstrat_ext.
########################################################################
px_figures <- function(main_dir, ext_dir, file = file.path(ext_dir, "popstrat_7p0_figures.pdf")) {
  rd <- function(f) { p <- file.path(f); if (file.exists(p)) utils::read.csv(p) else NULL }
  A1 <- rd(file.path(main_dir, "ps_A_gpc_dR2_by_diversity.csv"))
  A2 <- rd(file.path(main_dir, "ps_A_effect_reduction_by_diversity.csv"))
  X1 <- rd(file.path(ext_dir, "px_X1_map_similarity_by_diversity.csv"))
  mods <- c(smri = "Structural T1", dti = "Diffusion", rs = "Resting state", nback = "Task (n-back)")
  vlab <- c(ses = "COI", adi = "ADI", iq = "IQ", fiq = "Fluid IQ", internalizing = "Internalizing",
            height = "Height", weight = "Weight")
  pal <- c(ses = "#1b6ca8", adi = "#5aa9e6", iq = "#c0392b", fiq = "#e67e73", internalizing = "#7d3c98",
           height = "#2e8b57", weight = "#9acd32")
  grDevices::pdf(file, width = 11, height = 3.4)
  arrows <- function(...) suppressWarnings(graphics::arrows(...))
  on.exit(grDevices::dev.off())
  op <- graphics::par(mfrow = c(1, 4), mar = c(4, 4, 2.5, 0.5), oma = c(0, 0, 2, 8), xpd = NA)
  ## Fig 3: variance in IDPs explained by PCs 1-10
  if (!is.null(A1)) {
    yl <- c(0, max(A1$p99, na.rm = TRUE) * 1.05)
    for (m in names(mods)) { x <- A1[A1$modality == m, ]; x <- x[order(x$diversity), ]
      graphics::plot(NA, xlim = c(0, 1), ylim = yl, xlab = "Proportion not self-identified White",
                     ylab = if (m == "smri") "Variance explained by PCs 1-10" else "", main = mods[[m]])
      if (!nrow(x)) next
      graphics::lines(x$diversity, x$p50, lwd = 2, col = "grey60"); graphics::lines(x$diversity, x$p75, lwd = 2, col = "grey35")
      graphics::lines(x$diversity, x$p99, lwd = 2, col = "black")
      if (all(c("dR2_p99_lo", "dR2_p99_hi") %in% names(x))) arrows(x$diversity, x$dR2_p99_lo, x$diversity, x$dR2_p99_hi, angle = 90, code = 3, length = .03)
    }
    graphics::legend("topright", inset = c(-0.55, 0), legend = c("99th pct", "75th pct", "median"),
                     col = c("black", "grey35", "grey60"), lwd = 2, bty = "n")
    graphics::mtext("Variance in imaging measures explained by genetic PCs, by sample diversity (ABCD 7.0)", outer = TRUE, font = 2)
  }
  ## Fig 4: % reduction in association strength from adding PCs
  if (!is.null(A2)) {
    for (m in names(mods)) { x <- A2[A2$modality == m, ]
      graphics::plot(NA, xlim = c(0, 1), ylim = c(-50, 150), xlab = "Proportion not self-identified White",
                     ylab = if (m == "smri") "% reduction in association strength" else "", main = mods[[m]])
      graphics::abline(h = c(0, 100), lty = 3, col = "grey60", xpd = FALSE)
      for (k in names(vlab)) { z <- x[x$variable == k, ]; z <- z[order(z$diversity), ]; if (!nrow(z)) next
        graphics::lines(z$diversity, z$pct_reduction_mean, col = pal[[k]], lwd = 2)
        arrows(z$diversity, z$pct_reduction_lo, z$diversity, z$pct_reduction_hi, col = pal[[k]], angle = 90, code = 3, length = .02) }
    }
    graphics::legend("topright", inset = c(-0.6, 0), legend = vlab, col = pal, lwd = 2, bty = "n")
    graphics::mtext("Reduction in association strength after adding genetic PCs 1-10 (ABCD 7.0)", outer = TRUE, font = 2)
  }
  ## Fig 5: map similarity with the PC1 map, unadjusted vs adjusted
  if (!is.null(X1)) for (what in c("sim_unadj", "sim_adj")) {
    for (m in names(mods)) { x <- X1[X1$modality == m, ]
      graphics::plot(NA, xlim = c(0, 1), ylim = c(-1, 1), xlab = "Proportion not self-identified White",
                     ylab = if (m == "smri") "Map correlation" else "", main = mods[[m]])
      graphics::abline(h = 0, lty = 3, col = "grey60", xpd = FALSE)
      for (k in names(vlab)) { z <- x[x$variable == k, ]; z <- z[order(z$diversity), ]; if (!nrow(z)) next
        graphics::lines(z$diversity, z[[what]], col = pal[[k]], lwd = 2)
        arrows(z$diversity, z[[paste0(what, "_lo")]], z$diversity, z[[paste0(what, "_hi")]], col = pal[[k]], angle = 90, code = 3, length = .02)
        if (what == "sim_unadj") graphics::lines(z$diversity, z$sim_expect, col = pal[[k]], lty = 2)
      }
      if (what == "sim_adj" && nrow(x)) { nb <- stats::aggregate(sim_adj_null_q975 ~ diversity, x, mean)
        graphics::polygon(c(nb$diversity, rev(nb$diversity)), c(nb$sim_adj_null_q975, -rev(nb$sim_adj_null_q975)),
                          col = grDevices::adjustcolor("grey70", .35), border = NA) }
    }
    graphics::legend("topright", inset = c(-0.6, 0), legend = vlab, col = pal, lwd = 2, bty = "n")
    graphics::mtext(if (what == "sim_unadj") "Similarity of each variable's map to the PC1 map (dashed: expected if PC1 acts only through the variable)"
                    else "Similarity of each variable's map to the map of PC1 adjusted for the variable (grey band: 95% permutation null)",
                    outer = TRUE, font = 2, cex = .85)
  }
  graphics::par(op)
  invisible(file)
}
