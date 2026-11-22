########################################################################
##  sc_analysis.R — causal-inference ladder for Marek et al.'s exposure
##  (neighborhood SES, COI social and economic domain) on the brain and
##  on IQ. Mirrors the music-ed pipeline's methods where they are
##  identified for a neighborhood exposure; aggregate output only.
##
##  Rung (pipeline method number)            identifies the effect from
##  1  LME, pooled waves (M1)                 all variation, adjusted
##  3  AIPW, top vs bottom COI tertile (M3)   baseline, adjusted (DR)
##  4  Causal forest ATE (M4, if grf)         baseline, adjusted (ML-DR)
##  5  CLPM (M5)                              lagged waves
##  6  RI-CLPM (M6)                           within-person lagged change
##  7  Person fixed effects (M7)              within-child change = moves
##  7b First differences + placebo            moves; future move as placebo
##  8  Family FE (M8)                         reported only if COI varies
##                                            within family (it rarely does)
##  Each selection-on-observables rung is run under four adjustment sets:
##  cov (age, sex, site, scanner, motion), fam (+ baseline income,
##  caregiver education), gpc (+ genetic PCs 1-10), famgpc (+ both).
##  Person-FE rungs absorb every time-invariant confounder, including
##  genetic ancestry and baseline family SES.
##  Outcomes and exposure are in SD units (baseline SDs): estimates are
##  SD of outcome per SD of COI (AIPW: per top-vs-bottom tertile contrast).
########################################################################

sc_opts <- modifyList(list(
  seed = 20261121, n_boot = 200, outcomes = c("sig_rs", "sig_ct", "ms_rs", "iq", "fiq", "height"),
  use_grf = TRUE, move_eps = 0.01
), if (exists("sc_opts_user")) sc_opts_user else list())

sc_terms <- function(adj, outcome) {
  base <- c("age", "sex", "site", "manuf", if (outcome %in% c("sig_rs", "ms_rs")) "motion")
  switch(adj, cov = base, fam = c(base, "income", "edu"), gpc = c(base, gpc_names()),
         famgpc = c(base, "income", "edu", gpc_names()))
}

## cluster-robust (CR1) coefficient table for an lm.fit
cr_se <- function(X, e, cl) {
  XtXi <- solve(crossprod(X)); G <- length(unique(cl)); n <- nrow(X); k <- ncol(X)
  U <- rowsum(X * e, cl); meat <- crossprod(U)
  V <- XtXi %*% meat %*% XtXi * (G / (G - 1)) * ((n - 1) / (n - k))
  sqrt(pmax(diag(V), 0))
}
row_out <- function(method, outcome, adj, est, se, n, n_id, n_fam, note = "") {
  data.frame(method = method, outcome = outcome, adjustment = adj, estimate = est, se = se,
             lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * stats::pnorm(-abs(est / se)),
             n = n, n_id = n_id, n_fam = n_fam, note = note)
}

sc_run <- function(sp, out_dir, write = function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE) }) {
  set.seed(sc_opts$seed); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  P <- sp$panel; res <- list()
  ## the data-derived signatures are analysed in the replication half only
  sample_for <- function(y) if (y %in% c("sig_rs", "sig_ct")) P[P$split %in% "rep", ] else P
  ## standardize on baseline SDs
  z_by_base <- function(d, v) { s <- stats::sd(d[[v]][d$wave == 0], na.rm = TRUE); m <- mean(d[[v]][d$wave == 0], na.rm = TRUE); (d[[v]] - m) / s }
  outs <- intersect(sc_opts$outcomes, names(P))

  ## ------------------------------------------------------------ sample
  smp <- do.call(rbind, lapply(outs, function(y) {
    d <- sample_for(y)
    do.call(rbind, lapply(sort(unique(d$wave)), function(w) {
      x <- d[d$wave == w, ]
      data.frame(outcome = y, wave = w, term = c("outcome and COI", "COI", "outcome"),
                 n = c(sum(!is.na(x[[y]]) & !is.na(x$coi)), sum(!is.na(x$coi)), sum(!is.na(x[[y]]))))
    }))
  }))
  ## movers: COI changes between consecutive observed waves
  mv <- P[!is.na(P$coi), c("id", "wave", "coi")]; mv <- mv[order(mv$id, mv$wave), ]
  dco <- ave(mv$coi, mv$id, FUN = function(v) c(NA, diff(v)))
  movers <- data.frame(term = c("children with COI at >= 2 waves", "children whose COI changes (movers)",
                                "SD of COI change among movers (national z units)"),
                       value = c(length(unique(mv$id[duplicated(mv$id)])),
                                 length(unique(mv$id[abs(dco) > sc_opts$move_eps & !is.na(dco)])),
                                 stats::sd(dco[abs(dco) > sc_opts$move_eps], na.rm = TRUE)),
                       n = length(unique(mv$id)))
  write(smp, "sc_sample.csv"); write(movers, "sc_movers.csv")
  ## In ABCD 7.0 every linked-external (le_*) variable, COI and ADI included,
  ## exists at ses-00A only, so no child has COI at two waves. The
  ## within-child rungs (5, 6, 7, 7b, 15) then have no exposure variation
  ## and are reported as not identified instead of being attempted.
  x_long <- movers$value[1] > 0
  why_not <- "not identified: COI measured at baseline only (ABCD 7.0 linked-external data are ses-00A only)"
  lme_lab <- if (x_long) "1 LME (pooled waves)" else "1 LME (baseline; COI is baseline-only)"

  for (y in outs) {
    cat(sprintf("\n[SC] outcome: %s\n", y))
    d <- sample_for(y); d$Y <- z_by_base(d, y); d$X <- z_by_base(d, "coi")
    d <- d[!is.na(d$Y) & !is.na(d$X) & !is.na(d$fam), ]
    if (nrow(d) < 500) { cat("  too few observations; skipped\n"); next }

    ## ---- 1. LME, pooled waves, random intercepts for family and child
    for (adj in c("cov", "fam", "gpc", "famgpc")) {
      tt <- sc_terms(adj, y); dd <- d[stats::complete.cases(d[, c("Y", "X", tt)]), ]
      ## manufacturer is (nearly) nested in site: one combined factor avoids the
      ## rank deficiency nlme cannot handle
      dd$scanner <- droplevels(interaction(dd$site, dd$manuf, drop = TRUE))
      tt <- c(setdiff(tt, c("site", "manuf")), "scanner")
      if (length(unique(dd$wave)) > 1) tt <- c(tt, "factor(wave)")
      f <- stats::as.formula(paste("Y ~ X +", paste(tt, collapse = " + ")))
      fit <- mr_safe("lme", nlme::lme(f, random = ~ 1 | fam / id, data = dd, na.action = stats::na.omit,
                                       control = nlme::lmeControl(opt = "optim", returnObject = TRUE)))
      if (!is.null(fit)) { tb <- summary(fit)$tTable
        res[[length(res) + 1]] <- row_out(lme_lab, y, adj, tb["X", 1], tb["X", 2],
                                          nrow(dd), length(unique(dd$id)), length(unique(dd$fam))) }
    }

    ## ---- 3/4. AIPW and causal forest at baseline: top vs bottom COI tertile
    b0 <- d[d$wave == 0, ]
    if (nrow(b0) >= 500) {
      q <- stats::quantile(b0$X, c(1/3, 2/3), na.rm = TRUE)
      b0$A <- ifelse(b0$X >= q[2], 1, ifelse(b0$X <= q[1], 0, NA)); b0 <- b0[!is.na(b0$A), ]
      for (adj in c("cov", "fam", "gpc", "famgpc")) {
        tt <- sc_terms(adj, y); bb <- b0[stats::complete.cases(b0[, c("Y", "A", tt)]), ]
        if (nrow(bb) < 300) next
        Z <- make_Z(bb, tt)[, -1, drop = FALSE]
        aipw <- function(ix) {
          bx <- bb[ix, ]; Zx <- Z[ix, , drop = FALSE]
          ps <- stats::plogis(stats::glm.fit(cbind(1, Zx), bx$A, family = stats::binomial())$linear.predictors)
          ps <- pmin(pmax(ps, 0.01), 0.99)
          m1 <- stats::lm.fit(cbind(1, Zx[bx$A == 1, , drop = FALSE]), bx$Y[bx$A == 1])$coefficients
          m0 <- stats::lm.fit(cbind(1, Zx[bx$A == 0, , drop = FALSE]), bx$Y[bx$A == 0])$coefficients
          m1[is.na(m1)] <- 0; m0[is.na(m0)] <- 0
          mu1 <- cbind(1, Zx) %*% m1; mu0 <- cbind(1, Zx) %*% m0
          mean(mu1 - mu0 + bx$A * (bx$Y - mu1) / ps - (1 - bx$A) * (bx$Y - mu0) / (1 - ps))
        }
        ps_all <- mr_safe("ps", stats::plogis(stats::glm.fit(cbind(1, Z), bb$A, family = stats::binomial())$linear.predictors))
        est <- mr_safe("aipw", aipw(seq_len(nrow(bb))))
        if (is.null(est)) next
        fams <- unique(bb$fam); byf <- split(seq_len(nrow(bb)), bb$fam)
        bs <- replicate(sc_opts$n_boot, { g <- sample(fams, replace = TRUE)
          ix <- unlist(byf[as.character(g)], use.names = FALSE); nz(mr_safe("aipw boot", aipw(ix))) })
        res[[length(res) + 1]] <- row_out("3 AIPW (baseline, top vs bottom tertile)", y, adj, est,
                                          stats::sd(bs, na.rm = TRUE), nrow(bb), nrow(bb), length(fams),
                                          sprintf("PS in [.05,.95]: %.1f%%", if (is.null(ps_all)) NA else 100 * mean(ps_all > .05 & ps_all < .95)))
        if (isTRUE(sc_opts$use_grf) && requireNamespace("grf", quietly = TRUE) && adj %in% c("cov", "famgpc")) {
          cf <- mr_safe("grf", grf::causal_forest(Z, bb$Y, bb$A, clusters = as.integer(factor(bb$fam)), num.trees = 1000))
          if (!is.null(cf)) { a <- grf::average_treatment_effect(cf, target.sample = "overlap")
            res[[length(res) + 1]] <- row_out("4 Causal forest (baseline, overlap ATE)", y, adj, a[[1]], a[[2]],
                                              nrow(bb), nrow(bb), length(fams)) }
        }
      }
    }

    ## residualize Y within wave on age/scanner/motion for the panel models
    tt0 <- setdiff(sc_terms("cov", y), c("sex"))
    d$Yr <- NA_real_
    for (w in unique(d$wave)) { i <- which(d$wave == w & stats::complete.cases(d[, tt0]))
      if (length(i) > 50) d$Yr[i] <- as.numeric(resid_on(d$Y[i], make_Z(d[i, ], tt0))) }

    if (!x_long) {
      for (rg in c("5 CLPM", "6 RI-CLPM", "7 Person FE (within-child)", "7b First difference: COI change -> outcome change"))
        res[[length(res) + 1]] <- row_out(rg, y, "within-child", NA, NA, NA, NA, NA, why_not)
    } else {
    ## ---- 5/6. CLPM and RI-CLPM over the available waves
    ws <- sort(unique(d$wave[!is.na(d$Yr)]))
    if (length(ws) >= 3 && requireNamespace("lavaan", quietly = TRUE)) {
      wide <- stats::reshape(d[!is.na(d$Yr), c("id", "fam", "wave", "X", "Yr")], idvar = c("id", "fam"),
                             timevar = "wave", direction = "wide")
      xs <- paste0("X.", ws); ys <- paste0("Yr.", ws)
      clpm <- c(vapply(2:length(ws), function(t) sprintf("%s ~ ax*%s + byx*%s", xs[t], xs[t - 1], ys[t - 1]), ""),
                vapply(2:length(ws), function(t) sprintf("%s ~ ay*%s + bxy*%s", ys[t], ys[t - 1], xs[t - 1]), ""),
                vapply(seq_along(ws), function(t) sprintf("%s ~~ %s", xs[t], ys[t]), ""))
      ri <- c(sprintf("RIx =~ %s", paste0("1*", xs, collapse = " + ")), sprintf("RIy =~ %s", paste0("1*", ys, collapse = " + ")),
              vapply(seq_along(ws), function(t) sprintf("wx%d =~ 1*%s; wy%d =~ 1*%s; %s ~~ 0*%s; %s ~~ 0*%s",
                                                       t, xs[t], t, ys[t], xs[t], xs[t], ys[t], ys[t]), ""),
              vapply(2:length(ws), function(t) sprintf("wx%d ~ ax*wx%d + byx*wy%d; wy%d ~ ay*wy%d + bxy*wx%d",
                                                      t, t - 1, t - 1, t, t - 1, t - 1), ""),
              vapply(seq_along(ws), function(t) sprintf("wx%d ~~ wy%d", t, t), ""),
              "RIx ~~ RIy", "RIx ~~ 0*wx1 + 0*wy1", "RIy ~~ 0*wx1 + 0*wy1")
      for (spec in list(list("5 CLPM", clpm), list("6 RI-CLPM", ri))) {
        fit <- mr_safe(spec[[1]], lavaan::sem(paste(spec[[2]], collapse = "\n"), data = wide, estimator = "MLR",
                                               missing = "fiml", cluster = "fam"))
        if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) {
          cat(sprintf("  %s did not converge\n", spec[[1]])); next }
        pe <- lavaan::parameterEstimates(fit)
        for (lab in c("bxy", "byx")) { r <- pe[pe$label == lab, ][1, ]
          res[[length(res) + 1]] <- row_out(paste(spec[[1]], if (lab == "bxy") "COI -> outcome (2-yr lag)" else "outcome -> COI (falsification)"),
                                            y, "cov (outcome residualized)", r$est, r$se, nrow(wide), nrow(wide),
                                            length(unique(wide$fam)), sprintf("waves %s", paste(ws, collapse = ","))) }
      }
    }

    ## ---- 7. person fixed effects (within-child change in COI = moves)
    dp <- d
    tv <- intersect(c("age", "motion"), sc_terms("cov", y))
    dp <- dp[stats::complete.cases(dp[, c("Y", "X", tv, "manuf")]), ]
    dp <- dp[dp$id %in% dp$id[duplicated(dp$id)], ]
    if (nrow(dp) > 500) {
      Xm <- cbind(X = dp$X, as.matrix(dp[, tv, drop = FALSE]), stats::model.matrix(~ factor(wave) + manuf, dp)[, -1, drop = FALSE])
      dm <- function(M) M - apply(M, 2, function(v) stats::ave(v, dp$id))
      Xw <- dm(Xm); yw <- dm(matrix(dp$Y))[, 1]
      keep <- colSums(abs(Xw)) > 1e-8; Xw <- Xw[, keep, drop = FALSE]
      n_mov <- sum(tapply(dp$X, dp$id, function(v) diff(range(v)) > sc_opts$move_eps))
      if (!"X" %in% colnames(Xw) || n_mov < 30) {
        res[[length(res) + 1]] <- row_out("7 Person FE (within-child)", y, "all time-invariant confounders", NA, NA,
                                          nrow(dp), length(unique(dp$id)), length(unique(dp$fam)),
                                          sprintf("not identified: %d movers", n_mov))
      } else {
      fit <- stats::lm.fit(Xw, yw); e <- fit$residuals
      se <- cr_se(Xw, e, dp$fam); j <- which(colnames(Xw) == "X")
      res[[length(res) + 1]] <- row_out("7 Person FE (within-child)", y, "all time-invariant confounders",
                                        fit$coefficients[j], se[j], nrow(dp), length(unique(dp$id)),
                                        length(unique(dp$fam)), sprintf("movers contributing: %d", n_mov))
      }
    }

    ## ---- 7b. first differences with a future-move placebo
    fd <- d[order(d$id, d$wave), c("id", "fam", "wave", "Y", "X", tv)]
    lagv <- function(v) stats::ave(v, fd$id, FUN = function(z) c(NA, z[-length(z)]))
    fd$dY <- fd$Y - lagv(fd$Y); fd$dX <- fd$X - lagv(fd$X)
    for (v in tv) fd[[paste0("d", v)]] <- fd[[v]] - lagv(fd[[v]])
    fd$dX_next <- stats::ave(fd$dX, fd$id, FUN = function(z) c(z[-1], NA))
    dv <- paste0("d", tv)
    f1 <- fd[stats::complete.cases(fd[, c("dY", "dX", dv)]), ]
    if (nrow(f1) > 300 && sum(abs(f1$dX) > sc_opts$move_eps) >= 30) {
      M <- cbind(1, dX = f1$dX, as.matrix(f1[, dv, drop = FALSE]), stats::model.matrix(~ factor(wave), f1)[, -1, drop = FALSE])
      ft <- stats::lm.fit(M, f1$dY); se <- cr_se(M, ft$residuals, f1$fam)
      res[[length(res) + 1]] <- row_out("7b First difference: COI change -> outcome change", y, "time-invariant absorbed",
                                        ft$coefficients[2], se[2], nrow(f1), length(unique(f1$id)), length(unique(f1$fam)),
                                        sprintf("intervals with a move: %d", sum(abs(f1$dX) > sc_opts$move_eps)))
      f2 <- fd[stats::complete.cases(fd[, c("dY", "dX", "dX_next", dv)]), ]
      if (nrow(f2) > 300 && sum(abs(f2$dX_next) > sc_opts$move_eps) >= 30) {
        M <- cbind(1, dX_next = f2$dX_next, dX = f2$dX, as.matrix(f2[, dv, drop = FALSE]), stats::model.matrix(~ factor(wave), f2)[, -1, drop = FALSE])
        ft <- stats::lm.fit(M, f2$dY); se <- cr_se(M, ft$residuals, f2$fam)
        res[[length(res) + 1]] <- row_out("15 Placebo: FUTURE COI change -> current outcome change", y, "time-invariant absorbed",
                                          ft$coefficients[2], se[2], nrow(f2), length(unique(f2$id)), length(unique(f2$fam)),
                                          "should be ~0 if moves are not anticipated by outcome trends")
      }
    }

    }

    ## ---- 8. family FE at baseline: identified only if COI varies within family
    b0 <- d[d$wave == 0, ]; b0 <- b0[b0$fam %in% b0$fam[duplicated(b0$fam)], ]
    if (nrow(b0) > 100) {
      wX <- b0$X - stats::ave(b0$X, b0$fam); wsd <- stats::sd(wX)
      n_var <- sum(tapply(b0$X, b0$fam, function(v) diff(range(v)) > sc_opts$move_eps))
      if (n_var < 30) {
        res[[length(res) + 1]] <- row_out("8 Family FE (baseline)", y, "shared family factors", NA, NA, nrow(b0), nrow(b0),
                                          length(unique(b0$fam)), sprintf("not identified: COI differs within only %d families (within SD %.3f)", n_var, wsd))
      } else {
        tvf <- intersect(c("age", "motion"), sc_terms("cov", y)); bf <- b0[stats::complete.cases(b0[, c(tvf, "sex"), drop = FALSE]), ]
        Mf <- cbind(X = bf$X, as.matrix(bf[, tvf, drop = FALSE]), sexn = as.numeric(bf$sex))
        Mf <- Mf - apply(Mf, 2, function(v) stats::ave(v, bf$fam)); yf <- bf$Y - stats::ave(bf$Y, bf$fam)
        Mf <- Mf[, colSums(abs(Mf)) > 1e-8, drop = FALSE]
        ft <- stats::lm.fit(Mf, yf); se <- cr_se(Mf, ft$residuals, bf$fam); j <- which(colnames(Mf) == "X")
        res[[length(res) + 1]] <- row_out("8 Family FE (baseline)", y, "shared family factors", ft$coefficients[j], se[j],
                                          nrow(bf), nrow(bf), length(unique(bf$fam)),
                                          sprintf("COI differs within %d families (siblings at different addresses)", n_var))
      }
    }
  }
  ladder <- do.call(rbind, res)
  write(ladder, "sc_ladder.csv")
  cat("\nDone. Aggregate outputs in ", out_dir, "\n")
  invisible(ladder)
}
