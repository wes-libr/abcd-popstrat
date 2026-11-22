########################################################################
##  mr_analysis.R — analysis functions for the Marek et al. (2026,
##  Science 392:eaee6213) replication on ABCD 7.0, with the
##  population-structure (genetic PC) extension from the popstrat paper.
##
##  Pure functions on the cached analytic object `md` built by mr_data.R.
##  Nothing here reads or writes participant-level data; mr_run() writes
##  aggregate tables only (correlations, out-of-sample r, counts >= 10).
##
##  Adjustment sets used throughout:
##    raw : bivariate correlation (Marek's univariate BWAS)
##    cov : age, sex, site, scanner manufacturer (+ mean FD for RSFC)
##    gpc : cov + genetic PCs 1-10 (gn_y_popstruct_pc__01..10)
##  Extra adjustors (SES when the target is IQ, IQ when the target is SES)
##  are added on top of cov or gpc where noted.
########################################################################

mr_opts <- modifyList(list(
  seed          = 20261116,
  min_n         = 2000,     # Marek: continuous variables need n >= 2000
  min_cat       = 2000,     # Marek: >= 2000 per category for binary variables
  max_items_per_table = Inf,
  n_gpc         = 10,
  lambdas       = 10^seq(-1, 5, length.out = 31),
  folds         = 5,
  n_boot_map    = 200,      # participant bootstrap for key map correlations
  n_split_half  = 20,       # split-half resamples for the IQ/SES ridge models
  n_boot_strat  = 100,      # SES-stratified training resamples
  n_boot_grad   = 50,       # resamples per step of the SES-composition gradient
  high_q        = 0.75,     # high-SES group: top quartile of discovery SES
  low_q         = 0.25,     # low-SES group: bottom quartile
  iq_bins       = 10,       # IQ deciles for distribution matching
  do_ridge_all  = TRUE,     # ridge for every variable (slowest step)
  do_year2      = TRUE,
  ms_nets       = c("smh", "smm", "vis", "aud"),   # motor + sensory Gordon networks
  modalities    = c("rs", "ct")
), if (exists("mr_opts_user")) mr_opts_user else list())

gpc_names <- function() sprintf("gpc%02d", seq_len(mr_opts$n_gpc))

## ---- small helpers ----------------------------------------------------------
mr_safe <- function(label, expr) tryCatch(expr, error = function(e) {
  cat(sprintf("  [%s skipped: %s]\n", label, substr(conditionMessage(e), 1, 200))); NULL })

abs_q <- function(r, p) as.numeric(stats::quantile(abs(r), p, na.rm = TRUE, names = FALSE))

## Covariate model matrix for a data subset. `adj` in raw/cov/gpc; `extra` adds
## named columns of d. Returns NULL rows for incomplete covariates via `ok`.
cov_terms <- function(adj, mod, extra = NULL) {
  t <- character(0)
  if (adj %in% c("cov", "gpc")) t <- c("age", "sex", "site", "manuf", if (mod == "rs") "motion")
  if (adj == "gpc") t <- c(t, gpc_names())
  c(t, extra)
}

cov_ok <- function(d, terms) {
  if (!length(terms)) return(rep(TRUE, nrow(d)))
  stats::complete.cases(d[, terms, drop = FALSE])
}

make_Z <- function(d, terms) {
  if (!length(terms)) return(matrix(1, nrow(d), 1, dimnames = list(NULL, "(Intercept)")))
  dd <- d[, terms, drop = FALSE]
  for (v in names(dd)) if (is.factor(dd[[v]]) || is.character(dd[[v]])) dd[[v]] <- droplevels(factor(dd[[v]]))
  ## a factor with one level in this subset carries no information; drop it
  keep <- vapply(dd, function(x) !(is.factor(x) && nlevels(x) < 2), logical(1))
  dd <- dd[, keep, drop = FALSE]
  if (!ncol(dd)) return(matrix(1, nrow(d), 1, dimnames = list(NULL, "(Intercept)")))
  stats::model.matrix(~ ., data = dd)
}

## residualize the columns of M on Z (least squares, rank-safe)
resid_on <- function(M, Z) qr.resid(qr(Z), as.matrix(M))

## ---- univariate association map ---------------------------------------------
## Correlation of each brain feature with y after residualizing both on the
## adjustment set (partial correlation; raw = Pearson r). Rows: participants of
## `d` that pass the modality QC and have y and the covariates.
assoc_map <- function(d, feats, y, mod, adj = "raw", extra = NULL) {
  terms <- cov_terms(adj, mod, extra)
  ok <- !is.na(d[[y]]) & cov_ok(d, terms) & stats::complete.cases(d[, feats, drop = FALSE])
  if (sum(ok) < 50) return(list(r = stats::setNames(rep(NA_real_, length(feats)), feats), n = sum(ok)))
  dd <- d[ok, , drop = FALSE]
  X <- as.matrix(dd[, feats, drop = FALSE]); yy <- as.numeric(dd[[y]])
  if (stats::sd(yy) == 0) return(list(r = stats::setNames(rep(NA_real_, length(feats)), feats), n = sum(ok)))
  if (length(terms)) {
    Z <- make_Z(dd, terms); X <- resid_on(X, Z); yy <- as.numeric(resid_on(yy, Z))
  }
  r <- suppressWarnings(as.numeric(stats::cor(X, yy)))
  list(r = stats::setNames(r, feats), n = sum(ok))
}

map_summary <- function(r) c(max_abs_r = max(abs(r), na.rm = TRUE), q90_abs_r = abs_q(r, .90),
                             q99_abs_r = abs_q(r, .99), mean_abs_r = mean(abs(r), na.rm = TRUE))

## ---- ridge regression (closed form via SVD, lambda by k-fold CV) ------------
ridge_path_pred <- function(Xtr, ytr, Xte, lambdas) {
  mu <- colMeans(Xtr); s <- apply(Xtr, 2, stats::sd); s[s == 0 | is.na(s)] <- 1
  A <- sweep(sweep(Xtr, 2, mu), 2, s, "/"); B <- sweep(sweep(Xte, 2, mu), 2, s, "/")
  ym <- mean(ytr); sv <- svd(A); uty <- crossprod(sv$u, ytr - ym)
  sapply(lambdas, function(l) ym + B %*% (sv$v %*% ((sv$d / (sv$d^2 + l)) * uty)))
}

ridge_fit_predict <- function(Xtr, ytr, Xte, lambdas = mr_opts$lambdas, folds = mr_opts$folds) {
  n <- nrow(Xtr); f <- sample(rep_len(seq_len(folds), n))
  mse <- matrix(0, folds, length(lambdas))
  for (k in seq_len(folds)) {
    p <- ridge_path_pred(Xtr[f != k, , drop = FALSE], ytr[f != k], Xtr[f == k, , drop = FALSE], lambdas)
    mse[k, ] <- colMeans((p - ytr[f == k])^2)
  }
  lam <- lambdas[which.min(colMeans(mse))]
  as.numeric(ridge_path_pred(Xtr, ytr, Xte, lam))
}

## Residualize X and y on the adjustment set using TRAINING coefficients only,
## then apply them to the test set (no leakage). Returns the out-of-sample r
## and the predictions (for shortcut checks against another target).
prep_resid <- function(dtr, dte, feats, y, terms) {
  Xtr <- as.matrix(dtr[, feats]); Xte <- as.matrix(dte[, feats])
  ytr <- as.numeric(dtr[[y]]); yte <- as.numeric(dte[[y]])
  if (length(terms)) {
    Zall <- make_Z(rbind(dtr[, terms, drop = FALSE], dte[, terms, drop = FALSE]), terms)
    Ztr <- Zall[seq_len(nrow(dtr)), , drop = FALSE]; Zte <- Zall[-seq_len(nrow(dtr)), , drop = FALSE]
    keep <- colSums(Ztr != 0) > 0; Ztr <- Ztr[, keep, drop = FALSE]; Zte <- Zte[, keep, drop = FALSE]
    fitX <- stats::lm.fit(Ztr, cbind(Xtr, ytr)); Bc <- fitX$coefficients; Bc[is.na(Bc)] <- 0
    R <- cbind(Xtr, ytr) - Ztr %*% Bc; Rte <- cbind(Xte, yte) - Zte %*% Bc
    Xtr <- R[, seq_along(feats), drop = FALSE]; ytr <- R[, ncol(R)]
    Xte <- Rte[, seq_along(feats), drop = FALSE]; yte <- Rte[, ncol(Rte)]
  }
  list(Xtr = Xtr, ytr = ytr, Xte = Xte, yte = yte)
}

ridge_oos <- function(dtr, dte, feats, y, mod, adj = "cov", extra = NULL) {
  terms <- cov_terms(adj, mod, extra)
  oktr <- !is.na(dtr[[y]]) & cov_ok(dtr, terms) & stats::complete.cases(dtr[, feats])
  okte <- !is.na(dte[[y]]) & cov_ok(dte, terms) & stats::complete.cases(dte[, feats])
  if (sum(oktr) < 100 || sum(okte) < 100) return(list(r = NA_real_, n_train = sum(oktr), n_test = sum(okte)))
  dtr <- dtr[oktr, , drop = FALSE]; dte <- dte[okte, , drop = FALSE]
  P <- prep_resid(dtr, dte, feats, y, terms)
  if (stats::sd(P$ytr) == 0) return(list(r = NA_real_, n_train = nrow(dtr), n_test = nrow(dte)))
  pred <- ridge_fit_predict(P$Xtr, P$ytr, P$Xte)
  list(r = suppressWarnings(stats::cor(pred, P$yte)), pred = pred, test = dte,
       n_train = nrow(dtr), n_test = nrow(dte))
}

## ---- variable catalogue filter (Marek's n >= 2000 rule) ----------------------
filter_vis <- function(d, vis) {
  keep <- vapply(vis, function(v) {
    x <- d[[v]]; x <- x[!is.na(x)]
    if (length(x) < mr_opts$min_n) return(FALSE)
    u <- unique(x)
    if (length(u) < 2) return(FALSE)
    if (length(u) == 2) return(min(table(x)) >= min(mr_opts$min_cat, length(x) / 4))
    stats::sd(x) > 0 && max(table(x)) / length(x) < 0.95
  }, logical(1))
  vis[keep]
}

## ---- network profile for RSFC maps ------------------------------------------
## Node strength = sum of |r| over network pairs involving the network (a
## within-network pair counts once). Share in motor/sensory = strength of the
## motor + sensory networks / total strength (Marek Fig. 3B analogue).
net_profile <- function(r, net) {
  a <- abs(r[net$feature])
  nets <- sort(unique(c(net$net1, net$net2)))
  st <- vapply(nets, function(k) sum(a[net$net1 == k | net$net2 == k], na.rm = TRUE), numeric(1))
  list(strength = st, ms_share = sum(st[nets %in% mr_opts$ms_nets]) / sum(st))
}

########################################################################
##  mr_run(md, out_dir): every analysis; returns invisibly a list of the
##  tables written. `md` is the object from mr_build() (see mr_data.R).
########################################################################
mr_run <- function(md, out_dir, write = function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE) }) {
  set.seed(mr_opts$seed)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- md$base; disc <- d[d$split == "disc", ]; rep <- d[d$split == "rep", ]
  feats <- md$feat; out <- list()
  K <- md$key   # named vector: ses, iq, pfac, sleep, screen, adi, income, edu, height
  K <- K[K %in% names(d)]
  mods <- intersect(mr_opts$modalities, names(feats))
  sub_mod <- function(x, m) x[x[[paste0("ok_", m)]] %in% TRUE, , drop = FALSE]

  ## ---------------------------------------------------------------- sample
  cat("\n[1] Sample\n")
  smp <- do.call(rbind, lapply(c("disc", "rep"), function(s) {
    x <- d[d$split == s, ]
    data.frame(split = s, n_baseline = nrow(x), n_rs = sum(x$ok_rs %in% TRUE),
               n_ct = sum(x$ok_ct %in% TRUE), n_gpc = sum(stats::complete.cases(x[, gpc_names()])),
               n_ses = sum(!is.na(x[[K[["ses"]]]])), n_iq = sum(!is.na(x[[K[["iq"]]]])))
  }))
  write(smp, "sample.csv"); out$sample <- smp
  ## the SES ~ genetic PC1 correlation (the popstrat figure), aggregate only
  sg <- data.frame(term = c("r(SES, gPC1)", "r(SES, gPC2)", "r(IQ, gPC1)", "r(SES, IQ)"),
                   r = c(stats::cor(d[[K[["ses"]]]], d$gpc01, use = "pairwise"),
                         stats::cor(d[[K[["ses"]]]], d$gpc02, use = "pairwise"),
                         stats::cor(d[[K[["iq"]]]], d$gpc01, use = "pairwise"),
                         stats::cor(d[[K[["ses"]]]], d[[K[["iq"]]]], use = "pairwise")),
                   n = c(sum(stats::complete.cases(d[, c(K[["ses"]], "gpc01")])),
                         sum(stats::complete.cases(d[, c(K[["ses"]], "gpc02")])),
                         sum(stats::complete.cases(d[, c(K[["iq"]], "gpc01")])),
                         sum(stats::complete.cases(d[, c(K[["ses"]], K[["iq"]])]))))
  ## variance in SES explained by the 10 genetic PCs
  okg <- stats::complete.cases(d[, c(K[["ses"]], gpc_names())])
  sg <- rbind(sg, data.frame(term = "R2(SES ~ gPC1-10)",
                             r = summary(stats::lm(d[[K[["ses"]]]][okg] ~ as.matrix(d[okg, gpc_names()])))$r.squared,
                             n = sum(okg)))
  write(sg, "ses_gpc_correlations.csv"); out$ses_gpc <- sg

  ## ------------------------------------------------- variable catalogue
  vis <- filter_vis(sub_mod(disc, "rs"), md$vi$variable)
  vi <- md$vi[md$vi$variable %in% vis, ]
  cat(sprintf("  %d variables pass the n >= %d filter (of %d candidates)\n", nrow(vi), mr_opts$min_n, nrow(md$vi)))
  if (nrow(vi) < 10 || !all(c(K[["ses"]], K[["iq"]]) %in% vi$variable))
    stop("Too few variables pass the filter, or SES/IQ fail it: check QC counts in sample.csv ",
         "(n_rs per split must exceed mr_opts$min_n).", call. = FALSE)
  print(table(vi$category))
  dropped <- names(K)[!K %in% vi$variable]
  if (length(dropped)) cat("  key variables failing the filter (left out of key-map analyses):",
                           paste(dropped, collapse = ", "), "\n")
  K <- K[K %in% vi$variable]
  write(vi, "variables.csv")

  ## ------------------------------------------------ 2. univariate rankings
  cat("[2] Univariate BWAS for every variable (raw / cov / gpc; discovery and replication)\n")
  maps <- list(); rk <- list()
  for (m in mods) for (adj in c("raw", "cov", "gpc")) for (s in c("disc", "rep")) {
    ds <- sub_mod(if (s == "disc") disc else rep, m)
    M <- matrix(NA_real_, length(feats[[m]]), nrow(vi), dimnames = list(feats[[m]], vi$variable))
    rows <- vector("list", nrow(vi))
    for (j in seq_len(nrow(vi))) {
      a <- assoc_map(ds, feats[[m]], vi$variable[j], m, adj)
      M[, j] <- a$r
      rows[[j]] <- data.frame(modality = m, adjustment = adj, split = s, variable = vi$variable[j],
                              category = vi$category[j], n = a$n, t(map_summary(a$r)))
    }
    maps[[paste(m, adj, s, sep = "|")]] <- M
    rk[[paste(m, adj, s)]] <- do.call(rbind, rows)
    cat(sprintf("    %s %s %s done\n", m, adj, s))
  }
  rk <- do.call(rbind, rk)
  rk <- do.call(rbind, lapply(split(rk, paste(rk$modality, rk$adjustment, rk$split)), function(x) {
    x$rank <- rank(-x$q90_abs_r, ties.method = "min", na.last = "keep"); x }))
  write(rk, "rankings_univariate.csv"); out$rankings <- rk

  ## ranking summary: top-40 composition, discovery/replication agreement,
  ## and how the ranking moves with genetic PC adjustment
  rs_sum <- list()
  for (m in mods) {
    for (adj in c("raw", "cov", "gpc")) for (s in c("disc", "rep")) {
      x <- rk[rk$modality == m & rk$adjustment == adj & rk$split == s, ]
      x <- x[order(x$rank), ]; top <- utils::head(x, 40)
      tab <- table(factor(top$category, levels = sort(unique(vi$category))))
      av <- stats::anova(stats::lm(q90_abs_r ~ category, data = x))
      ses_rank <- x$rank[x$variable == K[["ses"]]]; iq_rank <- x$rank[x$variable == K[["iq"]]]
      rs_sum[[length(rs_sum) + 1]] <- data.frame(
        modality = m, adjustment = adj, split = s, n_variables = nrow(x),
        top_variable = x$variable[1], top_category = x$category[1], top_q90_abs_r = x$q90_abs_r[1],
        ses_rank = if (length(ses_rank)) ses_rank else NA, ses_q90_abs_r = x$q90_abs_r[x$variable == K[["ses"]]][1],
        iq_rank = if (length(iq_rank)) iq_rank else NA, iq_q90_abs_r = x$q90_abs_r[x$variable == K[["iq"]]][1],
        pfac_q90_abs_r = x$q90_abs_r[x$variable == (if ("pfac" %in% names(K)) K[["pfac"]] else NA)][1],
        category_F = av$`F value`[1], category_df1 = av$Df[1], category_df2 = av$Df[2], category_p = av$`Pr(>F)`[1],
        t(stats::setNames(as.numeric(tab), paste0("top40_", names(tab)))))
    }
    for (adj in c("raw", "cov", "gpc")) {
      a <- rk[rk$modality == m & rk$adjustment == adj & rk$split == "disc", ]
      b <- rk[rk$modality == m & rk$adjustment == adj & rk$split == "rep", ]
      b <- b[match(a$variable, b$variable), ]
      rs_sum[[length(rs_sum) + 1]] <- data.frame(modality = m, adjustment = adj, split = "disc_vs_rep",
        n_variables = nrow(a), rho_rank = stats::cor(a$q90_abs_r, b$q90_abs_r, method = "spearman", use = "pairwise"))
    }
    a <- rk[rk$modality == m & rk$adjustment == "cov" & rk$split == "disc", ]
    b <- rk[rk$modality == m & rk$adjustment == "gpc" & rk$split == "disc", ]; b <- b[match(a$variable, b$variable), ]
    rs_sum[[length(rs_sum) + 1]] <- data.frame(modality = m, adjustment = "cov_vs_gpc", split = "disc",
      n_variables = nrow(a), rho_rank = stats::cor(a$q90_abs_r, b$q90_abs_r, method = "spearman", use = "pairwise"),
      mean_pct_change_q90 = 100 * mean((b$q90_abs_r - a$q90_abs_r) / a$q90_abs_r, na.rm = TRUE),
      ses_pct_change_q90 = 100 * (b$q90_abs_r[a$variable == K[["ses"]]] - a$q90_abs_r[a$variable == K[["ses"]]]) /
        a$q90_abs_r[a$variable == K[["ses"]]])
  }
  rs_sum <- Reduce(function(x, y) merge(x, y, all = TRUE), rs_sum)
  write(rs_sum, "rankings_summary.csv"); out$rankings_summary <- rs_sum

  ## ---------------------------------------------- 3. multivariate (ridge)
  if (isTRUE(mr_opts$do_ridge_all)) {
    cat("[3] Ridge, every variable (train discovery, test replication)\n")
    rr <- list()
    for (m in mods) for (adj in c("cov", "gpc")) {
      dtr <- sub_mod(disc, m); dte <- sub_mod(rep, m)
      for (j in seq_len(nrow(vi))) {
        o <- mr_safe(vi$variable[j], ridge_oos(dtr, dte, feats[[m]], vi$variable[j], m, adj))
        if (is.null(o)) next
        rr[[length(rr) + 1]] <- data.frame(modality = m, adjustment = adj, variable = vi$variable[j],
                                           category = vi$category[j], r_oos = o$r, n_train = o$n_train, n_test = o$n_test)
      }
      cat(sprintf("    %s %s done\n", m, adj))
    }
    rr <- do.call(rbind, rr)
    rr <- do.call(rbind, lapply(split(rr, paste(rr$modality, rr$adjustment)), function(x) {
      x$rank <- rank(-x$r_oos, ties.method = "min", na.last = "keep"); x }))
    write(rr, "ridge_all.csv"); out$ridge_all <- rr
  }

  ## ------------------------------------------- 4. maps, exposome, networks
  cat("[4] Brain maps: similarity, exposome PC1, networks, genetic PC maps\n")
  ms <- list(); ex <- list(); np <- list(); gp <- list(); key_maps <- list()
  for (m in mods) {
    dm <- sub_mod(disc, m)
    ## genetic PC maps (adjusted for cov, not for other PCs) and the variance
    ## each brain feature shares with the 10 genetic PCs beyond cov
    gmap <- lapply(gpc_names()[1:3], function(g) assoc_map(dm, feats[[m]], g, m, "cov")$r)
    names(gmap) <- gpc_names()[1:3]
    okg <- cov_ok(dm, cov_terms("gpc", m)) & stats::complete.cases(dm[, feats[[m]]])
    X <- as.matrix(dm[okg, feats[[m]]])
    Z0 <- make_Z(dm[okg, ], cov_terms("cov", m)); Z1 <- make_Z(dm[okg, ], cov_terms("gpc", m))
    r0 <- resid_on(X, Z0); r1 <- resid_on(X, Z1)
    tss <- colSums(sweep(X, 2, colMeans(X))^2)
    dR2 <- (colSums(r0^2) - colSums(r1^2)) / tss
    gp[[m]] <- data.frame(modality = m, feature = feats[[m]], dR2_gpc = dR2, n = sum(okg),
                          r_gpc01 = gmap$gpc01, r_gpc02 = gmap$gpc02, r_gpc03 = gmap$gpc03)
    for (adj in c("raw", "cov", "gpc")) {
      M <- maps[[paste(m, adj, "disc", sep = "|")]]
      M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
      ## exposome map: PC1 of all variable maps (features as observations)
      pc <- stats::prcomp(M, center = TRUE, scale. = FALSE)
      expo <- pc$x[, 1]
      nonses <- vi$variable[vi$category != "socioeconomic"]
      M2 <- M[, intersect(colnames(M), nonses), drop = FALSE]
      expo_ns <- stats::prcomp(M2, center = TRUE, scale. = FALSE)$x[, 1]
      ses_map <- M[, K[["ses"]]]
      cs <- suppressWarnings(stats::cor(M, ses_map))[, 1]
      ce <- suppressWarnings(stats::cor(M, expo))[, 1]
      cg <- suppressWarnings(stats::cor(M, gmap$gpc01))[, 1]
      ex[[length(ex) + 1]] <- data.frame(
        modality = m, adjustment = adj, n_maps = ncol(M),
        pc1_var_explained = pc$sdev[1]^2 / sum(pc$sdev^2),
        r_pc1_ses = stats::cor(expo, ses_map), r_pc1_iq = stats::cor(expo, M[, K[["iq"]]]),
        r_pc1_exclSES_ses = stats::cor(expo_ns, ses_map), r_pc1_pc1exclSES = stats::cor(expo, expo_ns),
        r_pc1_gpc1map = stats::cor(expo, gmap$gpc01), r_ses_gpc1map = stats::cor(ses_map, gmap$gpc01),
        pct_nonSES_maps_absr_gt_.5_with_ses = 100 * mean(abs(cs[names(cs) %in% nonses]) > 0.5),
        pct_maps_absr_gt_.5_with_gpc1map = 100 * mean(abs(cg) > 0.5))
      ms[[length(ms) + 1]] <- data.frame(modality = m, adjustment = adj, variable = colnames(M),
        category = vi$category[match(colnames(M), vi$variable)],
        r_with_ses_map = cs, r_with_exposome_map = ce, r_with_gpc1_map = cg)
      for (k in names(K)) if (K[[k]] %in% colnames(M))
        key_maps[[length(key_maps) + 1]] <- data.frame(modality = m, adjustment = adj, key = k,
                                                       feature = rownames(M), r = M[, K[[k]]])
      key_maps[[length(key_maps) + 1]] <- data.frame(modality = m, adjustment = adj, key = "exposome_pc1",
                                                     feature = rownames(M), r = expo)
      if (m == "rs" && !is.null(md$rs_net)) {
        prof <- c(as.list(K[names(K) %in% c("ses", "iq", "pfac", "sleep", "screen", "adi", "height")]),
                  list(exposome_pc1 = NA))
        for (k in names(prof)) {
          r <- if (k == "exposome_pc1") expo else M[, prof[[k]]]
          p <- net_profile(r, md$rs_net)
          np[[length(np) + 1]] <- data.frame(adjustment = adj, map = k, network = names(p$strength),
                                             strength = p$strength, ms_share = p$ms_share)
        }
        if (adj == "cov") {
          p <- net_profile(gmap$gpc01, md$rs_net)
          np[[length(np) + 1]] <- data.frame(adjustment = "cov", map = "gpc01", network = names(p$strength),
                                             strength = p$strength, ms_share = p$ms_share)
        }
      }
    }
  }
  ex <- do.call(rbind, ex); ms <- do.call(rbind, ms); gp <- do.call(rbind, gp)
  write(ex, "exposome_summary.csv"); write(ms, "map_similarity_all.csv"); write(gp, "gpc_feature_maps.csv")
  write(do.call(rbind, key_maps), "key_maps.csv")
  if (length(np)) write(do.call(rbind, np), "network_profile.csv")
  out$exposome <- ex

  ## key map correlations with participant-bootstrap CIs
  cat("    bootstrap CIs for key map correlations\n")
  keys <- K[names(K) %in% c("ses", "iq", "pfac", "sleep", "screen", "adi", "height")]
  pairs <- c(lapply(setdiff(names(keys), "ses"), function(k) c("ses", k)),
             lapply(names(keys), function(k) c(k, "gpc01")))
  bt <- list()
  for (m in mods) {
    dm <- sub_mod(disc, m)
    one <- function(dd, adj) {
      mp <- lapply(keys, function(v) assoc_map(dd, feats[[m]], v, m, adj)$r)
      mp$gpc01 <- assoc_map(dd, feats[[m]], "gpc01", m, "cov")$r
      vapply(pairs, function(p) suppressWarnings(stats::cor(mp[[p[1]]], mp[[p[2]]])), numeric(1))
    }
    for (adj in c("cov", "gpc")) {
      est <- one(dm, adj)
      bs <- replicate(mr_opts$n_boot_map, one(dm[sample.int(nrow(dm), replace = TRUE), ], adj))
      bs <- matrix(bs, nrow = length(pairs))
      ## basic (reverse-percentile) bootstrap interval. Resampling with
      ## replacement adds noise to both maps, which pulls a map-to-map
      ## correlation toward 0; the percentile interval inherits that bias,
      ## the basic interval (2*est - quantiles) corrects it to first order.
      q_lo <- apply(bs, 1, stats::quantile, .025, na.rm = TRUE)
      q_hi <- apply(bs, 1, stats::quantile, .975, na.rm = TRUE)
      bt[[length(bt) + 1]] <- data.frame(modality = m, adjustment = adj,
        pair = vapply(pairs, paste, character(1), collapse = "~"), r = est,
        lo = pmax(-1, 2 * est - q_hi), hi = pmin(1, 2 * est - q_lo),
        boot_mean = rowMeans(bs, na.rm = TRUE), boot_sd = apply(bs, 1, stats::sd, na.rm = TRUE),
        n_boot = mr_opts$n_boot_map, n = nrow(dm))
    }
  }
  bt <- do.call(rbind, bt); write(bt, "key_map_correlations.csv"); out$key_map_cor <- bt

  ## ------------------------------------------- 5. IQ adjusted for SES
  cat("[5] Brain-IQ associations adjusted for SES (and SES for IQ)\n")
  ia <- list(); iu <- list()
  for (m in mods) {
    dm <- sub_mod(disc, m)
    gmap1 <- assoc_map(dm, feats[[m]], "gpc01", m, "cov")$r
    for (base in c("cov", "gpc")) {
      a0 <- assoc_map(dm, feats[[m]], K[["iq"]], m, base)
      a1 <- assoc_map(dm, feats[[m]], K[["iq"]], m, base, extra = K[["ses"]])
      s0 <- assoc_map(dm, feats[[m]], K[["ses"]], m, base)
      s1 <- assoc_map(dm, feats[[m]], K[["ses"]], m, base, extra = K[["iq"]])
      ## same rows for a fair before/after comparison
      iu[[length(iu) + 1]] <- data.frame(
        modality = m, base = base, n_iq = a0$n, n_iq_adj = a1$n,
        pct_features_abs_r_decrease = 100 * mean(abs(a1$r) < abs(a0$r), na.rm = TRUE),
        slope_change_vs_unadjusted = unname(stats::coef(stats::lm(I(a1$r - a0$r) ~ a0$r))[2]),
        r_change_vs_unadjusted = stats::cor(a1$r - a0$r, a0$r, use = "pairwise"),
        iq_q90_unadj = abs_q(a0$r, .9), iq_q90_sesadj = abs_q(a1$r, .9),
        ses_q90_unadj = abs_q(s0$r, .9), ses_q90_iqadj = abs_q(s1$r, .9),
        r_iqmap_sesmap_unadj = stats::cor(a0$r, s0$r), r_iqmap_sesmap_after = stats::cor(a1$r, s0$r),
        r_iqmap_gpc1map_unadj = stats::cor(a0$r, gmap1), r_iqmap_gpc1map_sesadj = stats::cor(a1$r, gmap1),
        ms_share_iq_unadj = if (m == "rs") net_profile(a0$r, md$rs_net)$ms_share else NA,
        ms_share_iq_sesadj = if (m == "rs") net_profile(a1$r, md$rs_net)$ms_share else NA)
      ## ridge, discovery -> replication, plus split-half resamples of the full sample
      dtr <- sub_mod(disc, m); dte <- sub_mod(rep, m); dall <- sub_mod(d, m)
      for (spec in list(c("iq", ""), c("iq", "ses"), c("ses", ""), c("ses", "iq"))) {
        tgt <- K[[spec[1]]]; ex_ <- if (nzchar(spec[2])) K[[spec[2]]] else NULL
        o <- ridge_oos(dtr, dte, feats[[m]], tgt, m, base, ex_)
        sh <- replicate(mr_opts$n_split_half, {
          fam <- unique(dall$fam); tr <- dall$fam %in% sample(fam, length(fam) %/% 2)
          ridge_oos(dall[tr, ], dall[!tr, ], feats[[m]], tgt, m, base, ex_)$r })
        ia[[length(ia) + 1]] <- data.frame(modality = m, base = base, target = spec[1],
          adjusted_for = if (nzchar(spec[2])) spec[2] else "none", r_oos_disc_rep = o$r,
          r_oos_splithalf_mean = mean(sh, na.rm = TRUE), r_oos_splithalf_sd = stats::sd(sh, na.rm = TRUE),
          n_train = o$n_train, n_test = o$n_test)
      }
    }
  }
  iu <- do.call(rbind, iu); ia <- do.call(rbind, ia)
  write(iu, "iq_ses_univariate.csv"); write(ia, "iq_ses_ridge.csv")
  out$iq_ses_uni <- iu; out$iq_ses_ridge <- ia

  ## ------------------------------- 6. SES-stratified generalization of IQ
  cat("[6] Generalization of brain-IQ models trained on SES strata\n")
  st <- list(); gr <- list()
  for (m in mods) {
    dtr0 <- sub_mod(disc, m); dte0 <- sub_mod(rep, m)
    need <- c(K[["iq"]], K[["ses"]])
    dtr0 <- dtr0[stats::complete.cases(dtr0[, need]), ]; dte0 <- dte0[stats::complete.cases(dte0[, need]), ]
    ses <- dtr0[[K[["ses"]]]]; iq <- dtr0[[K[["iq"]]]]
    hi <- which(ses >= stats::quantile(ses, mr_opts$high_q)); lo <- which(ses <= stats::quantile(ses, mr_opts$low_q))
    n_sub <- min(length(hi), length(lo))
    brks <- unique(stats::quantile(dte0[[K[["iq"]]]], seq(0, 1, length.out = mr_opts$iq_bins + 1)))
    pbin <- as.numeric(table(cut(dte0[[K[["iq"]]]], brks, include.lowest = TRUE))) / nrow(dte0)
    bin_tr <- cut(iq, brks, include.lowest = TRUE)
    draw_matched <- function(pool) {
      tgt <- round(n_sub * pbin); idx <- integer(0)
      for (b in seq_along(tgt)) {
        cand <- pool[as.integer(bin_tr[pool]) %in% b]
        if (length(cand)) idx <- c(idx, cand[sample.int(length(cand), min(tgt[b], length(cand)))])
      }
      idx
    }
    groups <- list(full = seq_len(nrow(dtr0)), high_ses = hi, low_ses = lo)
    for (adj in c("cov", "gpc")) for (g in names(groups)) for (matched in c(TRUE, FALSE)) {
      res <- replicate(mr_opts$n_boot_strat, {
        idx <- if (matched) draw_matched(groups[[g]]) else groups[[g]][sample.int(length(groups[[g]]), min(n_sub, length(groups[[g]])))]
        o <- ridge_oos(dtr0[idx, ], dte0, feats[[m]], K[["iq"]], m, adj)
        c(r_oos_iq = o$r,
          r_pred_ses = if (!is.null(o$pred)) suppressWarnings(stats::cor(o$pred, o$test[[K[["ses"]]]])) else NA,
          r_vars_train = stats::cor(dtr0[[K[["ses"]]]][idx], dtr0[[K[["iq"]]]][idx]),
          n_train = length(idx))
      })
      st[[length(st) + 1]] <- data.frame(modality = m, adjustment = adj, train_group = g, iq_matched = matched,
        n_train = round(mean(res["n_train", ])), n_test = nrow(dte0),
        r_oos_iq_mean = mean(res["r_oos_iq", ], na.rm = TRUE), r_oos_iq_sd = stats::sd(res["r_oos_iq", ], na.rm = TRUE),
        r_oos_iq_lo = stats::quantile(res["r_oos_iq", ], .025, na.rm = TRUE),
        r_oos_iq_hi = stats::quantile(res["r_oos_iq", ], .975, na.rm = TRUE),
        r_pred_ses_mean = mean(res["r_pred_ses", ], na.rm = TRUE),
        r_vars_train_mean = mean(res["r_vars_train", ], na.rm = TRUE), n_boot = mr_opts$n_boot_strat)
    }
    ## composition gradient: from high-SES only to the full SES range
    for (adj in c("cov", "gpc")) for (f in seq(0, 1, by = 0.1)) {
      res <- replicate(mr_opts$n_boot_grad, {
        nf <- round(f * n_sub)
        idx <- c(sample(seq_len(nrow(dtr0)), nf), sample(hi, n_sub - nf))
        o <- ridge_oos(dtr0[idx, ], dte0, feats[[m]], K[["iq"]], m, adj)
        o2 <- ridge_oos(dtr0[idx, ], dte0, feats[[m]], K[["ses"]], m, adj)
        c(r_oos_iq = o$r,
          r_pred_iq_on_ses = if (!is.null(o$pred)) suppressWarnings(stats::cor(o$pred, o$test[[K[["ses"]]]])) else NA,
          r_oos_ses = o2$r,
          r_vars_train = stats::cor(dtr0[[K[["ses"]]]][idx], dtr0[[K[["iq"]]]][idx]))
      })
      gr[[length(gr) + 1]] <- data.frame(modality = m, adjustment = adj, frac_full_range = f, n_train = n_sub,
        n_test = nrow(dte0), t(rowMeans(res, na.rm = TRUE)),
        r_oos_iq_sd = stats::sd(res["r_oos_iq", ], na.rm = TRUE), n_boot = mr_opts$n_boot_grad)
    }
  }
  st <- do.call(rbind, st); gr <- do.call(rbind, gr)
  write(st, "ses_stratified_generalization.csv"); write(gr, "ses_composition_gradient.csv")
  out$strat <- st; out$gradient <- gr

  ## ------------------------------------------------------ 7. Year 2
  if (isTRUE(mr_opts$do_year2) && !is.null(md$y2) && "rs" %in% mods) {
    cat("[7] Year-2 replication of the RSFC ranking\n")
    y2 <- sub_mod(md$y2, "rs")
    v2 <- filter_vis(y2, vi$variable)
    r2 <- do.call(rbind, lapply(v2, function(v) {
      a <- assoc_map(y2, feats$rs, v, "rs", "raw")
      data.frame(variable = v, n = a$n, q90_abs_r_y2 = abs_q(a$r, .9))
    }))
    b <- rk[rk$modality == "rs" & rk$adjustment == "raw" & rk$split == "disc", ]
    r2$q90_abs_r_baseline <- b$q90_abs_r[match(r2$variable, b$variable)]
    r2$category <- vi$category[match(r2$variable, vi$variable)]
    write(r2, "year2_rankings.csv")
    write(data.frame(n_variables = nrow(r2),
                     rho = stats::cor(r2$q90_abs_r_y2, r2$q90_abs_r_baseline, method = "spearman", use = "pairwise"),
                     n_participants = nrow(y2)), "year2_summary.csv")
  }

  ## ------------------------------------------------------ 8. figures
  mr_safe("figures", {
    grDevices::pdf(file.path(out_dir, "figures.pdf"), width = 10, height = 6)
    cats <- sort(unique(vi$category)); col <- grDevices::hcl.colors(length(cats), "Dark 3")
    for (m in mods) for (adj in c("raw", "gpc")) {
      x <- rk[rk$modality == m & rk$adjustment == adj & rk$split == "disc", ]; x <- x[order(x$rank), ]
      graphics::plot(seq_len(nrow(x)), x$q90_abs_r, col = col[match(x$category, cats)], pch = 19, cex = .6,
                     xlab = "variables ranked", ylab = "90th percentile |r| across features",
                     main = sprintf("%s univariate BWAS ranking (%s)", toupper(m), adj))
      graphics::legend("topright", cats, col = col, pch = 19, cex = .7, bty = "n")
    }
    for (m in mods) {
      g <- gr[gr$modality == m, ]
      graphics::plot(NA, xlim = c(0, 1), ylim = range(c(0, g$r_oos_iq, g$r_oos_ses), na.rm = TRUE),
                     xlab = "fraction of training sample drawn from the full SES range (0 = high SES only)",
                     ylab = "out-of-sample r", main = sprintf("%s: IQ and SES models by training SES composition", toupper(m)))
      for (adj in c("cov", "gpc")) {
        h <- g[g$adjustment == adj, ]; lt <- if (adj == "cov") 1 else 2
        graphics::lines(h$frac_full_range, h$r_oos_iq, col = "darkgreen", lty = lt, lwd = 2)
        graphics::lines(h$frac_full_range, h$r_oos_ses, col = "navy", lty = lt, lwd = 2)
        graphics::lines(h$frac_full_range, h$r_vars_train, col = "grey50", lty = lt)
      }
      graphics::legend("topleft", c("IQ", "SES", "r(SES, IQ) in training", "solid = cov, dashed = + genetic PCs"),
                       col = c("darkgreen", "navy", "grey50", NA), lty = c(1, 1, 1, NA), bty = "n", cex = .8)
    }
    grDevices::dev.off()
  })
  cat("\nDone. Aggregate outputs in ", out_dir, "\n")
  invisible(out)
}
