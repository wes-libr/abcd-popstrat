########################################################################
##  px_supp.R — data for the ABCD 7.0 versions of supplementary Figs S1
##  and S3 of the PopStrat manuscript. Aggregate output only.
##
##  S1  Distribution of genetic PCs 1-3 by self-identified race/ethnicity.
##      The 5.0 figure plotted one point per child; here we export
##        px_S1_gpc_quantiles.csv : per group and PC, n and the 2.5/25/50/
##                                  75/97.5th percentiles
##        px_S1_gpc_bins.csv      : per group and PC pair (1-2, 1-3, 2-3),
##                                  counts on a 25 x 25 grid over the
##                                  pooled 0.5-99.5th percentile range;
##                                  cells with < 10 children are suppressed
##  S3  Variance in each VI explained by genetic PCs 1-10 beyond age, sex
##      and site, by diversity proportion (n = 2,000; px_opts$n_div_boot
##      resamples), with 95% intervals and the null expectation 10 / n:
##        px_S3_vi_gpc_dR2_by_diversity.csv
##      and in the full baseline sample (px_S3_vi_gpc_dR2_full.csv).
########################################################################

px_supp_run <- function(pd, out_dir, write = function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE) }) {
  set.seed(px_opts$seed + 1); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  d <- pd$base
  d$nonwhite <- ifelse(is.na(d$race), NA, d$race != "White")
  V <- ps_vis[intersect(px_opts$vis, names(ps_vis))]; V <- V[V %in% names(d)]
  base <- c("age", "sex", "site")

  ## ------------------------------------------------------------- S1
  cat("\n[PX-S1] Genetic PC distributions by group (aggregate)\n")
  g <- d[!is.na(d$race) & stats::complete.cases(d[, c("gpc01", "gpc02", "gpc03")]), ]
  q <- do.call(rbind, lapply(sort(unique(g$race)), function(r) do.call(rbind, lapply(c("gpc01", "gpc02", "gpc03"), function(p) {
    x <- g[[p]][g$race == r]
    data.frame(group = r, pc = p, n = length(x),
               t(stats::setNames(stats::quantile(x, c(.025, .25, .5, .75, .975)), c("q025", "q25", "q50", "q75", "q975"))))
  }))))
  write(q, "px_S1_gpc_quantiles.csv")
  nb <- 25
  rng <- lapply(c(gpc01 = "gpc01", gpc02 = "gpc02", gpc03 = "gpc03"), function(p) stats::quantile(g[[p]], c(.005, .995)))
  brk <- lapply(rng, function(r) seq(r[1], r[2], length.out = nb + 1))
  bins <- do.call(rbind, lapply(list(c("gpc01", "gpc02"), c("gpc01", "gpc03"), c("gpc02", "gpc03")), function(pr) {
    do.call(rbind, lapply(sort(unique(g$race)), function(r) {
      x <- g[g$race == r, pr]
      ix <- cut(x[[1]], brk[[pr[1]]], include.lowest = TRUE, labels = FALSE)
      iy <- cut(x[[2]], brk[[pr[2]]], include.lowest = TRUE, labels = FALSE)
      ok <- !is.na(ix) & !is.na(iy)
      tb <- as.data.frame(table(factor(ix[ok], 1:nb), factor(iy[ok], 1:nb)))
      names(tb) <- c("bx", "by", "n"); tb <- tb[tb$n > 0, ]
      bx <- as.integer(as.character(tb$bx)); by <- as.integer(as.character(tb$by))
      data.frame(pair = paste(pr, collapse = "_"), group = r,
                 x_mid = (brk[[pr[1]]][bx] + brk[[pr[1]]][bx + 1]) / 2,
                 y_mid = (brk[[pr[2]]][by] + brk[[pr[2]]][by + 1]) / 2, n = tb$n)
    }))
  }))
  write(bins, "px_S1_gpc_bins.csv")

  ## ------------------------------------------------------------- S3
  cat("[PX-S3] Variance in VIs explained by genetic PCs, by diversity\n")
  r2 <- function(s, y) {
    s <- s[!is.na(s[[y]]), ]; if (nrow(s) < 200) return(NA_real_)
    e0 <- sum(resid_on(s[[y]], make_Z(s, base))^2); e1 <- sum(resid_on(s[[y]], make_Z(s, c(base, gpc_names())))^2)
    (e0 - e1) / e0
  }
  pool <- d[!is.na(d$nonwhite) & stats::complete.cases(d[, c(gpc_names(), base)]), ]
  W <- which(!pool$nonwhite); NW <- which(pool$nonwhite)
  n_use <- min(px_opts$n_div, length(W), length(NW))
  S3 <- list()
  for (p in px_opts$props) {
    n_nw <- round(p * n_use); n_w <- n_use - n_nw
    M <- t(replicate(px_opts$n_div_boot, { s <- pool[c(sample(W, n_w), sample(NW, n_nw)), ]
      vapply(V, function(v) r2(s, v), 0) }))
    for (k in names(V)) S3[[length(S3) + 1]] <- data.frame(diversity = p, variable = k, n = n_use,
      dR2_mean = mean(M[, k], na.rm = TRUE), dR2_lo = stats::quantile(M[, k], .025, na.rm = TRUE),
      dR2_hi = stats::quantile(M[, k], .975, na.rm = TRUE), null_expected = length(gpc_names()) / n_use,
      n_boot = sum(!is.na(M[, k])), row.names = NULL)
    cat(sprintf("    diversity %.2f done\n", p))
  }
  write(do.call(rbind, S3), "px_S3_vi_gpc_dR2_by_diversity.csv")
  full <- do.call(rbind, lapply(names(V), function(k) {
    s <- pool[!is.na(pool[[V[[k]]]]), ]
    data.frame(variable = k, n = nrow(s), dR2 = r2(s, V[[k]]), null_expected = length(gpc_names()) / nrow(s))
  }))
  write(full, "px_S3_vi_gpc_dR2_full.csv")
  cat("\nDone. Aggregate outputs in ", out_dir, "\n")
  invisible(TRUE)
}
