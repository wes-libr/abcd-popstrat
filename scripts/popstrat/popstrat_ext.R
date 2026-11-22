########################################################################
##  popstrat_ext.R — extension analyses for the PopStrat manuscript
##  (ABCD 7.0): corrected map similarity with permutation nulls, head and
##  body size, robustness of the population-structure adjustment, and the
##  7.0 versions of draft Figs 3-5. See README_popstrat.md, section
##  "Extensions (November 2026 review)".
##
##  Run from analyses/R, after popstrat_main.R (the figures read its ps_A
##  tables):
##      setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
##      source("scripts/popstrat/popstrat_ext.R")
##  Quick pass:  px_opts_user <- list(n_div_boot = 5, n_perm_div = 20, n_perm_full = 100,
##                                    n_cboot = 20, n_match_boot = 10)
##
##  Cache : dir_derived/popstrat_ext_data.rds (host-only)
##  Output: output/popstrat_ext/*.csv and popstrat_7p0_figures.pdf
##          (aggregate; counts < 10 suppressed)
########################################################################

cat("\n════════════════════════════════════════════\n")
cat("  PopStrat extension analyses, ABCD 7.0\n")
cat("════════════════════════════════════════════\n")

if (!exists("paths") || !exists("safe_write_csv")) source("00_utils.R")
.ps_dir <- file.path("scripts", "popstrat")
for (f in c("mr_analysis.R", "mr_data.R", "ps_data.R", "ps_analysis.R", "px_data.R", "px_analysis.R"))
  source(file.path(.ps_dir, f))

md <- mr_build()
pd <- ps_build(md); rm(md)
pd <- px_build(pd, rebuild = isTRUE(get0("px_rebuild", ifnotfound = FALSE)))
out_dir <- file.path(paths$output, "popstrat_ext")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

px_write <- function(df, name) {
  cc <- intersect(c("n", "n_families"), names(df))      # people; n_boot, n_perm, n_features are not
  safe_write_csv(df, file.path(out_dir, name), count_cols = cc,
                 key_cols = c("modality", "variable", "adjustment", "diversity", "group", "measure_type"))
}

t0 <- Sys.time()
px_run(pd, out_dir, write = px_write)
rm(pd)
px_figures(file.path(paths$output, "popstrat_main"), out_dir)
cat(sprintf("Figures: %s\nElapsed: %.1f min\n", file.path(out_dir, "popstrat_7p0_figures.pdf"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
