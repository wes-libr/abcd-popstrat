########################################################################
##  popstrat_main.R — the PopStrat paper's analyses on ABCD 7.0, with
##  extensions (family SES, siblings, genetically similar subsamples,
##  head motion). See README_popstrat.md.
##
##  Run from analyses/R:
##      setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
##      source("scripts/popstrat/popstrat_main.R")
##  Quick pass:  ps_opts_user <- list(n_div_boot = 5, n_match_boot = 5)   (defaults: 100 each)
##
##  Cache : dir_derived/marek_rep_data.rds and popstrat_data.rds (host-only)
##  Output: output/popstrat_main/*.csv (aggregate; counts < 10 suppressed)
########################################################################

cat("\n════════════════════════════════════════════\n")
cat("  PopStrat analyses, ABCD 7.0\n")
cat("════════════════════════════════════════════\n")

if (!exists("paths") || !exists("safe_write_csv")) source("00_utils.R")
.ps_dir <- file.path("scripts", "popstrat")
for (f in c("mr_analysis.R", "mr_data.R", "ps_data.R", "ps_analysis.R")) source(file.path(.ps_dir, f))

md <- mr_build()                                              # reuses the Marek cache
pd <- ps_build(md, rebuild = isTRUE(get0("ps_rebuild", ifnotfound = FALSE)))
rm(md)
out_dir <- file.path(paths$output, "popstrat_main")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ps_write <- function(df, name) {
  cc <- setdiff(grep("^n($|_)", names(df), value = TRUE), c("n_boot", "n_nonwhite"))   # resamples / design counts, not cells
  safe_write_csv(df, file.path(out_dir, name), count_cols = cc,
                 key_cols = c("modality", "variable", "adjustment", "diversity", "group", "sample",
                              "term", "analysis"))
}

t0 <- Sys.time()
ps_run(pd, out_dir, write = ps_write)
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rm(pd)
