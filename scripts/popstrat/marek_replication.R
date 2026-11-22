########################################################################
##  marek_replication.R — replicate Marek et al. (2026, Science
##  392:eaee6213, "Patterns of brain-wide associations reflect
##  socioeconomics") on ABCD 7.0, and test the popstrat question:
##  how much of the SES brain pattern is population structure?
##
##  Run from analyses/R:
##      setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
##      source("scripts/popstrat/marek_replication.R")
##
##  Optional settings, before sourcing:
##      mr_opts_user <- list(n_boot_strat = 20, do_ridge_all = FALSE)  # quick pass
##      mr_rebuild   <- TRUE                                           # re-pull data
##
##  Inputs : ABCD 7.0 tabulated release via NBDCtools (paths$abcd_data)
##  Cache  : dir_derived/marek_rep_data.rds (participant-level, HOST-ONLY)
##  Outputs: output/popstrat_marek/*.csv + figures.pdf (aggregate only;
##           counts < 10 suppressed by safe_write_csv)
########################################################################

cat("\n════════════════════════════════════════════\n")
cat("  Marek et al. (2026) replication, ABCD 7.0\n")
cat("════════════════════════════════════════════\n")

if (!exists("paths") || !exists("safe_write_csv")) source("00_utils.R")
.mr_dir <- file.path("scripts", "popstrat")
source(file.path(.mr_dir, "mr_analysis.R"))
source(file.path(.mr_dir, "mr_data.R"))

md <- mr_build(rebuild = isTRUE(get0("mr_rebuild", ifnotfound = FALSE)))
out_dir <- file.path(paths$output, "popstrat_marek")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## every table goes through the project's small-cell suppression
mr_write <- function(df, name) {
  cc <- setdiff(grep("^n($|_)|^n_", names(df), value = TRUE), "n_boot")   # n_boot counts resamples, not people
  safe_write_csv(df, file.path(out_dir, name), count_cols = cc,
                 key_cols = c("modality", "adjustment", "split", "variable", "category", "feature",
                              "key", "network", "map", "pair", "term", "base", "target",
                              "adjusted_for", "train_group", "iq_matched", "frac_full_range"))
}

t0 <- Sys.time()
res <- mr_run(md, out_dir, write = mr_write)
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rm(md)   # participant-level object; nothing individual-level is left in the workspace
