########################################################################
##  ses_causal.R — the causal-inference ladder applied to Marek et al.'s
##  exposure (neighborhood SES) with the brain and IQ as outcomes.
##
##  Run from analyses/R:
##      setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
##      source("scripts/popstrat/ses_causal.R")
##  Quick pass:  sc_opts_user <- list(n_boot = 50)
##
##  Cache : dir_derived/ses_causal_panel.rds (host-only)
##  Output: output/ses_causal/sc_ladder.csv, sc_sample.csv, sc_movers.csv
########################################################################

cat("\n════════════════════════════════════════════\n")
cat("  Causal ladder: neighborhood SES -> brain, IQ\n")
cat("════════════════════════════════════════════\n")

if (!exists("paths") || !exists("safe_write_csv")) source("00_utils.R")
.sc_dir <- file.path("scripts", "popstrat")
for (f in c("mr_analysis.R", "mr_data.R", "ps_analysis.R", "sc_data.R", "sc_analysis.R")) source(file.path(.sc_dir, f))

md <- mr_build()
sp <- sc_build(md, rebuild = isTRUE(get0("sc_rebuild", ifnotfound = FALSE)))
rm(md)
out_dir <- file.path(paths$output, "ses_causal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sc_write <- function(df, name) {
  cc <- intersect(c("n", "n_id", "n_fam", "n_coi", "n_outcome"), names(df))
  safe_write_csv(df, file.path(out_dir, name), count_cols = cc,
                 key_cols = c("method", "outcome", "adjustment", "wave", "term", "note"))
}

t0 <- Sys.time()
sc_run(sp, out_dir, write = sc_write)
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rm(sp)
