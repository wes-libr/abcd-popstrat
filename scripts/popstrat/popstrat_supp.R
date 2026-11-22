########################################################################
##  popstrat_supp.R — data for the ABCD 7.0 versions of supplementary
##  Figs S1 and S3 (see px_supp.R). Uses the cached extension object.
##
##  Run from analyses/R (about 1-2 minutes):
##      setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
##      source("scripts/popstrat/popstrat_supp.R")
##  Output: output/popstrat_ext/px_S1_*.csv, px_S3_*.csv
##          (aggregate; counts < 10 suppressed)
########################################################################

if (!exists("paths") || !exists("safe_write_csv")) source("00_utils.R")
.ps_dir <- file.path("scripts", "popstrat")
for (f in c("mr_analysis.R", "mr_data.R", "ps_data.R", "ps_analysis.R", "px_data.R", "px_analysis.R", "px_supp.R"))
  source(file.path(.ps_dir, f))
md <- mr_build(); pd <- ps_build(md); rm(md); pd <- px_build(pd)
out_dir <- file.path(paths$output, "popstrat_ext")
px_write <- function(df, name) {
  safe_write_csv(df, file.path(out_dir, name), count_cols = intersect("n", names(df)),
                 key_cols = c("group", "pc", "pair", "variable", "diversity", "x_mid", "y_mid"))
}
t0 <- Sys.time()
px_supp_run(pd, out_dir, write = px_write)
rm(pd)
cat(sprintf("Elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
