########################################################################
##  px_data.R — add what the PopStrat extension analyses (px_analysis.R)
##  need to the PopStrat analytic object; cached in
##  dir_derived/popstrat_ext_data.rds (participant-level, host-only).
##
##  Adds to pd$base (baseline, ses-00A):
##    icv      : intracranial volume (mr_y_smri__vol__aseg__icv_sum)
##    gpc11... : any genetic PCs beyond 10 that the release provides
##               (gn_y_popstruct_pc__11, ...); px_gpc_all() lists them
##    serial   : scanner serial number, if the release provides it
##    scanner  : site x manufacturer (always), site x serial (if serial)
##  Adds pd$feat_type: for every sMRI feature, "thickness", "area" or
##  "volume", used to split the head-size and subsample analyses.
########################################################################

px_build <- function(pd, cache_file = file.path(paths$cache, "popstrat_ext_data.rds"), rebuild = FALSE) {
  if (file.exists(cache_file) && !rebuild) { message("Using cached ", cache_file); return(readRDS(cache_file)) }
  suppressPackageStartupMessages(library(NBDCtools))
  nm <- as.data.frame(NBDCtools::get_dd("abcd", release = "7.0"))$name
  pcs <- grep("^gn_y_popstruct_pc__[0-9]+$", nm, value = TRUE)
  pcs <- pcs[as.integer(sub(".*__", "", pcs)) > 10]
  icv <- intersect("mr_y_smri__vol__aseg__icv_sum", nm)
  ser <- grep("^mr_y_adm__info__dev_serial", nm, value = TRUE)[1]
  want <- c(icv, pcs, if (!is.na(ser)) ser)
  message(sprintf("Pulling %d variables: ICV %s, %d extra genetic PCs, serial %s",
                  length(want), if (length(icv)) "yes" else "NOT FOUND", length(pcs), if (is.na(ser)) "not found" else ser))
  d <- pd$base
  if (length(want)) {
    raw <- as.data.frame(NBDCtools::create_dataset_abcd(dir_data = paths$abcd_data, vars = want, value_to_na = TRUE))
    raw$id <- as.character(raw$participant_id)
    num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    ## genetic PCs are static (session NA); take the first non-missing value per child
    first <- function(v) tapply(num(raw[[v]]), raw$id, function(x) x[!is.na(x)][1])
    for (v in pcs) d[[sprintf("gpc%02d", as.integer(sub(".*__", "", v)))]] <- as.numeric(first(v)[d$id])
    b <- raw[raw$session_id %in% "ses-00A", , drop = FALSE]; i <- match(d$id, b$id)
    if (length(icv)) d$icv <- num(b[[icv]][i])
    if (!is.na(ser)) d$serial <- as.character(b[[ser]][i])
  }
  if (is.null(d$icv)) d$icv <- NA_real_
  d$scanner <- droplevels(interaction(d$site, d$manuf, drop = TRUE))
  d$scanner_serial <- if (!is.null(d$serial)) droplevels(interaction(d$site, d$serial, drop = TRUE)) else d$scanner
  pd$base <- d
  pd$feat_type <- ifelse(grepl("__thk__", pd$feat$smri), "thickness",
                         ifelse(grepl("__area__", pd$feat$smri), "area", "volume"))
  names(pd$feat_type) <- pd$feat$smri
  saveRDS(pd, cache_file)
  message("Cached PopStrat extension object (participant-level, host-only): ", cache_file)
  pd
}

## every genetic PC available, in order
px_gpc_all <- function(d) {
  g <- grep("^gpc[0-9]+$", names(d), value = TRUE)
  g <- g[vapply(g, function(v) any(!is.na(d[[v]])), TRUE)]
  g[order(as.integer(sub("gpc", "", g)))]
}
