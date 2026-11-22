## tableS1_idps.R -- list of the ABCD Release 7.0 IDPs used in the PopStrat
## analyses (Table S1). Uses only the public data dictionary, no participant data.
## Same selection rules as ps_data.R (sMRI, DTI, n-back) and mr_data.R (RSFC).
dd <- as.data.frame(NBDCtools::get_dd("abcd", release = "7.0"))
nm <- dd$name
pick <- function(re) grep(re, nm, value = TRUE)
rs <- pick("^mr_y_rsfmri__corr__gpnet__[a-z]{3}__[a-z]{3}_mean$")
pr <- do.call(rbind, strsplit(sub("^mr_y_rsfmri__corr__gpnet__([a-z]{3})__([a-z]{3})_mean$", "\\1 \\2", rs), " "))
rs <- rs[pr[, 1] <= pr[, 2]]                                   # 169 ordered -> 91 unique pairs
L <- list(
  "Structural|Cortical thickness"     = pick("^mr_y_smri__thk__dsk__[a-z0-9]+__[lr]h_mean$"),
  "Structural|Cortical surface area"  = pick("^mr_y_smri__area__dsk__[a-z0-9]+__[lr]h_sum$"),
  "Structural|Subcortical volume"     = pick("^mr_y_smri__vol__aseg__(ab|ag|cd|hc|pl|pt|tp|vdc)__[lr]h_sum$"),
  "Diffusion|Fractional anisotropy"   = pick("^mr_y_dti__fs__fa__at__[a-z0-9]+__[lr]h_wmean$"),
  "Resting-state fMRI|Network connectivity" = rs,
  "Task fMRI (n-back)|2-back vs 0-back" = pick("^mr_y_tfmri__nback__2bv0b__dsk__[a-z0-9]+__[lr]h_beta$"))
out <- do.call(rbind, lapply(names(L), function(k) {
  v <- L[[k]]; p <- strsplit(k, "\\|")[[1]]
  data.frame(modality = p[1], measure = p[2], variable = v,
             label = dd$label[match(v, nm)], atlas = dd$atlas[match(v, nm)])
}))
stopifnot(nrow(out) == 405)
write.csv(out, "tableS1_idps_abcd7.csv", row.names = FALSE)
