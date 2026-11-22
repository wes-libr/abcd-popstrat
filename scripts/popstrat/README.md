# Marek et al. (2026) replication on ABCD 7.0, with genetic-PC adjustment

Replicates the main analyses of Marek et al., "Patterns of brain-wide associations
reflect socioeconomics" (*Science* 392:eaee6213, 2026). It also asks the popstrat
paper's question: how much of the socioeconomic brain pattern survives adjustment
for population structure (genetic PCs 1–10)?

## Run

```r
setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
source("scripts/popstrat/marek_replication.R")
```

The first run pulls about 1,600 variables with NBDCtools and caches the analytic
object in `dir_derived/marek_rep_data.rds`. That is participant-level, host-only
and outside the repo. Later runs reuse it. Set `mr_rebuild <- TRUE` to pull again.

For a quick first pass (minutes, not hours):

```r
mr_opts_user <- list(do_ridge_all = FALSE, n_boot_map = 20, n_split_half = 5,
                     n_boot_strat = 20, n_boot_grad = 10)
source("scripts/popstrat/marek_replication.R")
```

All outputs go to `output/popstrat_marek/`. They are aggregate only and go through
`safe_write_csv()`, so counts under 10 are suppressed.

## What is replicated, and how it differs from Marek

| Marek | Here (ABCD 7.0 tabulated) |
|---|---|
| ABCD 5.1, DCAN pipeline, 333 Gordon parcels, 55,278 edges | Gordon **network-pair** RSFC, 13 networks, 91 pairs (`mr_y_rsfmri__corr__gpnet`); parcel edges are not in the tabulated release |
| Cortical thickness, 59,412 vertices | Destrieux thickness, 148 regions |
| ≥ 5 min at FD < 0.08 mm | ABCD rsfMRI inclusion flag and ≥ 375 frames retained at FD < 0.2 mm |
| 649 curated variables | ~400 automated summary scores from the same table families, filtered by Marek's n ≥ 2000 rule |
| Matched ARMS discovery/replication split | Random half of families within site |
| Spin tests on parcel maps | Participant-bootstrap CIs (basic bootstrap) on map correlations; spin tests need surface parcels |
| CCA/ridge on brain PCs (20% variance) | Ridge on all features, λ by 5-fold CV |
| SES = COI social and economic domain | Same: `le_l_coi__addr1__se__total__national_zscore` |
| IQ = NIH Toolbox total | `nc_y_nihtb__comp__tot__uncor_score` (age is a covariate) |
| p-factor = CBCL total | `mh_p_cbcl_sum` |

Adjustment sets: **raw** (bivariate r, as in Marek), **cov** (age, sex, site,
scanner manufacturer, plus mean FD for RSFC), and **gpc** (cov plus genetic PCs
1–10, `gn_y_popstruct_pc__01..10`).

## Outputs, by Marek's claims

| Claim | File |
|---|---|
| SES has the strongest brain-wide associations (Fig. 1) | `rankings_univariate.csv`, `rankings_summary.csv` (top-40 composition, category ANOVA, SES rank, discovery vs replication ρ, change with genetic PCs) |
| Multivariate associations agree (Fig. 1D) | `ridge_all.csv` |
| Year-2 replication (fig. S2) | `year2_rankings.csv`, `year2_summary.csv` |
| One shared "exposome" map ≈ SES map (Fig. 2D–E) | `exposome_summary.csv`, `map_similarity_all.csv`, `key_maps.csv` |
| Motor/sensory networks dominate (Fig. 2F, 3B) | `network_profile.csv` (`ms_share` = share of association strength in smh, smm, vis, aud) |
| IQ map overlaps SES map; adjusting for SES weakens it (Fig. 4) | `iq_ses_univariate.csv`, `iq_ses_ridge.csv` |
| IQ models fail to generalize when trained on high-SES children (Fig. 5) | `ses_stratified_generalization.csv` |
| Shortcut learning as SES composition varies (Fig. 6) | `ses_composition_gradient.csv`, `figures.pdf` |
| **Popstrat:** SES is correlated with genetic PCs | `ses_gpc_correlations.csv` |
| **Popstrat:** variance in each brain feature explained by genetic PCs; genetic PC1 map | `gpc_feature_maps.csv` |
| **Popstrat:** each map's similarity to the genetic PC1 map | `map_similarity_all.csv` (`r_with_gpc1_map`), `key_map_correlations.csv` |

Each analysis is run with and without genetic PCs, so every Marek result has a
"gpc" counterpart.

## Reading the genetic-PC results

- Map similarity alone can't separate confounding from a shared true effect. If
  SES and genetic PC1 both affect the same networks, the SES map stays similar to
  the PC1 map even after adjustment. The decisive quantities are the change in
  association strength (`ses_pct_change_q90` in `rankings_summary.csv`) and the
  ridge r_oos with and without genetic PCs.
- Adjusting for genetic PCs can also remove real environmental signal that is
  structured by ancestry. That is the popstrat paper's main caveat: genetic PCs
  index population structure, not biology. Report the gpc results as a
  sensitivity bound, not as the "true" effect.

## Files

- `marek_replication.R`: runner (sources `00_utils.R`, builds, runs, writes)
- `mr_data.R`: `mr_build()`, the NBDCtools pull and caching
- `mr_analysis.R`: `mr_run()` and the analysis functions (no data I/O)

The analysis code was tested end to end on synthetic data in Claude's sandbox,
including the variable-selection logic against the real 7.0 dictionary names. It
has not yet been run on real data.
