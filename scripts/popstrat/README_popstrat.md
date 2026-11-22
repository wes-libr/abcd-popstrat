# PopStrat analyses on ABCD 7.0

The analysis behind "The confounding effects of population structure in large-scale
neurobehavioral studies" (draft of May 2024, ABCD 5.0), rebuilt for 7.0. It is
extended with the checks raised by the Marek et al. (2026) comparison.

## Run

```r
setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
ps_opts_user <- list(n_div_boot = 5, n_match_boot = 5)   # quick pass; omit for the full run
source("scripts/popstrat/popstrat_main.R")
```

It reuses the Marek cache (`marek_rep_data.rds`) and adds a second host-only cache
(`popstrat_data.rds`) with the extra imaging measures, race/ethnicity, twin and
sibling structure, ADMIXTURE proportions and motion. Set `ps_rebuild <- TRUE` to
pull again. Outputs go to `output/popstrat_main/` (aggregate only, counts < 10
suppressed).

## Imaging-derived phenotypes (baseline)

| Set | 7.0 source | Features | QC |
|---|---|---|---|
| smri | Desikan thickness + area, aseg subcortical volume | ~150 | T1 inclusion |
| dti | AtlasTrack FA | ~40 tracts | dMRI inclusion; dMRI mean FD covariate |
| rs | Gordon network-pair correlations | 91 | rsfMRI inclusion, ≥ 375 frames; mean FD covariate |
| nback | 2-back vs 0-back betas, Desikan | 68 | n-back inclusion; n-back mean FD covariate |

Variables of interest: COI social and economic domain (`ses`), ADI national
percentile, NIH Toolbox fluid and total composites, CBCL internalizing,
height and weight. Covariates (`cov`): age, sex, site, scanner manufacturer
and modality motion. Genetic PCs are `gn_y_popstruct_pc__01..10`.

Association strength is the sum over features of |partial r|, minus its
null expectation p·√(2/(πn)). Percentage reductions are reported only when the
reference strength exceeds 3 null SDs. Standardized coefficients were dropped
because they inflate with noise when an adjustment shrinks the residual variance
of the variable.

## Analyses and outputs

**A. Diversity-proportion resampling** (the draft's Figs 2–4). Samples of equal
size are drawn with 0, 25, 50, 75 and 100% non-White participants (legacy ABCD
ethno-racial variable, used only to set the sampling design).
- `ps_A_gpc_dR2_by_diversity.csv`: variance in each measure explained by
  genetic PCs 1–10 (mean, median, 75th and 99th percentile).
- `ps_A_effect_reduction_by_diversity.csv`: % drop in association strength when
  the PCs are added.
- `ps_A_gpc1_similarity_by_diversity.csv`: correlation of each variable's map
  with the genetic-PC1 map.

**B. Adjustment ladder** (full sample): cov → + family SES (household income,
highest caregiver education) → + genetic PCs → + both.
- `ps_B_adjustment_ladder.csv`: % reduction under each adjustment. The observed
  drop with genetic PCs is set against the drop expected if the PCs affected the
  brain only through the variable, 1 − √(1 − R²). It also reports the PC
  reduction after family SES (`pct_red_gpc_given_fam`) and map correlations with
  the genetic-PC1 map.
- `ps_B_structure_beyond_SES.csv`: variance in the brain measures explained by
  the genetic PCs beyond cov + income + education + COI + ADI, compared with the
  null 10/n. This is the direct test of whether population structure predicts the
  brain beyond every SES measure available.
- `ps_B_ridge_by_adjustment.csv`: out-of-sample r (discovery → replication) for
  SES, IQ and fluid IQ under each adjustment.

**C. Within-family maps** (`ps_C_within_family.csv`). Sibling and twin families
with two or more imaged children; brain measures, variable and person-level
covariates are demeaned within family. Siblings share family SES and, on average,
ancestry, so the within-family map is free of both. Reported with and without
families containing MZ twins. Only child-level variables are included, because
COI and ADI are address-level and essentially constant within family.

**D. Genetically similar subsamples** (`ps_D_genetically_similar_subsamples.csv`):
ADMIXTURE EUR ≥ 0.9 and AFR ≥ 0.5. Strength is compared with size-matched draws
from the full sample. Map similarity to the full-sample map is reported with and
without PCs. This is the ABCD analogue of Marek's mostly-white UK Biobank check.

**E. Head motion** (RSFC only).
- `ps_E_motion_selection.csv`: SES–motion correlations, and standardized
  differences in SES, genetic PC1 and IQ between children kept in and excluded
  from the RSFC sample (selection on motion).
- `ps_E_motion_ses_map.csv`: the SES map with extra motion covariates, in the
  low-motion half (with a size-matched comparison), and in the low-motion half
  with genetic PCs.

## What this cannot separate

Genetic PCs carry ancestry-patterned environment (segregation, discrimination,
neighborhood composition) as well as genetic and anthropometric variation. B and
D bound the confounding; C is the design that removes shared family environment
and, on average, shared ancestry. A local-ancestry analysis in admixed children
would separate genetic ancestry from social race further, but it needs genotype
data that are not in the tabulated release.

Tested end to end on synthetic data with the real 7.0 variable names. Not yet
run on real data.

## Extensions (November 2026 review): `popstrat_ext.R`

Run after `popstrat_main.R`, from `analyses/R`:

```r
source("scripts/popstrat/popstrat_ext.R")
```

Files: `px_data.R` (adds intracranial volume, any genetic PCs beyond 10,
scanner serial number if available; cache `dir_derived/popstrat_ext_data.rds`),
`px_analysis.R` (analyses and figures). Outputs in `output/popstrat_ext/`.
Defaults: 100 diversity resamples, 100 permutations per resample, 1,000
permutations and 200 family-bootstrap replicates in the full sample, 100
size-matched draws. Expect roughly 30-60 minutes.

| Output | Review item | What it answers |
|---|---|---|
| `px_X1_map_similarity_by_diversity.csv`, `px_X1_map_similarity_full.csv` | 2, 6 | Is a variable's map similar to the PC1 map beyond what is expected if PC1 reached the brain only through that variable? `sim_unadj` is the draft's metric; `sim_expect` is its value simulated under "PC1 acts only through the variable"; `sim_adj` correlates the variable's map with the map of PC1 residualized on the variable, with a permutation null (`sim_adj_null_q975`, `sim_adj_perm_p`) and, in the full sample, family-bootstrap 95% intervals (normal, from the bootstrap SD). `z_gadj` is the strength of the residualized-PC1 map in null SDs. |
| `px_X2_gpc_variance_vs_head_body_size.csv` | 4 | Variance in IDPs explained by PCs 1-10 before and after adjusting for ICV, height and weight, by sMRI measure type and modality. |
| `px_X2_attenuation_given_head_body_size.csv` | 4 | Does the PC attenuation of COI, ADI and IQ associations survive adjustment for head and body size? |
| `px_X2_subsamples_by_measure_type.csv` | 4 | EUR >= 0.9 and AFR >= 0.5 subsamples vs size-matched draws, by sMRI measure type (ratio reported only if the reference exceeds 3 null SDs). |
| `px_X3_robustness.csv` | 5 | % reduction in association strength under PCs 1-5, 1-10, all available PCs, ADMIXTURE proportions, PCs + ADMIXTURE, and site x scanner (x serial if available). PCs computed within ancestry groups need genotypes and are not run. |
| `popstrat_7p0_figures.pdf` | 1 | 7.0 versions of draft Figs 3-5 (from `ps_A_*.csv`), plus the adjusted map-similarity panel. |

Reading `sim_adj`: a significant value means PC1 reaches the brain other than
through the variable, along a pattern resembling the variable's map. That is
the signature of confounding only when the variable is itself correlated with
PC1 (`rho_vi_gpc1`); for a variable unrelated to PC1 it says nothing about that
variable. Tested end to end on synthetic data with planted structure: with no
direct ancestry effect (task fMRI), `sim_adj` was about 0 and `sim_unadj`
matched `sim_expect`; with an ancestry effect through head size, the PC
variance in area and volume fell by about 85% after ICV adjustment and in
thickness by about 20%, as planted.

`ps_analysis.R` defaults were raised to 100 diversity resamples and 100
size-matched draws (review item 6); `popstrat_main.R` no longer treats the
designed number of non-White participants per resample as a cell count (it
had blanked the diversity-0 rows of `ps_A_gpc_dR2_by_diversity.csv`).
