# abcd-popstrat

Analysis code for the ABCD Release 7.0 analyses in

> Fan CC, Loughnan R, Xu B, Zheng W, Misaki M, Martin A, Schork AJ, Paulus MP, Thompson WK.
> *The unintended racialization of neuroanatomy in the pursuit of external validity for neuroimaging studies.*

The code reads the ABCD Study tabulated data (Release 7.0) through
[NBDCtools](https://cran.r-project.org/package=NBDCtools). No data are included
here: ABCD data are available to approved researchers through the
[NBDC Data Hub](https://www.nbdc-datahub.org/abcd-release-7-0)
(doi:[10.82525/8f3w-5260](https://doi.org/10.82525/8f3w-5260)).

## Setup

1. Download the ABCD 7.0 tabulated release and set the environment variable
   `ABCD_DATA_ROOT` to the folder that contains it (see `00_utils.R` for the
   expected layout, including the ADMIXTURE ancestry file).
2. Run everything from the repository root (the folder that contains
   `00_utils.R`). Participant-level caches are written to the derived-data
   folder under `ABCD_DATA_ROOT`, never inside this repository;
   `00_utils.R` refuses to write `.rds`/`.RData` files inside it.
3. Every output table is aggregate, and any result based on fewer than 10
   participants is suppressed (`safe_write_csv()` in `00_utils.R`).

## Scripts (`scripts/popstrat/`)

| Runner | Analyses | Manuscript |
|---|---|---|
| `marek_replication.R` | Replication of Marek et al. (2026) with and without genetic PCs (`mr_*.R`) | Table S3 |
| `popstrat_main.R` | Diversity resampling, adjustment ladder, within-family, genetically similar subsamples, head motion (`ps_*.R`) | Figs 3-4, 6; Tables S4-S5 |
| `popstrat_ext.R` | Adjusted map similarity with permutation nulls, head and body size, robustness of the PC adjustment, 7.0 figures (`px_*.R`) | Fig 5; Tables S7-S8 |
| `popstrat_supp.R` | Genetic PC distributions and VI-PC associations (`px_supp.R`) | Figs S1, S3 |
| `ses_causal.R` | Causal analyses of neighborhood SES (`sc_*.R`) | Table S6 |

Run order: `marek_replication.R`, then `popstrat_main.R`, then
`popstrat_ext.R`, `popstrat_supp.R` and `ses_causal.R`. Each script's header
and the `README_*.md` files describe options, outputs and run times.

R packages: NBDCtools, nlme, lavaan, grf (optional), plus base R.
