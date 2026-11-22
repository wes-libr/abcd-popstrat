# Causal ladder for Marek et al.'s exposure (neighborhood SES)

Marek et al. (2026) propose that neighborhood SES (the Child Opportunity Index
social and economic domain) shapes children's brains, through arousal and sleep,
and thereby confounds brain–IQ associations. This module applies the causal
pipeline's identification ladder to that exposure, on ABCD 7.0 imaging waves 0, 2,
4 and 6.

```r
setwd("~/Dropbox/aaa_mac/research/brains/projects/abcd-analyses/analyses/R")
sc_opts_user <- list(n_boot = 50)          # quick pass; omit for the full run (200)
source("scripts/popstrat/ses_causal.R")
```

Output: `output/ses_causal/sc_ladder.csv` (one row per method × outcome ×
adjustment), `sc_sample.csv` and `sc_movers.csv`. Cache:
`dir_derived/ses_causal_panel.rds` (host-only).

## Outcomes

- **SES signature** of the brain (`sig_rs`, `sig_ct`): each child's resting-state
  network connectivity, or Destrieux thickness, projected onto the baseline SES map.
  The map is estimated in the discovery half only. These outcomes are analysed in
  the replication half only, so the outcome is defined independently of the
  estimation sample.
- **Motor/sensory RSFC** (`ms_rs`): mean connectivity among the hand, mouth,
  visual and auditory networks. A fixed, a-priori stand-in for Marek's map.
- **IQ** (`iq`, `fiq`): NIH Toolbox total and fluid composites.
- **Height**: a comparison outcome.

All variables are in baseline-SD units, so estimates are SD of outcome per SD of
COI.

## The ladder

| Rung | Method | Identifies the effect from | Removes ancestry confounding? |
|---|---|---|---|
| 1 | LME, pooled waves (family and child random intercepts) | all variation, adjusted | only with the genetic-PC adjustment |
| 3 | AIPW, top vs bottom COI tertile at baseline | baseline, doubly robust | only with the genetic-PC adjustment |
| 4 | Causal forest overlap ATE (`grf`) | baseline, ML doubly robust | only with the genetic-PC adjustment |
| 5 | CLPM | 2-year lags | no (stable confounding biases it) |
| 6 | RI-CLPM | within-person lagged change | yes (random intercepts) |
| 7 | Person fixed effects | within-child change in COI, i.e. **moves** | yes, along with every time-invariant confounder |
| 7b | First differences | change in COI → change in outcome | yes |
| 15 | Placebo: *future* COI change → current outcome change | should be ≈ 0 | — |
| 8 | Family fixed effects | siblings at different addresses (rare) | yes |

Rungs 1, 3 and 4 are run under four adjustment sets: **cov** (age, sex,
site/scanner, motion), **fam** (+ baseline income and caregiver education),
**gpc** (+ genetic PCs 1–10) and **famgpc** (+ both). Rungs 6, 7 and 7b need no
genetic adjustment, because ancestry is constant within a child. The RI-CLPM and
CLPM also report the reverse path (outcome → COI) as a falsification check.

How to read it: if neighborhood SES causes the brain pattern, the within-child
rungs (6, 7, 7b) should keep a non-zero effect and the placebo should be about 0.
If the pattern is between-family structure (population structure, family
resources), the effect should collapse from rung 1 to rungs 6–7.

## Limits

- **Moves are not random.** Families move for reasons (income shocks, family
  change) that can affect children directly. The future-move placebo checks
  whether outcomes were already on a different trend before the move.
- **Effects build slowly.** Within-child estimates capture effects over 2-year
  windows; neighborhood exposure accumulated before baseline is absorbed into the
  child's fixed effect.
- **Imaging coverage thins at waves 4 and 6**, and NIH Toolbox composites are not
  given at every wave. `sc_sample.csv` reports n by wave, one count per row so
  suppression blanks only the small cell.

## ABCD 7.0: the exposure is baseline-only (checked 2026-11-21)

All 1,242 linked-external variables in the 7.0 dictionary (`le_*` matching
coi/adi/addr), COI and ADI included, have values at `ses-00A` only. No child
has COI at two waves, so there are no movers and the within-child rungs
(5 CLPM, 6 RI-CLPM, 7 person FE, 7b first differences, 15 future-move placebo)
cannot be estimated. `sc_run()` now detects this and writes those rungs as
"not identified" rows; rung 1 is relabelled as a baseline model. The code runs
the full ladder unchanged if a later release adds follow-up address data.

Results of the full run (`n_boot = 200`, 2026-11-21), effect of COI in
baseline-SD units, adjusting for covariates / + family SES / + genetic PCs / + both:

- IQ, LME: 0.33 / 0.13 / 0.17 / 0.05; AIPW tertile contrast 0.68 / 0.20 / 0.29 / 0.00
  [-0.10, 0.10]. Fluid IQ the same (AIPW fully adjusted -0.03).
- Resting-state SES signature (replication half), LME: 0.22 / 0.13 / 0.08 / 0.04
  (about 2 SE); AIPW fully adjusted 0.15 [-0.01, 0.30], p = 0.07; causal
  forest 0.17 (p < 0.001). The only fully adjusted signal, and the methods
  disagree on it.
- Cortical-thickness SES signature: about 0 once genetic PCs are added.
- Motor/sensory connectivity: about 0 throughout.
- Height (negative control): -0.05 with covariates only, about 0 after genetic
  PCs; family SES alone does not remove it.
- Family FE: COI differs between siblings in only 18-57 families; not usable.

Earlier tests used synthetic data with the real 7.0 variable names.
