# Endemic-epidemic space-time model of daily COVID-19 cases by ABS, Catalunya

Built from the model of Slater, Brown, Rosenthal and Mateu (2025), *Leveraging cellphone-derived mobility networks to
assess COVID-19 travel risk*, AOAS 19(1), 56-74. Their code and data: `mobilityCovid/` (clone of
github.com/JustinJSlater/mobilityCovid, 4364ecc).

Render from this folder; Spain-data is expected at `../Spain-data`.

| File | What it does | Output |
|---|---|---|
| `0_reproduce.qmd` | `stan/weekly.stan` (their model 1 in our names) on their Castilla-León data | `outputs/reproduce_castilla.rds` |
| `1_barcelona.qmd` | `stan/daily.stan`, the replication on daily cases of the 67 Barcelona city ABS | `outputs/daily_barcelona_<from>_<to>_<lag>.rds` |
| `2_catalunya.qmd` | the same for all 374 ABS, 2020-03-22 to 2021-05-01 | `outputs/daily_catalunya_2020-03-22_2021-05-01_<lag>.rds` |
| `3_figures.qmd`, `4_results.qmd` | their figures and a summary of the replication fit | rendered html |
| `5_spacetime.qmd` | `stan/spacetime.stan`, the final model (`MODEL` in `helpers.R`) or one variant of it | `outputs/spacetime_<areas>_<from>_<to>_<variant>.rds` and `_loo.rds` |
| `6_cv.qmd` | cross-validation of a set of variants (`CV_SETS`): whole weeks held out, or a random share of the ABS of each day | `outputs/cv/` |
| `7_checks.qmd` | PIT histograms and posterior predictive checks of one fit | rendered html |
| `8_recovery.qmd` | fake-data recovery of the main quantities | `outputs/recovery/` |

The final model: daily negative binomial counts with a size per ABS; an endemic term proportional to population; an
epidemic part of autoregressive, first-order neighbour and mobility terms (trip-weighted prevalence at the origins, two
7-day trip windows), lagged 7 to 13 days and scaled by a weekly transmission factor (common RW1 plus a type IV ABS x
week interaction); weekday and holiday reporting factors; ABS coefficients as random effects with static covariates.

`helpers.R` builds the inputs from Spain-data (cases, lagged pressures, trips, holidays, covariates), holds the switches
of `stan/spacetime.stan` (`SWITCHES`), the final model (`MODEL`) and the cross-validation sets. Large fits run on PALMA
with the scripts in `palma/`, which are not part of the repository.
