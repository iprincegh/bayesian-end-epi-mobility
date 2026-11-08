# Bayesian endemic-epidemic model with mobility: replication of Slater et al. (2025) for Catalunya

A daily replication of model 1 of Slater et al. (2025), an endemic-epidemic model with an exogenous, an
autoregressive, a neighbourhood and a mobility component, fitted in Stan to the 374 basic health areas (ABS) of
Catalunya from 2020-03-22 to 2021-05-01, and to the 67 ABS of Barcelona city.

**Results:** [summary of the Catalunya fit](https://iprincegh.github.io/bayesian-end-epi-mobility/4_results.html) ·
[all figures](https://iprincegh.github.io/bayesian-end-epi-mobility/3_figures.html)

## Run

This branch holds the model for the article. The replication notebooks and their pages are on `main`. Here,
`1_barcelona.qmd` and `2_catalunya.qmd` fit the replicated model (`stan/daily.stan`), the comparison model. Render from
this folder; the data repository below is expected at `../Spain-data`.

## Space-time model (this branch, in progress)

An extension of the replicated model, fitted with `stan/spacetime.stan`: daily negative binomial counts with a size
per ABS; an endemic term proportional to population; an epidemic part of autoregressive, first-order neighbour and
mobility terms (trip-weighted prevalence at the origins, two 7-day trip windows), lagged 7 to 13 days and scaled by a
weekly transmission factor (common RW1 plus a type IV ABS x week interaction); weekday and holiday reporting factors;
ABS coefficients as random effects with static covariates. The final model is `MODEL` in `helpers.R`.

| Notebook | What it does |
|---|---|
| `5_spacetime.qmd` | fits the final model, or one variant of it (switches as `SWITCHES` in `helpers.R`) |
| `6_cv.qmd` | cross-validation of a set of variants (`CV_SETS`): whole weeks held out, or a random share of the ABS of each day |
| `7_checks.qmd` | PIT histograms and posterior predictive checks of one fit |
| `8_recovery.qmd` | fake-data recovery of the main quantities |

Results of this model are not published yet.

## Data

- Daily cases and populations per ABS: the panel of the
  [covid19-mobility-healthzones v1.0 release](https://github.com/iprincegh/covid19-mobility-healthzones/releases/tag/v1.0).
- Daily trips between ABS: Niraula, P., & Oppong Boakye, P. (2026). *Mobility data from MITMA, Spain reprojected at
  Zona Básica de Salud (ZBS) level for Madrid and Àrea Bàsica de Salut (ABS) for Catalunya* (Version 1.0) [Data set].
  Zenodo. [doi:10.5281/zenodo.23737062](https://doi.org/10.5281/zenodo.23737062)

## Method

Slater, J. J., Brown, P. E., Rosenthal, J. S., & Mateu, J. (2025). Leveraging cellphone-derived mobility networks to
assess COVID-19 travel risk. *The Annals of Applied Statistics*, 19(1), 56–74.
[doi:10.1214/24-AOAS1914](https://doi.org/10.1214/24-AOAS1914). Their code:
[JustinJSlater/mobilityCovid](https://github.com/JustinJSlater/mobilityCovid) (not included here).
