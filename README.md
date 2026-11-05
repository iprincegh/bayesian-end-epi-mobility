# Bayesian endemic-epidemic model with mobility: replication of Slater et al. (2025) for Catalunya

A daily replication of model 1 of Slater et al. (2025), an endemic-epidemic model with an exogenous, an
autoregressive, a neighbourhood and a mobility component, fitted in Stan to the 374 basic health areas (ABS) of
Catalunya from 2020-03-22 to 2021-05-01, and to the 67 ABS of Barcelona city.

**Results:** [summary of the Catalunya fit](https://iprincegh.github.io/bayesian-end-epi-mobility/4_results.html) ·
[all figures](https://iprincegh.github.io/bayesian-end-epi-mobility/3_figures.html)

## Run

Render the notebooks from this folder, in order: `1_barcelona.qmd`, `2_catalunya.qmd`, `3_figures.qmd`,
`4_results.qmd`. They expect the data repository below at `../Spain-data`. `0_reproduce.qmd` refits the original
Castilla-León model and needs Slater's repository cloned into `mobilityCovid/`.

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
