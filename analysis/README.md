# Analysis workflow (step-by-step)

This folder is the implementation of `../PROJECT_REQUIREMENTS.md`.

## What you run
From `work/doctor/analysis/` in R:

0) Prerequisite

- Install R (this environment currently doesn’t have it). On Ubuntu:

```bash
sudo apt update && sudo apt install -y r-base
```

1) Install packages (first time only)

```r
install.packages(c(
  "tidyverse", "readxl", "janitor", "stringr",
  "metafor",
  "brms", "cmdstanr",
  "posterior", "bayesplot",
  "gt", "here"
))
```

2) Run the pipeline scripts in order

```r
source("R/01_load_clean.R")
source("R/02_effect_sizes.R")
source("R/03_bayesian_meta.R")
source("R/04_meta_regressions.R")
source("R/05_outputs.R")
```

Outputs are written to:
- `outputs/figures/`
- `outputs/tables/`

## Notes
- This repo assumes `../data.xlsx` is the canonical input.
- If your Excel column names differ (common), adjust the mapping table in `R/01_load_clean.R`.

