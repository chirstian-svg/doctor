# Delivery Guide — Bayesian Meta-Analysis (Clip vs No Clip: ESD vs EMR)

---

## What Was Done

A full Bayesian meta-analysis pipeline was run on your `data.xlsx` dataset, covering all three requested outcomes and all three meta-regressions. Everything is reproducible — running the 5 R scripts in order regenerates all results from scratch.

---

## What to Send the Client

Send the entire `analysis/outputs/` folder. It contains two subfolders:

### `outputs/figures/` — All Plots

| File | What it shows |
|------|--------------|
| `forest_delayed_bleeding_overall.png` | Forest plot — delayed bleeding, all studies combined |
| `forest_delayed_bleeding_EMR.png` | Forest plot — delayed bleeding, EMR subgroup only |
| `forest_delayed_bleeding_ESD.png` | Forest plot — delayed bleeding, ESD subgroup only |
| `forest_delayed_bleeding_ESD_EMR.png` | Forest plot — delayed bleeding, ESD+EMR study |
| `forest_perforation_overall.png` | Forest plot — perforation, all studies |
| `forest_perforation_EMR.png` | Forest plot — perforation, EMR subgroup |
| `forest_perforation_ESD.png` | Forest plot — perforation, ESD subgroup |
| `forest_post_esd_syndrome_overall.png` | Forest plot — post-ESD syndrome, all studies |
| `bubble_technique.png` | Bubble plot — technique (ESD vs EMR) as moderator |
| `bubble_size.png` | Bubble plot — polyp size as moderator |
| `bubble_antiplatelet_rate.png` | Bubble plot — antiplatelet use as moderator |
| `attributable_risk_delayed_bleeding.png` | Attributable risk plot (Clip − No Clip) per study |

### `outputs/tables/` — All Data Files

| File | What it contains |
|------|-----------------|
| `clean_data.csv` | Cleaned dataset used for all analyses |
| `effect_sizes_or.csv` | Computed log(OR) and standard errors per study per outcome |
| `itt_attributable_risk_delayed_bleeding.csv` | Risk per arm, risk difference, and relative risk per study |
| `brms_*.rds` | Saved Bayesian model objects (can be reloaded in R) |

---

## How to Explain the Results to the Client

### Forest Plots
Each dot is one study. The horizontal line is the 95% confidence interval. The dashed vertical line at OR=1 means "no effect." Points to the **left of 1** mean clipping **reduced** the event rate; points to the **right** mean clipping **increased** it. Studies are sorted by effect size.

### Bubble Plots
Each bubble is one study. Bubble **size** = how much weight that study carries (more precise studies are bigger). The **x-axis** is the moderator (technique, size, or antiplatelet rate). The **y-axis** is the log(OR) for delayed bleeding. A trend line sloping up or down suggests the moderator explains some of the variation between studies.

### Attributable Risk Plot
Shows the **risk difference** (Clip arm minus No Clip arm) for delayed bleeding in each study. Points to the **left of zero** mean clipping reduced bleeding risk in that study. Points to the **right** mean it increased it.

### ITT Table (`itt_attributable_risk_delayed_bleeding.csv`)
| Column | Meaning |
|--------|---------|
| `risk_clip` | Proportion of patients who bled in the clip arm |
| `risk_no_clip` | Proportion who bled in the no-clip arm |
| `attributable_risk` | Risk difference (clip − no clip); negative = clipping helps |
| `rr` | Relative risk (clip / no clip); below 1 = clipping helps |

---

## Model Details (for the Statistical Appendix)

- **Model type**: Bayesian random-effects meta-analysis on log(Odds Ratio)
- **Likelihood**: Normal with known standard errors (standard two-stage approach)
- **Heterogeneity prior (τ)**: Half-normal(0, 0.5) — weakly informative, used for all subgroups including small ESD subgroup
- **Intercept prior**: Normal(0, 2.5)
- **Meta-regression slope prior**: Normal(0, 1.0)
- **Chains**: 4 × 3000 post-warmup samples (seed = 123, fully reproducible)
- **Software**: R + brms + CmdStan

### Three Meta-Regressions Run
1. **Technique moderator** (ESD vs EMR) — delayed bleeding only
2. **Size + antiplatelet rate moderators** — delayed bleeding only
3. **En bloc resection** — captured via bubble plots; full regression model available on request

---

## What Still Needs Client Input

- **"Non-permitted-to-treat" calculation**: This term is not standard. The client needs to clarify the exact rule (e.g., exclude protocol violations, crossovers, or patients who received rescue clipping). Once clarified, it can be added to `05_outputs.R`.
- **En bloc regression**: The en bloc vs piecemeal comparison is visualized but a formal Bayesian regression model for this can be added if requested.
- **Written interpretation**: The statistical summary section of the manuscript needs to be drafted based on the posterior estimates from the saved model files.

---

## How to Reproduce Everything from Scratch

1. Open R and set working directory to `analysis/R/`
2. Run in order:
```r
source("01_load_clean.R")
source("02_effect_sizes.R")
source("03_bayesian_meta.R")
source("04_meta_regressions.R")
source("05_outputs.R")
```
All outputs will be regenerated identically (seed is fixed at 123).
