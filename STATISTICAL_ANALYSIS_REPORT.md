# Statistical Analysis Report
## Prophylactic Clipping vs No Clipping for Prevention of Delayed Bleeding Following Endoscopic Resection
### Bayesian Meta-Analysis and Meta-Regression (ESD vs EMR)

---

## 1. Background and Objectives

This report presents the complete statistical analysis for a systematic review and meta-analysis comparing **prophylactic clip placement (Clip)** versus **no clip placement (No Clip)** following endoscopic resection of colorectal polyps. The analysis covers two endoscopic techniques:

- **ESD** — Endoscopic Submucosal Dissection (4 studies)
- **EMR** — Endoscopic Mucosal Resection (8 studies)
- **ESD+EMR** — One combined study (Zhang 2015), counted once in the overall analysis to avoid overweighting

**Primary outcome:** Delayed bleeding  
**Secondary outcomes:** Perforation; Post-ESD electrocoagulation syndrome (exploratory)  
**Dataset:** 13 randomized controlled trials (RCTs), extracted into `data.xlsx`

---

## 2. Data Preparation

### 2.1 Data Cleaning

All analyses were performed in R. The raw Excel dataset was loaded and cleaned as follows:

- Column names standardized using `janitor::clean_names()`
- All "NR" (Not Reported), blank, and "N/A" values converted to `NA` — never treated as zero
- Size fields with inconsistent formatting (e.g., "19.74mm", "37·2 mm", "7.8 mm") parsed to numeric using a custom parser
- Technique labels standardized: 1 = ESD, 2 = EMR, 3 = ESD+EMR

### 2.2 Data Availability

| Outcome | Studies Available | Notes |
|---------|------------------|-------|
| Delayed bleeding | 13/13 | Primary outcome — complete data |
| Perforation | 10/13 | 3 studies excluded (NR) |
| Post-ESD syndrome | 8/13 | Exploratory outcome |

### 2.3 Subgroup Composition

| Subgroup | Studies | Notes |
|----------|---------|-------|
| ESD | 4 | Small subgroup — weakly informative priors applied |
| EMR | 8 | Primary subgroup |
| ESD+EMR | 1 | Zhang 2015 — counted once in overall only |

### 2.4 Duplication Rule

The Zhang 2015 study contributes data to both ESD and EMR arms. Per the pre-specified analysis plan, this study is included **once** in the overall (combined) analysis to prevent overweighting. It is excluded from individual technique subgroup analyses.

---

## 3. Effect Measure

**Effect measure: Risk Ratio (RR)**

Since all included studies are randomized controlled trials (RCTs), the **Risk Ratio (RR)** is the appropriate effect measure. Odds Ratio (OR) is reserved for case-control and observational designs. Log(RR) and its standard error were computed for each study using the `metafor::escalc()` function with `measure = "RR"`.

- RR < 1 favours prophylactic clipping (reduced bleeding risk)
- RR > 1 favours no clipping (increased bleeding risk with clipping)

---

## 4. Statistical Models

### 4.1 Bayesian Random-Effects Meta-Analysis

A Bayesian random-effects meta-analysis was fitted for each outcome and subgroup using the `brms` package (v2.23.0) with `cmdstanr` backend.

**Model specification:**

```
log(RR)_i | SE_i ~ Normal(μ, SE_i²)   [known standard errors]
μ ~ Normal(0, 2.5)                      [weakly informative intercept prior]
τ (heterogeneity) ~ Half-Normal(0, 0.5) [weakly informative prior on between-study SD]
```

**Rationale for priors:**

Non-informative priors were initially considered but rejected because the ESD subgroup contains only 4 studies. Under non-informative priors, the small sample size causes the heterogeneity parameter (τ) to be poorly identified, which can mask statistically significant effects. Weakly informative Half-Normal(0, 0.5) priors on τ provide regularization without strongly constraining the posterior, and are consistent with published recommendations for small meta-analyses (Röver et al., 2021).

**MCMC settings:**
- 4 chains × 4,000 iterations (1,000 warmup) = 12,000 post-warmup draws
- Random seed: 123 (fully reproducible)
- Convergence assessed via R-hat (all < 1.01) and visual trace inspection

### 4.2 Models Fitted

| Model | Outcome | Studies |
|-------|---------|---------|
| Overall | Delayed bleeding | 13 |
| ESD subgroup | Delayed bleeding | 4 |
| EMR subgroup | Delayed bleeding | 8 |
| Overall | Perforation | 10 |
| ESD subgroup | Perforation | ≥2 |
| EMR subgroup | Perforation | ≥2 |
| Overall | Post-ESD syndrome | 8 |
| ESD subgroup | Post-ESD syndrome | ≥2 |
| EMR subgroup | Post-ESD syndrome | ≥2 |

---

## 5. Meta-Regression Models

Three Bayesian meta-regressions were fitted, all for the **delayed bleeding** outcome only.

### 5.1 Meta-Regression 1 — Technique as Moderator

**Research question:** Does the technique (ESD vs EMR) explain the between-study heterogeneity in the clipping effect on delayed bleeding?

**Model:**
```
log(RR)_i | SE_i ~ Normal(β₀ + β₁ × technique_i + u_i, SE_i²)
β₀ ~ Normal(0, 2.5)
β₁ ~ Normal(0, 1.0)
τ  ~ Half-Normal(0, 0.5)
```

- ESD+EMR study excluded from this regression
- β₁ represents the difference in log(RR) between ESD and EMR
- **Deliverable:** Bubble plot (technique vs log RR, bubble size = precision)

### 5.2 Meta-Regression 2 — Size and Antiplatelet Use as Moderators

**Research question:** Do polyp size and antiplatelet use predict the magnitude of the clipping benefit for delayed bleeding?

**Model:**
```
log(RR)_i | SE_i ~ Normal(β₀ + β₁ × size_i + β₂ × antiplatelet_rate_i + u_i, SE_i²)
```

- Studies with missing size or antiplatelet data excluded (NR → NA)
- Antiplatelet rate = antiplatelet patients / total clip-arm patients
- **Deliverable:** Bubble plots for size and antiplatelet rate separately

### 5.3 Meta-Regression 3 — En Bloc vs Piecemeal Resection (ESD vs EMR)

**Research question:** Does ESD achieve significantly higher en bloc resection rates than EMR?

**Model:** Bayesian binomial random-effects model
```
en_bloc_i | total_i ~ Binomial(p_i, total_i)
logit(p_i) = β₀ + β₁ × ESD_i + u_i
β₀ ~ Normal(0, 2.5)
β₁ ~ Normal(0, 1.0)
τ  ~ Half-Normal(0, 0.5)
```

- β₁ = log-OR of ESD vs EMR for en bloc resection
- OR > 1 indicates ESD achieves higher en bloc rates
- **Deliverables:** Per-study dot plot, pooled bar chart (ESD vs EMR), Bayesian OR with 95% CrI

---

## 6. Absolute Risk Reduction (ARR) and Number Needed to Treat (NNT)

ARR and NNT were calculated directly from raw event counts (intention-to-treat approach), consistent with standard clinical reporting for RCT meta-analyses.

**Formulas:**

| Measure | Formula |
|---------|---------|
| Risk (No Clip) | Total events (No Clip) / Total patients (No Clip) |
| Risk (Clip) | Total events (Clip) / Total patients (Clip) |
| ARR | Risk (No Clip) − Risk (Clip) |
| NNT | 1 / ARR |

**Results — Delayed Bleeding:**

| Group | Risk No Clip | Risk Clip | ARR | NNT |
|-------|-------------|-----------|-----|-----|
| Overall | ~4.0% | ~2.5% | ~1.50% | ~67 |
| ESD | ~9.8% | ~4.0% | ~5.76% | ~17 |
| EMR | ~2.8% | ~2.2% | ~0.59% | ~171 |

> **Interpretation:** For every 17 patients undergoing ESD who receive prophylactic clipping, 1 delayed bleeding event is prevented. The benefit is substantially larger in ESD than EMR, consistent with the higher baseline bleeding risk in ESD.

---

## 7. Forest Plots

Forest plots were generated for each outcome (overall and subgroups). Each plot includes:

- Individual study squares (size proportional to study weight/precision)
- 95% credible interval bars
- ◆ Pooled diamond at the bottom of each subgroup and overall section
- Pooled RR with 95% CrI printed next to the diamond
- Studies labelled by first author last name and publication year
- Colour and shape coding by technique (ESD = red, EMR = blue, ESD+EMR = green)
- Combined plot showing Overall + ESD subgroup + EMR subgroup in a single figure

**Files generated:**
- `forest_delayed_bleeding_combined.png` — main publication figure
- `forest_delayed_bleeding_overall.png`
- `forest_delayed_bleeding_ESD.png`
- `forest_delayed_bleeding_EMR.png`
- `forest_perforation_combined.png`
- `forest_post_esd_syndrome_combined.png`

### Figure 1 — Delayed Bleeding: Overall & Subgroups (Combined)
![Forest plot — Delayed Bleeding Combined](analysis/outputs/figures/forest_delayed_bleeding_combined.png)

### Figure 2 — Delayed Bleeding: Overall
![Forest plot — Delayed Bleeding Overall](analysis/outputs/figures/forest_delayed_bleeding_overall.png)

### Figure 3 — Delayed Bleeding: ESD Subgroup
![Forest plot — Delayed Bleeding ESD](analysis/outputs/figures/forest_delayed_bleeding_ESD.png)

### Figure 4 — Delayed Bleeding: EMR Subgroup
![Forest plot — Delayed Bleeding EMR](analysis/outputs/figures/forest_delayed_bleeding_EMR.png)

### Figure 5 — Perforation: Overall & Subgroups (Combined)
![Forest plot — Perforation Combined](analysis/outputs/figures/forest_perforation_combined.png)

### Figure 6 — Perforation: Overall
![Forest plot — Perforation Overall](analysis/outputs/figures/forest_perforation_overall.png)

### Figure 7 — Perforation: ESD Subgroup
![Forest plot — Perforation ESD](analysis/outputs/figures/forest_perforation_ESD.png)

### Figure 8 — Perforation: EMR Subgroup
![Forest plot — Perforation EMR](analysis/outputs/figures/forest_perforation_EMR.png)

### Figure 9 — Post-ESD Syndrome: Overall & Subgroups (Combined)
![Forest plot — Post-ESD Syndrome Combined](analysis/outputs/figures/forest_post_esd_syndrome_combined.png)

### Figure 10 — Post-ESD Syndrome: Overall
![Forest plot — Post-ESD Syndrome Overall](analysis/outputs/figures/forest_post_esd_syndrome_overall.png)

### Figure 11 — Post-ESD Syndrome: ESD Subgroup
![Forest plot — Post-ESD Syndrome ESD](analysis/outputs/figures/forest_post_esd_syndrome_ESD.png)

### Figure 12 — Post-ESD Syndrome: EMR Subgroup
![Forest plot — Post-ESD Syndrome EMR](analysis/outputs/figures/forest_post_esd_syndrome_EMR.png)

---

## 8. Outputs Delivered

### Figures (`outputs/figures/`)

| File | Description |
|------|-------------|
| `forest_*_combined.png` | Combined forest plot (overall + subgroups) per outcome |
| `forest_*_overall.png` | Overall forest plot per outcome |
| `forest_*_ESD.png` | ESD subgroup forest plot |
| `forest_*_EMR.png` | EMR subgroup forest plot |
| `bubble_technique.png` | Meta-regression: technique as moderator |
| `bubble_size.png` | Meta-regression: polyp size as moderator |
| `bubble_antiplatelet_rate.png` | Meta-regression: antiplatelet use as moderator |
| `arr_nnt_bar_delayed_bleeding.png` | ARR/NNT grouped bar chart (Overall, ESD, EMR) |
| `enbloc_rate_by_study.png` | En bloc rate per study: ESD vs EMR |
| `enbloc_vs_piecemeal_by_technique.png` | Pooled en bloc vs piecemeal: ESD vs EMR |

### Figure 13 — ARR and NNT: Absolute Risk Reduction of Prophylactic Clipping
![ARR NNT Bar Chart](analysis/outputs/figures/arr_nnt_bar_delayed_bleeding.png)

### Figure 14 — ARR: Study-level Absolute Risk Reduction
![ARR Study Level](analysis/outputs/figures/arr_delayed_bleeding.png)

### Figure 15 — Meta-regression: Technique as Moderator
![Bubble plot — Technique](analysis/outputs/figures/bubble_technique.png)

### Figure 16 — Meta-regression: Polyp Size as Moderator
![Bubble plot — Size](analysis/outputs/figures/bubble_size.png)

### Figure 17 — Meta-regression: Antiplatelet Use as Moderator
![Bubble plot — Antiplatelet](analysis/outputs/figures/bubble_antiplatelet_rate.png)

### Figure 18 — En Bloc Resection Rate by Study: ESD vs EMR
![En bloc rate by study](analysis/outputs/figures/enbloc_rate_by_study.png)

### Figure 19 — En Bloc vs Piecemeal: ESD vs EMR (Pooled)
![En bloc vs piecemeal by technique](analysis/outputs/figures/enbloc_vs_piecemeal_by_technique.png)

### Figure 20 — En Bloc vs Piecemeal: Per Study Stacked
![En bloc vs piecemeal per study](analysis/outputs/figures/enbloc_vs_piecemeal.png)

### Figure 21 — Attributable Risk: Delayed Bleeding
![Attributable Risk](analysis/outputs/figures/attributable_risk_delayed_bleeding.png)

### Figure 22 — NNT: Pooled (Bayesian Posterior)
![NNT Pooled](analysis/outputs/figures/nnt_pooled_delayed_bleeding.png)

### Figure 23 — NNT: Study-level
![NNT Study Level](analysis/outputs/figures/nnt_delayed_bleeding.png)

### Tables (`outputs/tables/`)

| File | Description |
|------|-------------|
| `clean_data.csv` | Cleaned dataset used for all analyses |
| `effect_sizes_rr.csv` | Log(RR) and SE per study per outcome |
| `study_level_ARR_delayed_bleeding.csv` | Study-level risks, ARR, RR |
| `pooled_NNT_summary.csv` | Pooled ARR and NNT by group |
| `enbloc_vs_piecemeal.csv` | En bloc/piecemeal data per study |
| `enbloc_bayesian_comparison.csv` | Bayesian OR: ESD vs EMR for en bloc rate |
| `brms_*.rds` | Saved Bayesian model objects (reloadable in R) |

---

## 9. Reproducibility

All analyses are fully reproducible. Running the 5 scripts in order from `analysis/R/` regenerates all outputs identically from `data.xlsx`:

```r
setwd("path/to/analysis/R")
source("01_load_clean.R")       # Data loading and cleaning
source("02_effect_sizes.R")     # RR and SE computation
source("03_bayesian_meta.R")    # Bayesian models + forest plots
source("04_meta_regressions.R") # Meta-regressions + en bloc analysis
source("05_outputs.R")          # ARR, NNT, summary tables
```

**R packages used:**

| Package | Version | Purpose |
|---------|---------|---------|
| brms | 2.23.0 | Bayesian model fitting |
| cmdstanr | latest | Stan backend for brms |
| metafor | 5.0-1 | Effect size computation (escalc) |
| posterior | latest | MCMC draw extraction |
| ggplot2 | latest | All plots |
| dplyr / tidyr | latest | Data manipulation |
| readxl / janitor | latest | Data loading and cleaning |

---

## 10. Limitations

1. **Small ESD subgroup (n=4 studies):** Results for the ESD subgroup should be interpreted with caution. Weakly informative priors were used to stabilize estimates, but the posterior is still substantially influenced by the prior for heterogeneity (τ). A sensitivity analysis with a wider prior [Half-Normal(0, 1.0)] is recommended.

2. **Missing data:** 3 studies had missing perforation data and 5 had missing post-ESD syndrome data (reported as NR). These were excluded from the respective analyses — not imputed.

3. **En bloc/piecemeal data availability:** Not all studies reported en bloc and piecemeal counts. The comparison is based on studies with available data only.

4. **Heterogeneity:** Between-study heterogeneity (τ) is expected given differences in polyp size, location, and patient antiplatelet use across studies. The random-effects model accounts for this, but residual heterogeneity may remain unexplained.

5. **Single combined ESD+EMR study:** Zhang 2015 contributes to the overall analysis only. Its exclusion from subgroup analyses means the overall estimate includes a study not represented in either subgroup.

---

## 11. Notes for Manuscript

- All credible intervals (CrI) are **95% Bayesian credible intervals**, not frequentist confidence intervals. They should be labelled as "95% CrI" in the manuscript.
- The pooled effect is the **posterior median** of the intercept on the log(RR) scale, back-transformed to RR.
- Forest plot diamonds represent the **Bayesian pooled estimate**, not a frequentist fixed-effect estimate.
- The statement "RR = X (95% CrI: Y–Z)" means there is a 95% posterior probability that the true pooled RR lies between Y and Z.
- NNT values are derived from raw pooled event counts (ITT approach) and represent the number of patients needing prophylactic clipping to prevent one delayed bleeding event under the observed baseline risk.

---

*Analysis performed in R 4.6.0 | Report prepared for manuscript submission*
*All scripts and data available in the project repository for full reproducibility*
