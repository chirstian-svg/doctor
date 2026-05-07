# Project Requirements (Freelance) — Bayesian Meta-analysis + Meta-regression (ESD vs EMR)

## 1) Background / Goal

The client needs a **Bayesian meta-analysis** of clinical trial data comparing **prophylactic intervention** (Clip vs No Clip) across **two endoscopic techniques**:

- **ESD** (Endoscopic Submucosal Dissection)
- **EMR** (Endoscopic Mucosal Resection)
- One study includes **ESD+EMR** and must be handled to avoid overweighting in subgroup analyses.

Primary interest is **side effects**, with focus on **delayed bleeding** as the main endpoint.

## 2) Stakeholders

- **Client / project lead**: provides extracted Excel dataset, manuscript context, submission timeline.
- **Freelancer (you)**: performs Bayesian models, meta-regressions, diagnostics, plots/tables, reproducible code & handoff.

## 3) Inputs (what is provided)

- **Excel dataset**: `data.xlsx` (study-level rows; columns include events, sample sizes, moderators).
- Supporting notes:
  - `projectdetail.md` (job post)
  - `chat.md` (chat history, timing expectations, priors note)
  - `call.md` (analysis plan and deliverables)

## 4) Outcomes / Endpoints

### 4.1 Primary outcome

- **Delayed bleeding** (full availability noted as 13/13 studies in chat).

### 4.2 Secondary outcomes

- **Perforation** (noted as 10/13 studies in chat).
- **Post-ESD electrocoagulation syndrome** (noted as 8/13 studies in chat; exploratory).

## 5) Analyses Required (Functional Requirements)

### 5.1 Core meta-analyses (Bayesian)

- Run **Bayesian random-effects meta-analysis** for:
  - Overall (all eligible studies)
  - Subgroups by technique (ESD vs EMR; ESD+EMR handled as specified in §5.4)
- Provide posterior estimates (effect size + uncertainty) and heterogeneity summaries.

### 5.2 Meta-regressions (Bayesian)

Per call notes, **three meta-regressions** are required:

1. **Technique as moderator** (main meta-regression)
  - Moderator: technique (ESD vs EMR)
  - Outcome: **delayed bleeding** only
  - Goal: quantify how much of the effect difference is explained by technique.
2. **Size and antiplatelet moderators** (risk-factor meta-regression)
  - Moderator(s): **size** and **size + antiplatelet use**
  - Outcome: delayed bleeding
  - Deliverable: bubble plots + regression summary.
3. **En bloc resection comparison** (technique comparison)
  - Compare en bloc resection rates between techniques.
  - Source fields referenced in call: “U and B” columns / en bloc vs piecemeal indicators.
  - Deliverable: model output + plot/table.

### 5.3 “Regular” meta-analyses for other outcomes

For outcomes other than delayed bleeding (per call):

- Perforation: Bayesian meta-analysis (no technique meta-regression required unless requested later).
- Post-ESD syndrome: Bayesian meta-analysis (exploratory).

### 5.4 Study weighting / duplication rule

- If a single article contributes data to both technique arms (ESD+EMR), it should be counted **once in the combined analysis** (to avoid overweighting), while still allowing the planned subgroup/overall strategy described by the client.
  - **Implementation detail**: define an explicit rule in code and document it in the analysis report.

### 5.5 Priors & small-subgroup stability

- Initial plan considered non-informative priors, but client observed ESD subgroup (≈4 studies) is unstable under non-informative priors.
- Requirement:
  - Use **weakly informative priors** for small subgroups (e.g., half-normal / half-t for heterogeneity), with **sensitivity checks**.
  - Document priors, rationale, and show that results are not artifacts of a single prior choice.

### 5.6 Additional effect measures / calculations

Client requested:

- **Intention-to-treat** calculations
- **Non-permitted-to-treat** calculations (as requested wording)
- **Attributable risk** calculations
Deliverables include both numeric outputs and a **graphic** summarizing these measures.

## 6) Outputs / Deliverables

### 6.1 Reproducible code package

- R project (or structured scripts) that run end-to-end:
  - Data loading/cleaning
  - Model fitting
  - Diagnostics
  - Plot generation
  - Table generation
- A single command/script entrypoint to reproduce results on another machine.

### 6.2 Publication-ready results

For each outcome and relevant subgroup/regression:

- **Forest plots** (Bayesian)
- **Summary tables** suitable for manuscript/statistical appendix
- **Bubble plots** for each meta-regression
- Diagnostics summary (convergence and fit checks)

### 6.3 Written analysis summary

Short report including:

- Model specifications
- Priors used + sensitivity analysis
- Interpretation of key findings (statistical portion for publication)
- Notes on limitations (e.g., small ESD subgroup, missingness)

## 7) Non-functional Requirements (Quality)

- **Correctness**: event counts and sample sizes must match the provided Excel.
- **Reproducibility**: same inputs reproduce same outputs (set seeds, record package versions).
- **Transparency**: clear documentation of:
  - inclusion/exclusion rules
  - missing data handling (e.g., “NR” treated as missing, not zero)
  - subgroup definitions
  - duplication handling for ESD+EMR study
- **Robustness**: sensitivity analyses for priors and influence of small subgroups.

## 8) Data Requirements / Field Mapping (from `data.xlsx`)

The Excel appears to contain columns for:

- Sample sizes (intervention “clip” arm; control “no clip” arm)
- Antiplatelet and anticoagulation counts
- Location categories (cecum, ascending, transverse, descending, sigmoid, rectum) — at least for control shown
- Event counts for:
  - delayed bleeding (control)
  - perforation (control)
  - post-ESD syndrome (control)
  - (and corresponding intervention fields are expected / required; if not present, must be clarified/derived)
- Morphology / histology categories (adenoma, serrated, hyperplastic; flat/non-polypoid; sessile)
- Average size (messy formatting possible, e.g., “19.74mm”, “30 mm”, “7.8 mm”, “NR”)
- Resection outcomes: en bloc, piecemeal
- Technique label + encoded technique (1=ESD, 2=EMR, 3=ESD+EMR)
- Size cut-off field (optional / may be blank)

If any required paired fields (intervention vs control) are missing for an outcome, that becomes a **blocker** for that endpoint until resolved.

## 9) Project Phases (Work Plan)

### Phase 1 — Data review & preparation

- Validate dataset schema and consistency
- Clean/parse size fields
- Standardize technique labels
- Missingness handling (NR → NA)
- Confirm duplication rule for ESD+EMR article

### Phase 2 — Modeling & visualization

- Fit Bayesian random-effects models per outcome and subgroup
- Fit meta-regressions as specified in §5.2
- Generate forest plots, bubble plots, and tables
- Run convergence + sensitivity checks

### Phase 3 — Final report & handoff

- Final written results summary + interpretation notes
- Deliver full reproducible code bundle + outputs

## 10) Timeline

Target: **~7–10 days after receiving the final, duplicated/validated dataset**, consistent with chat/call expectations and submission goal (mid-May).

## 11) Acceptance Criteria (Definition of Done)

- All requested analyses in §5 completed, with outputs in §6 produced.
- Forest plots and summary tables are publication-ready.
- Bubble plots delivered for all requested regressions.
- ITT / non-permitted-to-treat / attributable risk delivered + one combined graphic.
- Clear documentation of priors + at least one sensitivity analysis.
- A reviewer can run the provided scripts and regenerate the same figures/tables from `data.xlsx`.