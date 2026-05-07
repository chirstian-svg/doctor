library(dplyr)
library(ggplot2)
library(brms)

es <- readRDS("outputs/tables/effect_sizes_or.rds")

# Duplication / overweighting rule placeholder:
# If a single paper contributes multiple technique arms, ensure it is counted once for "overall".
# Current Excel appears to be one row per study; if not, add a `paper_id` column in 01_load_clean.R
es <- es |> mutate(overall_inclusion = TRUE)

fit_meta <- function(df, model_name, weakly_informative_tau = TRUE) {
  # Normal likelihood meta-analysis on log(OR):
  # yi | se(sei) ~ 1
  # weakly informative prior for intercept; half-normal for residual SD is not used here since sei is known.
  priors <- c(
    set_prior("normal(0, 2.5)", class = "Intercept")
  )

  # Heterogeneity via group-level random intercept
  # Equivalent to random-effects meta-analysis on yi
  # Use weakly informative prior on sd (tau) especially for small subgroups
  if (weakly_informative_tau) {
    priors <- c(priors, set_prior("normal(0, 0.5)", class = "sd")) # half-normal implied by sd>0
  } else {
    priors <- c(priors, set_prior("student_t(3, 0, 2.5)", class = "sd"))
  }

  brm(
    yi | se(sei) ~ 1 + (1 | study_id),
    data = df,
    family = gaussian(),
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 1000,
    seed = 123,
    backend = "cmdstanr",
    file = file.path("outputs", "tables", paste0("brms_", model_name))
  )
}

plot_forest <- function(df, title, file) {
  df |> 
    mutate(or = exp(yi), lo = exp(yi - 1.96 * sei), hi = exp(yi + 1.96 * sei)) |>
    ggplot(aes(y = reorder(paste0("Study ", study_id), yi), x = or)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(title = title, x = "Odds Ratio (log scale)", y = NULL) +
    theme_minimal(base_size = 12)

  ggsave(filename = file, width = 8, height = 5, dpi = 200)
}

outcomes <- unique(es$outcome)
for (oc in outcomes) {
  df_oc <- es |> filter(outcome == oc, overall_inclusion)

  if (nrow(df_oc) < 2) next

  # Overall
  fit_meta(df_oc, model_name = paste0(oc, "_overall"), weakly_informative_tau = TRUE)
  plot_forest(
    df_oc,
    title = paste0("Forest plot (approx.) — ", oc, " — Overall"),
    file = file.path("outputs/figures", paste0("forest_", oc, "_overall.png"))
  )

  # Subgroups
  for (tech in sort(unique(df_oc$technique))) {
    df_t <- df_oc |> filter(technique == tech)
    if (nrow(df_t) < 2) next

    small_group <- nrow(df_t) <= 4
    fit_meta(
      df_t,
      model_name = paste0(oc, "_", gsub("[^A-Za-z0-9]+", "_", tech)),
      weakly_informative_tau = small_group
    )
    plot_forest(
      df_t,
      title = paste0("Forest plot (approx.) — ", oc, " — ", tech),
      file = file.path("outputs/figures", paste0("forest_", oc, "_", gsub("[^A-Za-z0-9]+", "_", tech), ".png"))
    )
  }
}

message("Bayesian meta-analysis models saved under outputs/tables/. Forest plots saved to outputs/figures/.")

