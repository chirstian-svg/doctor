library(dplyr)
library(ggplot2)
library(brms)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
es <- readRDS(file.path(tbl_dir, "effect_sizes_or.rds"))

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
    file = file.path(tbl_dir, paste0("brms_", model_name))
  )
}

plot_forest <- function(df, title, file) {
  # Outcome label for legend
  outcome_label <- unique(df$outcome)
  outcome_pretty <- dplyr::case_when(
    outcome_label == "delayed_bleeding"  ~ "Delayed Bleeding",
    outcome_label == "perforation"       ~ "Perforation",
    outcome_label == "post_esd_syndrome" ~ "Post-ESD Electrocoagulation Syndrome",
    TRUE ~ outcome_label
  )

  df |>
    mutate(
      or  = exp(yi),
      lo  = exp(yi - 1.96 * sei),
      hi  = exp(yi + 1.96 * sei),
      label = if ("study_label" %in% names(df)) study_label else paste0("Study ", study_id)
    ) |>
    ggplot(aes(y = reorder(label, yi), x = or, color = technique, shape = technique)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_color_manual(
      name = "Technique",
      values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087"),
      drop = FALSE
    ) +
    scale_shape_manual(
      name = "Technique",
      values = c("ESD" = 16, "EMR" = 17, "ESD+EMR" = 15),
      drop = FALSE
    ) +
    labs(
      title    = title,
      subtitle = paste0("Outcome: ", outcome_pretty, "  |  Effect measure: Odds Ratio (log scale)"),
      x        = "Odds Ratio (log scale)  [OR < 1 favours Clip]",
      y        = "Study (First Author, Year)",
      caption  = "Points show OR; horizontal lines show 95% CI.\nDashed line at OR = 1 indicates no effect."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.caption  = element_text(size = 8,  color = "grey50"),
      legend.position = "bottom",
      legend.title  = element_text(face = "bold")
    )

  ggsave(filename = file, width = 9, height = max(5, nrow(df) * 0.5 + 2), dpi = 200)
}

outcomes <- unique(es$outcome)
for (oc in outcomes) {
  df_oc <- es |> filter(outcome == oc, overall_inclusion)

  if (nrow(df_oc) < 2) next

  # Overall
  fit_meta(df_oc, model_name = paste0(oc, "_overall"), weakly_informative_tau = TRUE)
  plot_forest(
    df_oc,
    title = paste0("Clip vs No Clip — ", gsub("_", " ", tools::toTitleCase(oc)), " — All Studies (Overall)"),
    file = file.path(fig_dir, paste0("forest_", oc, "_overall.png"))
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
      title = paste0("Clip vs No Clip — ", gsub("_", " ", tools::toTitleCase(oc)), " — ", tech, " Subgroup"),
      file = file.path(fig_dir, paste0("forest_", oc, "_", gsub("[^A-Za-z0-9]+", "_", tech), ".png"))
    )
  }

  # Combined: overall + subgroups in one plot
  techs_available <- sort(unique(df_oc$technique))
  subgroup_dfs <- lapply(techs_available, function(tech) {
    d <- df_oc |> filter(technique == tech)
    if (nrow(d) < 2) return(NULL)
    d |> mutate(subgroup = paste0(tech, " subgroup"))
  })
  subgroup_dfs <- Filter(Negate(is.null), subgroup_dfs)

  if (length(subgroup_dfs) > 0) {
    df_combined <- bind_rows(
      df_oc |> mutate(subgroup = "Overall"),
      bind_rows(subgroup_dfs)
    ) |> mutate(subgroup = factor(subgroup, levels = c("Overall", paste0(techs_available, " subgroup"))))

    outcome_pretty_combined <- dplyr::case_when(
      oc == "delayed_bleeding"  ~ "Delayed Bleeding",
      oc == "perforation"       ~ "Perforation",
      oc == "post_esd_syndrome" ~ "Post-ESD Electrocoagulation Syndrome",
      TRUE ~ oc
    )

    label_col <- if ("study_label" %in% names(df_combined)) "study_label" else "study_id"

    df_combined |>
      mutate(
        or    = exp(yi),
        lo    = exp(yi - 1.96 * sei),
        hi    = exp(yi + 1.96 * sei),
        label = if (label_col == "study_label") study_label else paste0("Study ", study_id)
      ) |>
      ggplot(aes(y = reorder(label, yi), x = or, color = technique, shape = technique)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.6) +
      geom_point(size = 3) +
      scale_x_log10() +
      scale_color_manual(
        name = "Technique",
        values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087"),
        drop = FALSE
      ) +
      scale_shape_manual(
        name = "Technique",
        values = c("ESD" = 16, "EMR" = 17, "ESD+EMR" = 15),
        drop = FALSE
      ) +
      facet_wrap(~ subgroup, scales = "free_y", ncol = 1) +
      labs(
        title    = paste0("Clip vs No Clip — ", outcome_pretty_combined, " — Overall & Subgroups"),
        subtitle = "Effect measure: Odds Ratio (log scale)",
        x        = "Odds Ratio (log scale)  [OR < 1 favours Clip]",
        y        = "Study (First Author, Year)",
        caption  = "Points show OR; horizontal lines show 95% CI.\nDashed line at OR = 1 indicates no effect."
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10, color = "grey40"),
        plot.caption    = element_text(size = 8,  color = "grey50"),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold"),
        strip.text      = element_text(face = "bold", size = 11)
      )

    ggsave(
      filename = file.path(fig_dir, paste0("forest_", oc, "_combined.png")),
      width = 10,
      height = max(8, nrow(df_combined) * 0.45 + 3),
      dpi = 200
    )
  }
}

message("Bayesian meta-analysis models saved under outputs/tables/. Forest plots saved to outputs/figures/.")

