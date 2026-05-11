library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
es  <- readRDS(file.path(tbl_dir, "effect_sizes_or.rds"))
dat <- readRDS(file.path(tbl_dir, "clean_data.rds"))

# Join moderators
mod <- dat |>
  transmute(
    study_id,
    technique,
    technique_code,
    avg_size,
    antiplatelet_clip
  )

df <- es |>
  filter(outcome == "delayed_bleeding") |>
  left_join(mod, by = c("study_id", "technique"))

if (nrow(df) < 2) stop("Not enough delayed bleeding studies for meta-regression.")

# Technique moderator (ESD vs EMR) — exclude ESD+EMR unless explicitly desired
df_tech <- df |> filter(technique %in% c("ESD", "EMR"))

fit_tech <- brm(
  yi | se(sei) ~ 1 + technique + (1 | study_id),
  data = df_tech,
  family = gaussian(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "Intercept"),
    set_prior("normal(0, 1.0)", class = "b"),
    set_prior("normal(0, 0.5)", class = "sd")
  ),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 123,
  backend = "cmdstanr",
  file = file.path(tbl_dir, "brms_meta_reg_technique")
)

# Size (+ antiplatelet) moderators
df_size <- df |> filter(!is.na(avg_size))
df_size <- df_size |>
  mutate(
    antiplatelet_rate = antiplatelet_clip / n_clip,
    antiplatelet_rate = ifelse(is.finite(antiplatelet_rate), antiplatelet_rate, NA_real_)
  )

fit_size <- brm(
  yi | se(sei) ~ 1 + avg_size + antiplatelet_rate + (1 | study_id),
  data = df_size,
  family = gaussian(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "Intercept"),
    set_prior("normal(0, 1.0)", class = "b"),
    set_prior("normal(0, 0.5)", class = "sd")
  ),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 123,
  backend = "cmdstanr",
  file = file.path(tbl_dir, "brms_meta_reg_size_antiplatelet")
)

# Bubble plots (approximate): yi vs moderator with size by precision
bubble <- function(d, x, xlab, file) {
  d <- d |> mutate(weight = 1 / (sei^2))
  ggplot(d, aes(x = .data[[x]], y = yi, size = weight, color = technique)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.75) +
    scale_size_continuous(range = c(2, 10)) +
    labs(x = xlab, y = "log(OR) for delayed bleeding", title = "Bubble plot (study-level)") +
    theme_minimal(base_size = 12)
  ggsave(file, width = 8, height = 5, dpi = 200)
}

bubble(df_tech, "technique_code", "Technique (1=ESD, 2=EMR)", file.path(fig_dir, "bubble_technique.png"))
bubble(df_size |> filter(!is.na(avg_size)), "avg_size", "Average polyp size (mm)", file.path(fig_dir, "bubble_size.png"))
bubble(df_size |> filter(!is.na(antiplatelet_rate)), "antiplatelet_rate", "Antiplatelet rate (clip arm)", file.path(fig_dir, "bubble_antiplatelet_rate.png"))

# ---- Piecemeal vs En Bloc analysis ----
# Pull en bloc / piecemeal columns from clean data
enbloc_cols <- c("end_block_resection_21", "piecemeal_22")
enbloc_available <- all(enbloc_cols %in% names(dat))

if (enbloc_available) {
  enbloc_dat <- dat |>
    transmute(
      study_id,
      study_label = if ("study_label" %in% names(dat)) study_label else paste0("Study ", study_id),
      technique,
      n_clip,
      en_bloc   = suppressWarnings(as.numeric(end_block_resection_21)),
      piecemeal = suppressWarnings(as.numeric(piecemeal_22)),
      total_resections = en_bloc + piecemeal,
      pct_en_bloc   = en_bloc   / total_resections * 100,
      pct_piecemeal = piecemeal / total_resections * 100
    ) |>
    filter(!is.na(pct_en_bloc), is.finite(pct_en_bloc))

  if (nrow(enbloc_dat) >= 2) {
    # Summary table
    write.csv(enbloc_dat, file.path(tbl_dir, "enbloc_vs_piecemeal.csv"), row.names = FALSE)

    # Bar chart: en bloc vs piecemeal rate per study
    enbloc_long <- enbloc_dat |>
      mutate(en_bloc_order = pct_en_bloc) |>
      tidyr::pivot_longer(
        cols = c(pct_en_bloc, pct_piecemeal),
        names_to  = "resection_type",
        values_to = "percent"
      ) |>
      mutate(resection_type = dplyr::recode(resection_type,
        pct_en_bloc   = "En Bloc",
        pct_piecemeal = "Piecemeal"
      ))

    ggplot(enbloc_long,
           aes(x = percent,
               y = reorder(study_label, en_bloc_order),
               fill = resection_type)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(
        name   = "Resection Type",
        values = c("En Bloc" = "#4DBBD5", "Piecemeal" = "#E64B35")
      ) +
      facet_wrap(~ technique, scales = "free_y", ncol = 1) +
      labs(
        title    = "En Bloc vs Piecemeal Resection by Study",
        subtitle = "Proportion of resections by type, grouped by technique",
        x        = "Percentage (%)",
        y        = "Study (First Author, Year)",
        caption  = "En Bloc = complete single-piece resection; Piecemeal = multiple fragments."
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

    ggsave(file.path(fig_dir, "enbloc_vs_piecemeal.png"), width = 9, height = max(6, nrow(enbloc_dat) * 0.5 + 3), dpi = 200)
    message("En bloc vs piecemeal plot saved.")
  } else {
    message("Not enough en bloc/piecemeal data to plot.")
  }
} else {
  message("En bloc/piecemeal columns not found in clean data — skipping.")
}

message("Meta-regressions saved under outputs/tables/. Bubble plots saved to outputs/figures/.")

