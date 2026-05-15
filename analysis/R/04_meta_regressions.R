library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(brms)
library(posterior)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
es  <- readRDS(file.path(tbl_dir, "effect_sizes_rr.rds"))
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

# Bubble plots: yi vs moderator with regression line and credible band
bubble <- function(d, x, xlab, title, file) {
  d <- d |> mutate(weight = 1 / (sei^2))
  p <- ggplot(d, aes(x = .data[[x]], y = yi, size = weight, color = technique)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_smooth(aes(x = .data[[x]], y = yi, weight = weight),
                method = "lm", se = TRUE,
                color = "grey30", fill = "grey80", alpha = 0.3,
                linewidth = 0.8, inherit.aes = FALSE) +
    geom_point(alpha = 0.85, shape = 16) +
    scale_size_continuous(
      name   = "Precision (1/SE\u00b2)",
      range  = c(2, 10)
    ) +
    scale_color_manual(
      name   = "Technique",
      values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087")
    ) +
    labs(
      title    = title,
      subtitle = "Outcome: Delayed Bleeding  |  Effect measure: log(RR)",
      x        = xlab,
      y        = "log(RR) for delayed bleeding  [< 0 favours Clip]",
      caption  = "Bubble size proportional to study precision (1/SE\u00b2).\nLine = weighted regression; shaded band = 95% CI."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 10, color = "grey40"),
      plot.caption    = element_text(size = 8,  color = "grey50", hjust = 0),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold")
    )
  ggsave(file, plot = p, width = 8, height = 6, dpi = 200)
}

bubble(df_tech, "technique_code",
       xlab  = "Technique (1 = ESD, 2 = EMR)",
       title = "Meta-regression of Clipping Effect by Technique",
       file  = file.path(fig_dir, "bubble_technique.png"))

bubble(df_size |> filter(!is.na(avg_size)), "avg_size",
       xlab  = "Average polyp size (mm)",
       title = "Meta-regression of Clipping Effect by Polyp Size",
       file  = file.path(fig_dir, "bubble_size.png"))

bubble(df_size |> filter(!is.na(antiplatelet_rate)), "antiplatelet_rate",
       xlab  = "Antiplatelet use rate (clip arm)",
       title = "Meta-regression of Clipping Effect by Antiplatelet Use",
       file  = file.path(fig_dir, "bubble_antiplatelet_rate.png"))

# ---- En Bloc vs Piecemeal: ESD vs EMR comparison ----
enbloc_cols <- c("end_block_resection_21", "piecemeal_22")
enbloc_available <- all(enbloc_cols %in% names(dat))

if (enbloc_available) {

  enbloc_dat <- dat |>
    transmute(
      study_id,
      study_label = paste0(
        stringr::str_to_title(stringr::str_trim(as.character(first_author_last_name))),
        " ", publication_year
      ),
      technique,
      en_bloc   = suppressWarnings(as.numeric(end_block_resection_21)),
      piecemeal = suppressWarnings(as.numeric(piecemeal_22)),
      total_resections = en_bloc + piecemeal,
      pct_en_bloc   = en_bloc   / total_resections * 100,
      pct_piecemeal = piecemeal / total_resections * 100
    ) |>
    filter(!is.na(pct_en_bloc), is.finite(pct_en_bloc),
           technique %in% c("ESD", "EMR"))

  if (nrow(enbloc_dat) >= 2) {
    write.csv(enbloc_dat, file.path(tbl_dir, "enbloc_vs_piecemeal.csv"), row.names = FALSE)

    # ── Plot 1: Per-study dot plot comparing en bloc rate by technique ──────
    ggplot(enbloc_dat,
           aes(x = pct_en_bloc,
               y = reorder(study_label, pct_en_bloc),
               color = technique, shape = technique)) +
      geom_vline(xintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_point(size = 3.5) +
      scale_color_manual(
        name   = "Technique",
        values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5")
      ) +
      scale_shape_manual(
        name   = "Technique",
        values = c("ESD" = 16, "EMR" = 17)
      ) +
      scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
      labs(
        title    = "En Bloc Resection Rate by Study: ESD vs EMR",
        subtitle = "Each point = en bloc rate for that study; dashed line at 50%",
        x        = "En Bloc Resection Rate (%)",
        y        = "Study (First Author, Year)",
        caption  = "En Bloc = complete single-piece resection.\nESD is expected to achieve higher en bloc rates than EMR."
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10, color = "grey40"),
        plot.caption    = element_text(size = 8,  color = "grey50", hjust = 0),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold"),
        legend.text     = element_text(size = 11)
      )
    ggsave(file.path(fig_dir, "enbloc_rate_by_study.png"),
           width = 9, height = max(5, nrow(enbloc_dat) * 0.5 + 2), dpi = 200)

    # ── Plot 2: Grouped bar — pooled en bloc vs piecemeal per technique ─────
    tech_summary <- enbloc_dat |>
      group_by(technique) |>
      summarise(
        total_en_bloc   = sum(en_bloc,   na.rm = TRUE),
        total_piecemeal = sum(piecemeal, na.rm = TRUE),
        total           = total_en_bloc + total_piecemeal,
        pct_en_bloc     = total_en_bloc   / total * 100,
        pct_piecemeal   = total_piecemeal / total * 100,
        .groups = "drop"
      )

    tech_long <- tech_summary |>
      tidyr::pivot_longer(
        cols      = c(pct_en_bloc, pct_piecemeal),
        names_to  = "resection_type",
        values_to = "percent"
      ) |>
      mutate(resection_type = dplyr::recode(resection_type,
        pct_en_bloc   = "En Bloc",
        pct_piecemeal = "Piecemeal"
      ))

    ggplot(tech_long,
           aes(x = technique, y = percent, fill = resection_type)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.5) +
      geom_text(
        aes(label = paste0(round(percent, 1), "%")),
        position = position_dodge(width = 0.6),
        vjust = -0.4, size = 3.5, fontface = "bold"
      ) +
      scale_fill_manual(
        name   = "Resection Type",
        values = c("En Bloc" = "#4DBBD5", "Piecemeal" = "#E64B35")
      ) +
      scale_y_continuous(
        limits = c(0, 110),
        labels = function(x) paste0(x, "%")
      ) +
      labs(
        title    = "En Bloc vs Piecemeal Resection: ESD vs EMR",
        subtitle = "Pooled resection rates across all studies per technique",
        x        = "Resection Technique",
        y        = "Resection Rate (%)",
        caption  = "En Bloc = complete single-piece resection; Piecemeal = multiple fragments.\nPooled from all studies with available resection data."
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(size = 10, color = "grey40"),
        plot.caption       = element_text(size = 8,  color = "grey50", hjust = 0),
        legend.position    = "bottom",
        legend.title       = element_text(face = "bold"),
        legend.text        = element_text(size = 11),
        axis.text.x        = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank()
      )
    ggsave(file.path(fig_dir, "enbloc_vs_piecemeal_by_technique.png"),
           width = 7, height = 6, dpi = 200)

    # ── Bayesian comparison: en bloc rate ESD vs EMR ────────────────────────
    # Model: logit(en_bloc / total) ~ technique, binomial
    enbloc_model_dat <- enbloc_dat |>
      filter(!is.na(en_bloc), !is.na(total_resections), total_resections > 0) |>
      mutate(
        technique_bin = ifelse(technique == "ESD", 1, 0),
        en_bloc_int   = as.integer(round(en_bloc)),
        total_int     = as.integer(round(total_resections))
      )

    if (nrow(enbloc_model_dat) >= 4) {
      fit_enbloc <- brm(
        en_bloc_int | trials(total_int) ~ technique_bin + (1 | study_id),
        data    = enbloc_model_dat,
        family  = binomial(link = "logit"),
        prior   = c(
          set_prior("normal(0, 2.5)", class = "Intercept"),
          set_prior("normal(0, 1.0)", class = "b"),
          set_prior("normal(0, 0.5)", class = "sd")
        ),
        chains  = 4, iter = 4000, warmup = 1000, seed = 123,
        backend = "cmdstanr",
        file    = file.path(tbl_dir, "brms_enbloc_esd_vs_emr")
      )

      # Summary
      draws      <- posterior::as_draws_df(fit_enbloc)
      beta_draws <- draws$b_technique_bin   # log-OR of ESD vs EMR for en bloc
      or_med <- round(exp(median(beta_draws)), 2)
      or_lo  <- round(exp(quantile(beta_draws, 0.025)), 2)
      or_hi  <- round(exp(quantile(beta_draws, 0.975)), 2)
      prob_esd_higher <- mean(beta_draws > 0)

      enbloc_summary <- tibble(
        Comparison          = "ESD vs EMR (en bloc rate)",
        OR_median           = or_med,
        OR_95CrI            = paste0(or_lo, "\u2013", or_hi),
        Prob_ESD_higher     = paste0(round(prob_esd_higher * 100, 1), "%"),
        Interpretation      = ifelse(or_med > 1,
          "ESD achieves higher en bloc rate than EMR",
          "EMR achieves higher en bloc rate than ESD")
      )
      write.csv(enbloc_summary,
                file.path(tbl_dir, "enbloc_bayesian_comparison.csv"),
                row.names = FALSE)
      message("En bloc Bayesian comparison: OR = ", or_med,
              " (", or_lo, "\u2013", or_hi, "), P(ESD>EMR) = ",
              round(prob_esd_higher * 100, 1), "%")
    }

    message("En bloc vs piecemeal plots saved.")
  } else {
    message("Not enough en bloc/piecemeal data.")
  }
} else {
  message("En bloc/piecemeal columns not found — skipping.")
}

message("Meta-regressions saved under outputs/tables/. Bubble plots saved to outputs/figures/.")

