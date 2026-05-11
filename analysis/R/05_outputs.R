library(dplyr)
library(ggplot2)
library(stringr)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
dat <- readRDS(file.path(tbl_dir, "clean_data.rds"))

# --- ITT / ARR / NNT calculations (study-level, delayed bleeding) ---
calc_risks <- function(events, n) {
  risk <- events / n
  ifelse(is.finite(risk), risk, NA_real_)
}

if (all(c("delayed_bleeding_clip", "delayed_bleeding_no_clip", "n_clip", "n_no_clip") %in% names(dat))) {
  out <- dat |>
    transmute(
      study_id,
      study_label = paste0(
        str_to_title(str_trim(as.character(first_author_last_name))),
        " ", publication_year
      ),
      technique,
      events_clip    = suppressWarnings(as.numeric(delayed_bleeding_clip)),
      n_clip         = suppressWarnings(as.numeric(n_clip)),
      events_no_clip = suppressWarnings(as.numeric(delayed_bleeding_no_clip)),
      n_no_clip      = suppressWarnings(as.numeric(n_no_clip)),
      risk_clip      = calc_risks(suppressWarnings(as.numeric(delayed_bleeding_clip)),
                                  suppressWarnings(as.numeric(n_clip))),
      risk_no_clip   = calc_risks(suppressWarnings(as.numeric(delayed_bleeding_no_clip)),
                                  suppressWarnings(as.numeric(n_no_clip))),
      # ARR: positive = clipping reduces risk
      arr            = risk_no_clip - risk_clip,
      # NNT: 1/ARR; positive = benefit, negative = harm (NNH)
      nnt            = ifelse(arr != 0 & !is.na(arr), 1 / arr, NA_real_),
      rr             = risk_clip / risk_no_clip
    )

  write.csv(out, file.path(tbl_dir, "itt_attributable_risk_nnt_delayed_bleeding.csv"), row.names = FALSE)

  # ---- ARR plot with author/year labels ----
  p_arr <- out |>
    filter(!is.na(arr)) |>
    ggplot(aes(x = arr, y = reorder(study_label, arr), color = technique, shape = technique)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(
      name   = "Technique",
      values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087"),
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "Technique",
      values = c("ESD" = 16, "EMR" = 17, "ESD+EMR" = 15),
      drop   = FALSE
    ) +
    labs(
      title    = "Absolute Risk Reduction (ARR) — Delayed Bleeding",
      subtitle = "ARR = Risk (No Clip) \u2212 Risk (Clip)  |  Positive values favour Clip",
      x        = "Absolute Risk Reduction (ARR)",
      y        = "Study (First Author, Year)",
      caption  = "Dashed line at ARR = 0 indicates no difference.\nPositive ARR = clipping reduces bleeding risk."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 10, color = "grey40"),
      plot.caption    = element_text(size = 8,  color = "grey50"),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold")
    )

  ggsave(file.path(fig_dir, "arr_delayed_bleeding.png"), p_arr,
         width = 9, height = max(5, nrow(out |> filter(!is.na(arr))) * 0.5 + 2), dpi = 200)

  # ---- NNT plot ----
  p_nnt <- out |>
    filter(!is.na(nnt), is.finite(nnt)) |>
    ggplot(aes(x = nnt, y = reorder(study_label, nnt), color = technique, shape = technique)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(
      name   = "Technique",
      values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087"),
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "Technique",
      values = c("ESD" = 16, "EMR" = 17, "ESD+EMR" = 15),
      drop   = FALSE
    ) +
    labs(
      title    = "Number Needed to Treat (NNT) \u2014 Delayed Bleeding",
      subtitle = "NNT = 1 / ARR  |  Positive = benefit from Clip; Negative = harm (NNH)",
      x        = "NNT (positive = benefit, negative = harm)",
      y        = "Study (First Author, Year)",
      caption  = "NNT: number of patients needing Clip to prevent 1 delayed bleeding event.\nNegative values (NNH) indicate clipping may increase risk in that study."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 10, color = "grey40"),
      plot.caption    = element_text(size = 8,  color = "grey50"),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold")
    )

  ggsave(file.path(fig_dir, "nnt_delayed_bleeding.png"), p_nnt,
         width = 9, height = max(5, nrow(out |> filter(!is.na(nnt), is.finite(nnt))) * 0.5 + 2), dpi = 200)

  message("ARR and NNT plots saved to outputs/figures/.")
  message("Summary table saved to outputs/tables/itt_attributable_risk_nnt_delayed_bleeding.csv")
}

message("Wrote summary tables/figures to outputs/.")
