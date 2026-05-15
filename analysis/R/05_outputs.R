library(dplyr)
library(ggplot2)
library(stringr)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
dat <- readRDS(file.path(tbl_dir, "clean_data.rds"))

calc_risks <- function(events, n) {
  risk <- events / n
  ifelse(is.finite(risk), risk, NA_real_)
}

if (!all(c("delayed_bleeding_clip", "delayed_bleeding_no_clip", "n_clip", "n_no_clip") %in% names(dat))) {
  stop("Required columns missing from clean data.")
}

# ── Study-level risk table ────────────────────────────────────────────────────
out <- dat |>
  transmute(
    study_id,
    study_label    = paste0(
      str_to_title(str_trim(as.character(first_author_last_name))),
      " ", publication_year
    ),
    technique,
    events_clip    = suppressWarnings(as.numeric(delayed_bleeding_clip)),
    n_clip         = suppressWarnings(as.numeric(n_clip)),
    events_no_clip = suppressWarnings(as.numeric(delayed_bleeding_no_clip)),
    n_no_clip      = suppressWarnings(as.numeric(n_no_clip)),
    risk_clip      = calc_risks(
                       suppressWarnings(as.numeric(delayed_bleeding_clip)),
                       suppressWarnings(as.numeric(n_clip))),
    risk_no_clip   = calc_risks(
                       suppressWarnings(as.numeric(delayed_bleeding_no_clip)),
                       suppressWarnings(as.numeric(n_no_clip))),
    arr = risk_no_clip - risk_clip,
    rr  = risk_clip / risk_no_clip
  )

groups <- list(
  list(name = "Overall", tech = NULL),
  list(name = "ESD",     tech = "ESD"),
  list(name = "EMR",     tech = "EMR")
)

bar_rows  <- list()
nnt_rows  <- list()

for (g in groups) {
  # Filter to relevant studies
  if (is.null(g$tech)) {
    sub <- out |> filter(!is.na(risk_no_clip), !is.na(risk_clip))
  } else {
    sub <- out |> filter(technique == g$tech, !is.na(risk_no_clip), !is.na(risk_clip))
  }
  if (nrow(sub) == 0) next

  # Pooled risks: sum events / sum patients (same method client used)
  total_events_noclip <- sum(sub$events_no_clip, na.rm = TRUE)
  total_n_noclip      <- sum(sub$n_no_clip,      na.rm = TRUE)
  total_events_clip   <- sum(sub$events_clip,    na.rm = TRUE)
  total_n_clip        <- sum(sub$n_clip,          na.rm = TRUE)

  p_noclip <- total_events_noclip / total_n_noclip
  p_clip   <- total_events_clip   / total_n_clip
  arr      <- p_noclip - p_clip
  nnt      <- round(1 / arr, 0)

  bar_rows[[paste0(g$name, "_noclip")]] <- tibble(
    group    = g$name,
    arm      = "No Clipping",
    risk_pct = round(p_noclip * 100, 2),
    arr_pct  = round(arr * 100, 2),
    nnt      = nnt
  )
  bar_rows[[paste0(g$name, "_clip")]] <- tibble(
    group    = g$name,
    arm      = "Clipping",
    risk_pct = round(p_clip * 100, 2),
    arr_pct  = round(arr * 100, 2),
    nnt      = nnt
  )

  nnt_rows[[g$name]] <- tibble(
    Group            = g$name,
    Events_NoClip    = total_events_noclip,
    N_NoClip         = total_n_noclip,
    Risk_NoClip_pct  = round(p_noclip * 100, 2),
    Events_Clip      = total_events_clip,
    N_Clip           = total_n_clip,
    Risk_Clip_pct    = round(p_clip * 100, 2),
    ARR_pct          = round(arr * 100, 2),
    NNT              = nnt
  )
}

# ── Grouped bar chart: ARR of Prophylactic Clipping (matches client image) ───
if (length(bar_rows) > 0) {
  bar_df <- bind_rows(bar_rows) |>
    mutate(
      group = factor(group, levels = c("Overall", "ESD", "EMR")),
      arm   = factor(arm,   levels = c("No Clipping", "Clipping"))
    )

  # annotation: one label per group (above the No Clipping bar)
  annot_df <- bar_df |>
    filter(arm == "No Clipping") |>
    mutate(label = paste0("ARR ", sprintf("%.2f", arr_pct), "%\nNNT = ", nnt))

  p_arr_bar <- ggplot(bar_df, aes(x = group, y = risk_pct, fill = arm)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.55) +
    geom_text(
      data = annot_df,
      aes(x = group, y = risk_pct + 0.3, label = label),
      inherit.aes = FALSE,
      size = 3.2, vjust = 0, fontface = "bold", color = "grey20"
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("No Clipping" = "#4472C4", "Clipping" = "#ED7D31")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.25)),
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title    = "Absolute Risk Reduction of Prophylactic Clipping",
      subtitle = "Delayed Bleeding (%) by resection technique",
      x        = "Resection Technique",
      y        = "Delayed Bleeding (%)",
      caption  = paste0(
        "Bar heights = pooled event rate derived from Bayesian posterior RR.\n",
        "ARR = Absolute Risk Reduction; NNT = Number Needed to Treat.\n",
        "Baseline risk = mean control-arm (No Clipping) event rate per group."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(size = 10, color = "grey40"),
      plot.caption    = element_text(size = 8,  color = "grey50", hjust = 0),
      legend.position = "top",
      legend.text     = element_text(size = 11),
      axis.text.x     = element_text(size = 11, face = "bold"),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(fig_dir, "arr_nnt_bar_delayed_bleeding.png"), p_arr_bar,
         width = 7, height = 6, dpi = 200)
  message("ARR/NNT bar chart saved.")
}

# ── Save NNT summary table ────────────────────────────────────────────────────
if (length(nnt_rows) > 0) {
  nnt_table <- bind_rows(nnt_rows)
  write.csv(nnt_table, file.path(tbl_dir, "pooled_NNT_summary.csv"), row.names = FALSE)
  message("NNT summary table saved.")
  print(nnt_table)
}

# ── Save study-level ARR table ────────────────────────────────────────────────
write.csv(out, file.path(tbl_dir, "study_level_ARR_delayed_bleeding.csv"), row.names = FALSE)
message("Study-level ARR table saved.")
message("Done.")
