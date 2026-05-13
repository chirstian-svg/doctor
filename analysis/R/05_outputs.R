library(dplyr)
library(ggplot2)
library(stringr)
library(brms)
library(posterior)

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

# ── Pooled NNT from Bayesian posterior ───────────────────────────────────────
# For RCT meta-analysis: use RR directly (not OR-to-risk conversion)
# p_clip = RR * p0
# ARR = p0 - p_clip = p0 * (1 - RR)
# NNT = 1 / ARR

compute_nnt_rr <- function(fit, p0) {
  draws     <- as_draws_df(fit)$b_Intercept
  rr_draws  <- exp(draws)                  # posterior RR draws
  arr_draws <- p0 * (1 - rr_draws)         # ARR = p0 * (1 - RR)
  nnt_draws <- 1 / arr_draws               # NNT = 1 / ARR
  list(
    p0        = p0,
    p0_pct    = round(p0 * 100, 2),
    RR_med    = round(median(rr_draws), 3),
    RR_lo     = round(quantile(rr_draws, 0.025), 3),
    RR_hi     = round(quantile(rr_draws, 0.975), 3),
    ARR_med   = round(median(arr_draws), 4),
    ARR_pct   = round(median(arr_draws) * 100, 2),
    NNT_med   = abs(round(median(nnt_draws), 0)),
    NNT_lo    = abs(round(quantile(nnt_draws, 0.025), 0)),
    NNT_hi    = abs(round(quantile(nnt_draws, 0.975), 0)),
    # clip arm predicted risk
    p_clip_med = round(median(rr_draws * p0) * 100, 2)
  )
}

groups <- list(
  list(name = "Overall", file = "brms_delayed_bleeding_overall.rds", tech = NULL),
  list(name = "ESD",     file = "brms_delayed_bleeding_ESD.rds",     tech = "ESD"),
  list(name = "EMR",     file = "brms_delayed_bleeding_EMR.rds",     tech = "EMR")
)

bar_rows  <- list()
nnt_rows  <- list()

for (g in groups) {
  mfile <- file.path(tbl_dir, g$file)
  if (!file.exists(mfile)) { message("Skipping (not found): ", g$file); next }
  fit <- readRDS(mfile)

  if (is.null(g$tech)) {
    p0 <- mean(out$risk_no_clip, na.rm = TRUE)   # pooled baseline = mean across all studies
  } else {
    p0 <- out |> filter(technique == g$tech) |>
      summarise(m = mean(risk_no_clip, na.rm = TRUE)) |> pull(m)
    if (is.na(p0) || length(p0) == 0) p0 <- mean(out$risk_no_clip, na.rm = TRUE)
  }

  res <- compute_nnt_rr(fit, p0)

  # rows for grouped bar chart
  bar_rows[[paste0(g$name, "_noclip")]] <- tibble(
    group = g$name, arm = "No Clipping",
    risk_pct = res$p0_pct,
    arr_pct  = res$ARR_pct,
    nnt      = res$NNT_med
  )
  bar_rows[[paste0(g$name, "_clip")]] <- tibble(
    group = g$name, arm = "Clipping",
    risk_pct = res$p_clip_med,
    arr_pct  = res$ARR_pct,
    nnt      = res$NNT_med
  )

  nnt_rows[[g$name]] <- tibble(
    Group         = g$name,
    Baseline_Risk = paste0(res$p0_pct, "%"),
    Pooled_RR     = paste0(res$RR_med, " (", res$RR_lo, "\u2013", res$RR_hi, ")"),
    ARR_pct       = paste0(res$ARR_pct, "%"),
    NNT           = res$NNT_med,
    NNT_95CrI     = paste0(res$NNT_lo, "\u2013", res$NNT_hi)
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
