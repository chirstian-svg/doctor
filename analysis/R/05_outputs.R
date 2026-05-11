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
    # ARR: positive = clipping reduces risk
    arr            = risk_no_clip - risk_clip,
    rr             = risk_clip / risk_no_clip
  )

# ── ARR plot (study-level, with author/year labels) ───────────────────────────
p_arr <- out |>
  filter(!is.na(arr)) |>
  ggplot(aes(x = arr, y = reorder(study_label, arr),
             color = technique, shape = technique)) +
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
    title    = "Absolute Risk Reduction (ARR) \u2014 Delayed Bleeding",
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
       width = 9,
       height = max(5, nrow(out |> filter(!is.na(arr))) * 0.5 + 2),
       dpi = 200)
message("ARR plot saved.")

# ── Pooled NNT from Bayesian posterior ───────────────────────────────────────
# NNT = 1 / ARR_pooled
# ARR_pooled = p0 - p_clip
# p_clip = (OR * p0) / (1 - p0 + OR * p0)   [standard OR-to-risk conversion]
# p0 = baseline risk = median control-arm risk across studies

compute_nnt <- function(fit, p0) {
  draws     <- as_draws_df(fit)$b_Intercept
  or_draws  <- exp(draws)
  p_clip    <- (or_draws * p0) / (1 - p0 + or_draws * p0)
  arr_draws <- p0 - p_clip
  nnt_draws <- 1 / arr_draws
  list(
    p0         = p0,
    OR_med     = round(median(or_draws), 2),
    OR_lo      = round(quantile(or_draws, 0.025), 2),
    OR_hi      = round(quantile(or_draws, 0.975), 2),
    ARR_med    = round(median(arr_draws), 4),
    ARR_lo     = round(quantile(arr_draws, 0.025), 4),
    ARR_hi     = round(quantile(arr_draws, 0.975), 4),
    NNT_med    = round(median(nnt_draws), 0),
    NNT_lo     = round(quantile(nnt_draws, 0.025), 0),
    NNT_hi     = round(quantile(nnt_draws, 0.975), 0)
  )
}

# Groups to compute NNT for
groups <- list(
  list(name = "Overall",      file = "brms_delayed_bleeding_overall.rds",  tech = NULL),
  list(name = "ESD subgroup", file = "brms_delayed_bleeding_ESD.rds",      tech = "ESD"),
  list(name = "EMR subgroup", file = "brms_delayed_bleeding_EMR.rds",      tech = "EMR")
)

nnt_rows <- list()

for (g in groups) {
  mfile <- file.path(tbl_dir, g$file)
  if (!file.exists(mfile)) {
    message("Model file not found, skipping: ", g$file)
    next
  }
  fit <- readRDS(mfile)

  # Baseline risk: median control-arm risk for this group
  if (is.null(g$tech)) {
    p0 <- median(out$risk_no_clip, na.rm = TRUE)
  } else {
    p0 <- out |> filter(technique == g$tech) |>
      summarise(m = median(risk_no_clip, na.rm = TRUE)) |> pull(m)
    if (is.na(p0) || length(p0) == 0) p0 <- median(out$risk_no_clip, na.rm = TRUE)
  }

  res <- compute_nnt(fit, p0)

  nnt_rows[[g$name]] <- tibble(
    Group            = g$name,
    Baseline_Risk    = round(res$p0, 4),
    Pooled_OR        = paste0(res$OR_med,  " (", res$OR_lo,  "\u2013", res$OR_hi,  ")"),
    ARR              = paste0(res$ARR_med, " (", res$ARR_lo, "\u2013", res$ARR_hi, ")"),
    NNT              = res$NNT_med,
    NNT_95CrI        = paste0(res$NNT_lo, "\u2013", res$NNT_hi),
    NNT_med_raw      = res$NNT_med,
    NNT_lo_raw       = res$NNT_lo,
    NNT_hi_raw       = res$NNT_hi
  )
}

if (length(nnt_rows) > 0) {
  nnt_table <- bind_rows(nnt_rows)

  # Save table
  write.csv(
    nnt_table |> select(-NNT_med_raw, -NNT_lo_raw, -NNT_hi_raw),
    file.path(tbl_dir, "pooled_NNT_summary.csv"),
    row.names = FALSE
  )

  # ── NNT summary plot ────────────────────────────────────────────────────────
  p_nnt <- nnt_table |>
    mutate(Group = factor(Group, levels = rev(Group))) |>
    ggplot(aes(y = Group, x = NNT_med_raw)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(
      aes(xmin = NNT_lo_raw, xmax = NNT_hi_raw),
      width = 0.25, linewidth = 0.9, color = "#4DBBD5"
    ) +
    geom_point(size = 5, color = "#E64B35", shape = 18) +
    geom_text(
      aes(label = paste0("NNT = ", NNT_med_raw, "\n(", NNT_lo_raw, "\u2013", NNT_hi_raw, ")")),
      hjust = -0.15, size = 3.2, color = "grey20"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
    labs(
      title    = "Pooled Number Needed to Treat (NNT) \u2014 Delayed Bleeding",
      subtitle = paste0(
        "NNT = 1 / ARR  |  ARR derived from Bayesian posterior pooled OR\n",
        "Baseline risk = median control-arm (No Clip) event rate per group"
      ),
      x        = "NNT  (positive = benefit from Clip)",
      y        = NULL,
      caption  = paste0(
        "NNT: number of patients who need prophylactic clipping to prevent 1 delayed bleeding event.\n",
        "Point = posterior median NNT; bars = 95% credible interval.\n",
        "Negative NNT (NNH) would indicate clipping increases risk."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      plot.caption  = element_text(size = 8, color = "grey50", hjust = 0),
      axis.text.y   = element_text(face = "bold", size = 11)
    )

  ggsave(file.path(fig_dir, "nnt_pooled_delayed_bleeding.png"), p_nnt,
         width = 10, height = max(4, nrow(nnt_table) * 1.2 + 2), dpi = 200)

  message("NNT plot saved.")
  message("NNT summary:")
  print(nnt_table |> select(Group, Baseline_Risk, Pooled_OR, NNT, NNT_95CrI))
}

# Save study-level table
write.csv(out, file.path(tbl_dir, "study_level_ARR_delayed_bleeding.csv"), row.names = FALSE)
message("Study-level ARR table saved.")
message("Done.")
