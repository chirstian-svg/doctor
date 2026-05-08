library(dplyr)
library(ggplot2)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
dat <- readRDS(file.path(tbl_dir, "clean_data.rds"))

# --- ITT / attributable risk style calculations (study-level) ---
# Using delayed bleeding as default; extend similarly for other outcomes if needed.
calc_risks <- function(events, n) {
  risk <- events / n
  ifelse(is.finite(risk), risk, NA_real_)
}

if (all(c("delayed_bleeding_clip", "delayed_bleeding_no_clip", "n_clip", "n_no_clip") %in% names(dat))) {
  out <- dat |>
    transmute(
      study_id,
      technique,
      events_clip = delayed_bleeding_clip,
      n_clip,
      events_no_clip = delayed_bleeding_no_clip,
      n_no_clip,
      risk_clip = calc_risks(events_clip, n_clip),
      risk_no_clip = calc_risks(events_no_clip, n_no_clip),
      attributable_risk = risk_clip - risk_no_clip,
      rr = risk_clip / risk_no_clip
    )

  write.csv(out, file.path(tbl_dir, "itt_attributable_risk_delayed_bleeding.csv"), row.names = FALSE)

  p <- out |>
    filter(!is.na(attributable_risk)) |>
    ggplot(aes(x = attributable_risk, y = reorder(paste0("Study ", study_id), attributable_risk), color = technique)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2) +
    labs(
      title = "Attributable risk (Clip − No Clip) — Delayed bleeding",
      x = "Risk difference",
      y = NULL
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(fig_dir, "attributable_risk_delayed_bleeding.png"), p, width = 8, height = 5, dpi = 200)
}

# NOTE: “Non-permitted-to-treat” was requested in the call notes but is not a standard label.
# If the client clarifies the exact rule (e.g., exclude protocol violations, crossovers, or rescue clipping),
# implement it here as an alternative denominator/numerator and output a second table/figure.

message("Wrote summary tables/figures to outputs/.")

