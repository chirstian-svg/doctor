library(dplyr)
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

message("Meta-regressions saved under outputs/tables/. Bubble plots saved to outputs/figures/.")

