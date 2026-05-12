library(dplyr)
library(ggplot2)
library(brms)
library(posterior)

tbl_dir <- file.path(dirname(getwd()), "outputs", "tables")
fig_dir <- file.path(dirname(getwd()), "outputs", "figures")
es <- readRDS(file.path(tbl_dir, "effect_sizes_rr.rds"))

# ESD+EMR study counts once in overall (not duplicated)
es <- es |> mutate(overall_inclusion = TRUE)

# ── helper: pretty outcome name ──────────────────────────────────────────────
pretty_outcome <- function(oc) {
  dplyr::case_when(
    oc == "delayed_bleeding"  ~ "Delayed Bleeding",
    oc == "perforation"       ~ "Perforation",
    oc == "post_esd_syndrome" ~ "Post-ESD Electrocoagulation Syndrome",
    TRUE ~ oc
  )
}

# ── fit or reload a brms model ───────────────────────────────────────────────
fit_meta <- function(df, model_name, weakly_informative_tau = TRUE) {
  priors <- c(set_prior("normal(0, 2.5)", class = "Intercept"))
  if (weakly_informative_tau) {
    priors <- c(priors, set_prior("normal(0, 0.5)", class = "sd"))
  } else {
    priors <- c(priors, set_prior("student_t(3, 0, 2.5)", class = "sd"))
  }
  brm(
    yi | se(sei) ~ 1 + (1 | study_id),
    data    = df,
    family  = gaussian(),
    prior   = priors,
    chains  = 4, iter = 4000, warmup = 1000, seed = 123,
    backend = "cmdstanr",
    file    = file.path(tbl_dir, paste0("brms_", model_name))
  )
}

# ── extract pooled OR (median + 95% CrI) from a fitted model ─────────────────
pooled_or <- function(fit) {
  draws <- as_draws_df(fit)$b_Intercept
  list(
    or_med = exp(median(draws)),
    or_lo  = exp(quantile(draws, 0.025)),
    or_hi  = exp(quantile(draws, 0.975))
  )
}

# ── publication-style combined forest plot ───────────────────────────────────
# Layout (top → bottom):
#   ESD subgroup header
#     ESD studies (sorted by OR)
#     ◆ Pooled ESD
#   [blank spacer]
#   Combined ESD+EMR study header  (if present)
#     ESD+EMR study row
#   [blank spacer]
#   EMR subgroup header
#     EMR studies (sorted by OR)
#     ◆ Pooled EMR
#   [blank spacer]
#   ◆ Overall pooled effect
plot_combined_forest <- function(df_all, fits, oc, file) {

  op <- pretty_outcome(oc)

  # ── build rows ──────────────────────────────────────────────────────────────
  rows <- list()
  row_idx <- 0   # manual y position (higher = top of plot)

  add_row <- function(label, or, lo, hi, technique, row_type, weight = NA) {
    row_idx <<- row_idx + 1
    tibble(
      y         = row_idx,
      label     = label,
      or        = or,
      lo        = lo,
      hi        = hi,
      technique = technique,
      row_type  = row_type,   # "study" | "pooled" | "header" | "spacer"
      weight    = weight
    )
  }

  # overall pooled — added last (bottom), so build first then prepend
  overall_pooled <- NULL
  if (!is.null(fits[["overall"]])) {
    p <- pooled_or(fits[["overall"]])
    overall_pooled <- add_row(
      label     = "Overall pooled effect",
      or = p$or_med, lo = p$or_lo, hi = p$or_hi,
      technique = "Overall", row_type = "pooled"
    )
  }

  # spacer before overall
  rows[["sp_overall"]] <- add_row("", NA, NA, NA, NA, "spacer")

  # subgroups: ESD, ESD+EMR, EMR
  subgroup_order <- c("ESD", "ESD+EMR", "EMR")
  for (tech in subgroup_order) {
    df_t <- df_all |> filter(technique == tech) |> arrange(yi)
    if (nrow(df_t) == 0) next

    # spacer between subgroups
    rows[[paste0("sp_", tech)]] <- add_row("", NA, NA, NA, NA, "spacer")

    # pooled diamond for this subgroup (only if ≥2 studies)
    tech_key <- gsub("[^A-Za-z0-9]+", "_", tech)
    if (!is.null(fits[[tech_key]]) && nrow(df_t) >= 2) {
      p <- pooled_or(fits[[tech_key]])
      rows[[paste0("pool_", tech)]] <- add_row(
        label     = paste0("Pooled ", tech),
        or = p$or_med, lo = p$or_lo, hi = p$or_hi,
        technique = tech, row_type = "pooled"
      )
    }

    # individual studies
    for (i in seq_len(nrow(df_t))) {
      s <- df_t[i, ]
      lbl <- if ("study_label" %in% names(s)) s$study_label else paste0("Study ", s$study_id)
      rows[[paste0(tech, "_", i)]] <- add_row(
        label     = lbl,
        or  = exp(s$yi),
        lo  = exp(s$yi - 1.96 * s$sei),
        hi  = exp(s$yi + 1.96 * s$sei),
        technique = tech,
        row_type  = "study",
        weight    = round(1 / s$sei^2, 1)
      )
    }

    # subgroup header (printed above studies — added after so y is higher)
    header_label <- if (tech == "ESD+EMR") "Combined ESD+EMR study" else paste0(tech, " subgroup")
    rows[[paste0("hdr_", tech)]] <- add_row(
      label = header_label, or = NA, lo = NA, hi = NA,
      technique = tech, row_type = "header"
    )
  }

  # prepend overall pooled at the very bottom (lowest y)
  if (!is.null(overall_pooled)) {
    overall_pooled$y <- row_idx + 1
    row_idx <<- row_idx + 1
    rows[["overall_pooled"]] <- overall_pooled
  }

  df_plot <- bind_rows(rows) |>
    mutate(
      y         = max(y) - y + 1,   # flip so top of list = top of plot
      is_pooled = row_type == "pooled",
      is_header = row_type == "header",
      label_bold = is_pooled | is_header
    )

  # colour palette
  pal <- c(
    "ESD"     = "#E64B35",
    "EMR"     = "#4DBBD5",
    "ESD+EMR" = "#00A087",
    "Overall" = "#3C5488"
  )

  # x range for log scale
  x_min <- 0.05; x_max <- 50

  p <- ggplot(df_plot) +
    # reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +

    # CI bars for studies
    geom_errorbar(
      data = df_plot |> filter(row_type == "study"),
      aes(y = y, xmin = lo, xmax = hi, color = technique),
      width = 0.3, linewidth = 0.6
    ) +
    # study points
    geom_point(
      data = df_plot |> filter(row_type == "study"),
      aes(y = y, x = or, color = technique, size = weight),
      shape = 15
    ) +

    # pooled diamonds
    geom_errorbar(
      data = df_plot |> filter(row_type == "pooled"),
      aes(y = y, xmin = lo, xmax = hi, color = technique),
      width = 0.5, linewidth = 1.0
    ) +
    geom_point(
      data = df_plot |> filter(row_type == "pooled"),
      aes(y = y, x = or, color = technique),
      shape = 18, size = 5
    ) +

    # y-axis labels
    scale_y_continuous(
      breaks = df_plot$y,
      labels = df_plot$label,
      expand = expansion(add = 0.8)
    ) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10")
    ) +
    scale_color_manual(
      name   = "Technique",
      values = pal,
      breaks = c("ESD", "EMR", "ESD+EMR", "Overall")
    ) +
    scale_size_continuous(range = c(2, 6), guide = "none") +

    labs(
      title    = paste0("Figure. Prophylactic Clipping vs No Clipping \u2014 ", op),
      subtitle = "Bayesian random-effects meta-analysis  |  Effect measure: Risk Ratio (log scale)",
      x        = "Risk Ratio (95% CrI)  [< 1 favours Clipping]",
      y        = NULL,
      caption  = paste0(
        "Squares and horizontal lines show study-level RR and 95% CrI.\n",
        "\u25c6 Diamonds show Bayesian random-effects pooled estimates with 95% credible intervals.\n",
        "Weakly informative Half-Normal(0, 0.5) prior on \u03c4; Normal(0, 2.5) prior on intercept."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      plot.caption    = element_text(size = 7, color = "grey50", hjust = 0),
      axis.text.y     = element_text(
        size  = 9,
        face  = ifelse(df_plot$label_bold[order(df_plot$y)], "bold", "plain"),
        color = ifelse(df_plot$is_header[order(df_plot$y)], "grey30", "black")
      ),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )

  n_rows <- nrow(df_plot |> filter(row_type != "spacer"))
  ggsave(filename = file, plot = p,
         width = 11, height = max(7, n_rows * 0.45 + 3), dpi = 200)
  message("Saved: ", file)
}

# ── main loop ────────────────────────────────────────────────────────────────
outcomes <- unique(es$outcome)

for (oc in outcomes) {
  df_oc <- es |> filter(outcome == oc, overall_inclusion)
  if (nrow(df_oc) < 2) next

  fits <- list()

  # Fit overall
  fits[["overall"]] <- fit_meta(df_oc,
    model_name = paste0(oc, "_overall"), weakly_informative_tau = TRUE)

  # Fit subgroups
  for (tech in sort(unique(df_oc$technique))) {
    df_t <- df_oc |> filter(technique == tech)
    if (nrow(df_t) < 2) next
    tech_key <- gsub("[^A-Za-z0-9]+", "_", tech)
    small_group <- nrow(df_t) <= 4
    fits[[tech_key]] <- fit_meta(df_t,
      model_name = paste0(oc, "_", tech_key),
      weakly_informative_tau = small_group)
  }

  # Combined publication-style forest plot
  plot_combined_forest(
    df_all = df_oc,
    fits   = fits,
    oc     = oc,
    file   = file.path(fig_dir, paste0("forest_", oc, "_combined.png"))
  )

  # Individual subgroup plots (separate files, with pooled diamond)
  plot_forest_simple <- function(df, fit, title, file) {
    op <- pretty_outcome(unique(df$outcome))

    # study rows
    study_rows <- df |>
      mutate(
        or    = exp(yi),
        lo    = exp(yi - 1.96 * sei),
        hi    = exp(yi + 1.96 * sei),
        label = if ("study_label" %in% names(df)) study_label else paste0("Study ", study_id),
        row_type = "study"
      ) |>
      arrange(yi)

    # pooled row from model
    pooled_row <- NULL
    if (!is.null(fit)) {
      p <- pooled_or(fit)
      tech_val <- unique(df$technique)
      pooled_row <- tibble(
        or       = p$or_med,
        lo       = p$or_lo,
        hi       = p$or_hi,
        label    = paste0("Pooled ", paste(tech_val, collapse = "/")),
        technique = tech_val[1],
        row_type = "pooled"
      )
    }

    plot_df <- bind_rows(study_rows, pooled_row) |>
      mutate(label = factor(label, levels = c(
        if (!is.null(pooled_row)) pooled_row$label,
        rev(study_rows$label)
      )))

    ggplot(plot_df, aes(y = label, x = or, color = technique)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      # CI bars
      geom_errorbar(aes(xmin = lo, xmax = hi,
                        linewidth = ifelse(row_type == "pooled", 1.0, 0.6)),
                    width = 0.25, show.legend = FALSE) +
      # study squares
      geom_point(data = ~ filter(., row_type == "study"),
                 aes(shape = technique), size = 3) +
      # pooled diamond
      geom_point(data = ~ filter(., row_type == "pooled"),
                 shape = 18, size = 6) +
      # CI text label on pooled row
      geom_text(data = ~ filter(., row_type == "pooled"),
                aes(label = paste0(round(or, 2), " (", round(lo, 2), "\u2013", round(hi, 2), ")")),
                hjust = -0.15, size = 3, color = "grey20") +
      scale_x_log10(
        breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
        labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10")
      ) +
      scale_linewidth_identity() +
      scale_color_manual(
        name   = "Technique",
        values = c("ESD" = "#E64B35", "EMR" = "#4DBBD5", "ESD+EMR" = "#00A087"),
        drop   = FALSE
      ) +
      scale_shape_manual(
        name   = "Technique",
        values = c("ESD" = 15, "EMR" = 15, "ESD+EMR" = 15),
        drop   = FALSE
      ) +
      labs(
        title    = title,
        subtitle = paste0("Outcome: ", op, "  |  Effect measure: Odds Ratio (log scale)"),
        x        = "Odds Ratio (95% CrI)  [OR < 1 favours Clip]",
        y        = "Study (First Author, Year)",
        caption  = paste0(
          "Squares = study-level OR; horizontal lines = 95% CI.\n",
          "\u25c6 Diamond = Bayesian pooled estimate with 95% credible interval.\n",
          "Dashed line at OR = 1 indicates no effect."
        )
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 10, color = "grey40"),
        plot.caption    = element_text(size = 8,  color = "grey50", hjust = 0),
        legend.position = "bottom",
        legend.title    = element_text(face = "bold"),
        axis.text.y     = element_text(
          face = ifelse(levels(plot_df$label) %in% (if (!is.null(pooled_row)) pooled_row$label else ""), "bold", "plain")
        )
      )

    ggsave(filename = file, width = 9,
           height = max(5, nrow(plot_df) * 0.55 + 2), dpi = 200)
  }

  plot_forest_simple(df_oc, fits[["overall"]],
    title = paste0("Clip vs No Clip \u2014 ", gsub("_", " ", tools::toTitleCase(oc)), " \u2014 All Studies"),
    file  = file.path(fig_dir, paste0("forest_", oc, "_overall.png")))

  for (tech in sort(unique(df_oc$technique))) {
    df_t <- df_oc |> filter(technique == tech)
    if (nrow(df_t) < 2) next
    tech_key <- gsub("[^A-Za-z0-9]+", "_", tech)
    plot_forest_simple(df_t, fits[[tech_key]],
      title = paste0("Clip vs No Clip \u2014 ", gsub("_", " ", tools::toTitleCase(oc)), " \u2014 ", tech, " Subgroup"),
      file  = file.path(fig_dir, paste0("forest_", oc, "_", tech_key, ".png")))
  }
}

message("All forest plots saved to outputs/figures/.")
