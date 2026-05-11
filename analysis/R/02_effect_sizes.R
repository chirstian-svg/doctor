library(dplyr)
library(stringr)
library(metafor)

out_dir  <- file.path(dirname(getwd()), "outputs", "tables")
dat <- readRDS(file.path(out_dir, "clean_data.rds"))

require_cols <- function(cols) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(
      "Missing required columns for effect sizes: ",
      paste(missing, collapse = ", "),
      "\nFix: update `col_map` in R/01_load_clean.R to match your Excel headers."
    )
  }
}

compute_or <- function(events_t, n_t, events_c, n_c) {
  tmp <- tibble(
    ai = events_t,
    bi = n_t - events_t,
    ci = events_c,
    di = n_c - events_c
  )
  out <- metafor::escalc(measure = "OR", ai = ai, bi = bi, ci = ci, di = di, data = tmp)
  tibble(yi = out$yi, sei = sqrt(out$vi))
}

make_es <- function(outcome_prefix) {
  e_t <- paste0(outcome_prefix, "_clip")
  e_c <- paste0(outcome_prefix, "_no_clip")

  require_cols(c("study_id", "technique", "n_clip", "n_no_clip", e_t, e_c))

  dat |>
    transmute(
      study_id,
      study_label = paste0(
        str_to_title(str_trim(as.character(first_author_last_name))),
        " ",
        publication_year
      ),
      technique,
      n_clip,
      n_no_clip,
      events_clip = .data[[e_t]],
      events_no_clip = .data[[e_c]],
      outcome = outcome_prefix
    ) |>
    filter(
      !is.na(n_clip), !is.na(n_no_clip),
      !is.na(events_clip), !is.na(events_no_clip)
    ) |>
    rowwise() |>
    mutate(
      across(c(events_clip, events_no_clip), ~ as.numeric(.x)),
      across(c(n_clip, n_no_clip), ~ as.numeric(.x)),
      # continuity correction will be handled by escalc internally if needed
      tmp = list(compute_or(events_clip, n_clip, events_no_clip, n_no_clip)),
      yi = tmp$yi,
      sei = tmp$sei
    ) |>
    ungroup() |>
    select(-tmp)
}

es_delayed <- make_es("delayed_bleeding")
es_perforation <- make_es("perforation")
es_post_esd <- make_es("post_esd_syndrome")

es_all <- bind_rows(es_delayed, es_perforation, es_post_esd)

saveRDS(es_all, file.path(out_dir, "effect_sizes_or.rds"))
write.csv(es_all, file.path(out_dir, "effect_sizes_or.csv"), row.names = FALSE)

message("Saved effect sizes to outputs/tables/effect_sizes_or.{rds,csv}")

