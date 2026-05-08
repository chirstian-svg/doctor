library(readxl)
library(dplyr)
library(stringr)
library(janitor)

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_nr_to_na <- function(x) {
  if (is.numeric(x)) return(x)
  x_chr <- as.character(x)
  x_chr <- str_trim(x_chr)
  x_chr[x_chr %in% c("", "NA", "N/A", "NR", "na", "n/a", "nr")] <- NA_character_
  x_chr
}

parse_numeric_loose <- function(x) {
  # Handles: "19.74mm", "30 mm", "37·2 mm", "7.8 mm", "NR"
  x_chr <- parse_nr_to_na(x)
  if (is.numeric(x_chr)) return(x_chr)
  x_chr <- str_replace_all(x_chr, "\u00b7", ".")
  x_chr <- str_replace_all(x_chr, ",", ".")
  x_chr <- str_replace_all(x_chr, "[^0-9.\\-]+", "")
  suppressWarnings(as.numeric(x_chr))
}

excel_path <- file.path(dirname(dirname(getwd())), "data.xlsx")
raw <- readxl::read_excel(excel_path) |> janitor::clean_names()

# ---- Column mapping (edit this if Excel headers change) ----
# We keep this explicit so later scripts can rely on consistent names.
col_map <- c(
  # sample sizes
  n_clip = "number_of_patients_assessed_clip_intervention",
  n_no_clip = "number_of_patients_assessed_control_no_clip",

  # outcomes (events) — control/no-clip side
  delayed_bleeding_no_clip = "no_of_delayed_bleendings_control_no_clip_12",
  perforation_no_clip = "no_of_perforation_control_no_clip_13",
  post_esd_syndrome_no_clip = "no_of_post_esd_electrocoagulation_syndrome_control_no_clip_14",

  # outcomes (events) — clip/intervention side (duplicate-named columns, suffixed _32/_33/_34)
  delayed_bleeding_clip = "no_of_delayed_bleendings_control_no_clip_32",
  perforation_clip = "no_of_perforation_control_no_clip_33",
  post_esd_syndrome_clip = "no_of_post_esd_electrocoagulation_syndrome_control_no_clip_34",

  # moderators / covariates
  antiplatelet_clip = "no_percent_treated_with_antiplate_medications_cclip_intervention",
  anticoag_clip = "no_of_patient_treated_with_anticuagulation_clip_intervention",
  avg_size = "average_size_20",

  technique_label = "technique",
  technique_code = "techinque_1_esd_2_emr_3_esd_emr"
)

missing_cols <- setdiff(unname(col_map), names(raw))
if (length(missing_cols) > 0) {
  message("Missing expected columns in Excel (update col_map):")
  message(paste0("- ", missing_cols, collapse = "\n"))
}

dat <- raw

rename_existing <- function(df, new, old) {
  if (old %in% names(df)) dplyr::rename(df, !!new := !!sym(old)) else df
}
for (nm in names(col_map)) {
  dat <- rename_existing(dat, nm, col_map[[nm]])
}

dat <- dat |>
  mutate(
    across(everything(), parse_nr_to_na),
    across(
      c(
        n_clip, n_no_clip,
        delayed_bleeding_clip, delayed_bleeding_no_clip,
        perforation_clip, perforation_no_clip,
        post_esd_syndrome_clip, post_esd_syndrome_no_clip,
        antiplatelet_clip, anticoag_clip
      ),
      parse_numeric_loose
    ),
    avg_size = parse_numeric_loose(avg_size),
    technique_code = parse_numeric_loose(technique_code),
    technique = case_when(
      technique_code == 1 ~ "ESD",
      technique_code == 2 ~ "EMR",
      technique_code == 3 ~ "ESD+EMR",
      TRUE ~ as.character(technique_label)
    ),
    technique = str_trim(technique)
  )

# Add a stable study id if not present (Excel seems to be one row per study)
dat <- dat |> mutate(study_id = row_number())

saveRDS(dat, file = "outputs/tables/clean_data.rds")
write.csv(dat, file = "outputs/tables/clean_data.csv", row.names = FALSE)

message("Saved cleaned data to outputs/tables/clean_data.{rds,csv}")

