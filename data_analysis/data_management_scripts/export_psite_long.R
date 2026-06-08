library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
processed_tsv    <- args[1]
conservation_csv <- args[2]
output_tsv       <- args[3]

raw <- read_tsv(processed_tsv, show_col_types = FALSE,
                col_types = cols(.default = "c"))

MULTI_VAL_COLS <- c(
  "P-site positions",
  "P-site Exact Conservation (%)",
  "P-site Functional STY (%)",
  "P-site pLDDT Score",
  "P-site Structural State",
  "Metapredict Disorder Score",
  "Metapredict State",
  "P-site 3D Location (SASA)"
)

sasa_col <- grep("P-site SASA", names(raw), value = TRUE)[1]
if (!is.na(sasa_col)) MULTI_VAL_COLS <- c(MULTI_VAL_COLS, sasa_col)

psite_long <- raw %>%
  mutate(across(all_of(MULTI_VAL_COLS), ~ strsplit(trimws(as.character(.)), ",\\s*"))) %>%
  unnest(cols = all_of(MULTI_VAL_COLS)) %>%
  mutate(across(all_of(MULTI_VAL_COLS), trimws))

# Add residue type (pS / pT / pY) from conservation results
residue_types <- read_csv(conservation_csv, show_col_types = FALSE) %>%
  transmute(
    `Standard gene name` = Gene,
    `P-site positions`   = as.character(Original_Position),
    `Residue Type`       = paste0("p", trimws(Reference_AA))
  )

psite_long <- left_join(psite_long, residue_types,
                        by = c("Standard gene name", "P-site positions"))

write_tsv(psite_long, output_tsv)
message(sprintf("Saved %d P-sites to %s", nrow(psite_long), output_tsv))
