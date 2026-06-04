library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
processed_tsv    <- args[1]
conservation_csv <- args[2]
output_tsv       <- args[3]

raw <- read_tsv(processed_tsv, show_col_types = FALSE,
                check.names = FALSE, col_types = cols(.default = "c"))

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
  mutate(across(all_of(MULTI_VAL_COLS), trimws)) %>%
  rename(
    systematic       = `Systematic gene name`,
    gene             = `Standard gene name`,
    uniprot_id       = `Uniprot ID`,
    n_psites         = `Number of P-sites`,
    position         = `P-site positions`,
    exact_cons       = `P-site Exact Conservation (%)`,
    functional_sty   = `P-site Functional STY (%)`,
    plddt            = `P-site pLDDT Score`,
    structural_state = `P-site Structural State`,
    disorder_score   = `Metapredict Disorder Score`,
    disorder_state   = `Metapredict State`,
    sasa_location    = `P-site 3D Location (SASA)`,
    annotation       = Annotation
  ) %>%
  mutate(
    position       = as.integer(position),
    n_psites       = as.integer(n_psites),
    exact_cons     = as.numeric(exact_cons),
    functional_sty = as.numeric(functional_sty),
    plddt          = as.numeric(plddt),
    disorder_score = as.numeric(disorder_score)
  )

sasa_col_actual <- grep("P-site SASA", names(psite_long), value = TRUE)[1]
if (!is.na(sasa_col_actual) && sasa_col_actual != "sasa") {
  psite_long <- rename(psite_long, sasa = !!sasa_col_actual)
}
psite_long <- mutate(psite_long, sasa = as.numeric(sasa))

residue_types <- read_csv(conservation_csv, show_col_types = FALSE) %>%
  transmute(
    gene         = Gene,
    position     = Original_Position,
    residue_type = case_when(
      trimws(Reference_AA) == "S" ~ "pS",
      trimws(Reference_AA) == "T" ~ "pT",
      trimws(Reference_AA) == "Y" ~ "pY",
      TRUE ~ NA_character_
    )
  )

psite_long <- left_join(psite_long, residue_types, by = c("gene", "position"))

write_tsv(psite_long, output_tsv)
message(sprintf("Saved %d P-sites to %s", nrow(psite_long), output_tsv))
