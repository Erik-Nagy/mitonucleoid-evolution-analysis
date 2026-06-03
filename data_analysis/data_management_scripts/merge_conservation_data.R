library(readr)
library(dplyr)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
tsv_file <- args[1]
results_csv <- args[2]
output_tsv <- args[3]

# Load tables
protein_data <- read_tsv(tsv_file)
conservation_data <- read_csv(results_csv)

# Group and format
conservation_collapsed <- conservation_data %>%
  arrange(Gene, Original_Position) %>%
  group_by(Gene) %>%
  summarise(
    `P-site Exact Conservation (%)` = paste(Exact_Match_Perc, collapse = ","),
    `P-site Functional STY (%)` = paste(Functional_STY_Perc, collapse = ",")
  )

# Clean
updated_protein_data <- protein_data %>%
  left_join(conservation_collapsed, by = c("Standard gene name" = "Gene")) %>%
  
  mutate(across(where(is.character), ~ str_replace_all(., "[^[:print:]]", ""))) %>%
  mutate(across(where(is.character), ~ str_squish(.))) %>%
  
  relocate(Annotation, .after = last_col())

write_tsv(updated_protein_data, output_tsv)
