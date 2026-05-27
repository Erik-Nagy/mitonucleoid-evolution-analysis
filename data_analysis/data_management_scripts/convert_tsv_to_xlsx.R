library(readr)
library(writexl)

tsv_file <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.tsv"
output_file <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.xlsx"

protein_data <- read_tsv(tsv_file, show_col_types = FALSE)

write_xlsx(protein_data, output_file)
