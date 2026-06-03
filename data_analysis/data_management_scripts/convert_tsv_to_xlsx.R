library(readr)
library(writexl)

args <- commandArgs(trailingOnly = TRUE)
tsv_file <- args[1]
output_file <- args[2]

protein_data <- read_tsv(tsv_file, show_col_types = FALSE, col_types = cols(.default = "c"))

write_xlsx(protein_data, output_file)
