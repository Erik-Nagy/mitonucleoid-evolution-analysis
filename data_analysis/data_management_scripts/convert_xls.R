#install.packages("readxl")
library(readxl)

table_xls <- read_excel("./data/mt_nucleoid_PTMs_list_P-sites_20260414.xlsx")[1:37, ]
output_file <- "./data/mt_nucleoid_PTMs_list_P-sites_20260414.tsv"

# Export to TSV
write.table(table_xls, 
            file = output_file, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
