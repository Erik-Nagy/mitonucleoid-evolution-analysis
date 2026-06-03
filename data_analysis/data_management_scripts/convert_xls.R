#install.packages("readxl")
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
table_xls <- read_excel(args[1])[1:37, ]
output_file <- args[2]

# Export to TSV
write.table(table_xls, 
            file = output_file, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
