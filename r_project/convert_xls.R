#install.packages("readxl")
library(readxl)

table_xls <- read_excel("../data/mt_nucleoid_PTMs_list_P-sites_20260414.xlsx")[1:37, ]

# 1. Extract all unique positions across the whole column first
all_raw_pos <- unlist(strsplit(as.character(table_xls$`P-site positions`), ","))
unique_pos  <- sort(as.numeric(unique(trimws(all_raw_pos))))

# 2. Create a brand new matrix of FALSE values for all unique positions
# This is much faster than adding columns one-by-one
pos_matrix <- matrix(FALSE, nrow = nrow(table_xls), ncol = length(unique_pos))
colnames(pos_matrix) <- as.character(unique_pos)
pos_df <- as.data.frame(pos_matrix)

# 3. Fill the TRUEs
#for (i in 1:nrow(table_xls)) {
#  row_pos <- trimws(unlist(strsplit(as.character(table_xls$`P-site positions`[i]), ",")))
#  pos_df[i, row_pos] <- TRUE
#}

# 4. Combine with original metadata
final_table <- cbind(table_xls, pos_df)

# Define the output filename
output_file <- "../data/mt_nucleoid_processed_20260414.tsv"

# Export to TSV
write.table(table_xls, 
            file = output_file, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

