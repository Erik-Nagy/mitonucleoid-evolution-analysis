library(readr)
library(dplyr)
library(bio3d)

# Routes (Do after merging conservation data)
tsv_file <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.tsv" 
pdb_dir <- "./data/proteins_pdb"
output_tsv <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.tsv"

# Load table
protein_data <- read_tsv(tsv_file, show_col_types = FALSE)

# Vectors for new columns
plddt_scores <- character(nrow(protein_data))
structural_state <- character(nrow(protein_data))

# Main cycle for each protein
for (i in 1:nrow(protein_data)) {
  uniprot_id <- protein_data$`Uniprot ID`[i]
  p_sites_raw <- protein_data$`P-site positions`[i]
  gene_name <- protein_data$`Standard gene name`[i]
  
  if (is.na(p_sites_raw) || is.na(uniprot_id)) next
  
  p_sites <- as.numeric(trimws(unlist(strsplit(as.character(p_sites_raw), ","))))
  p_sites <- p_sites[!is.na(p_sites)]
  
  # Hľadanie PDB súboru podľa Uniprot ID
  pdb_files <- list.files(path = pdb_dir, pattern = uniprot_id, full.names = TRUE)
  pdb_files <- pdb_files[grep("\\.pdb$", pdb_files, ignore.case = TRUE)] # Iba .pdb súbory
  
  if (length(pdb_files) == 0) {
    cat(" PDB not found for:", uniprot_id, "(", gene_name, ")\n")
    plddt_scores[i] <- NA
    structural_state[i] <- NA
    next
  }

  pdb <- read.pdb(pdb_files[1], verbose = FALSE)

  scores <- numeric(length(p_sites))
  states <- character(length(p_sites))

  # Cycle for each p-site
  for (j in seq_along(p_sites)) {
    site <- p_sites[j]

    # Find CA for each position
    atom_idx <- which(pdb$atom$resno == site & pdb$atom$elety == "CA")

    if (length(atom_idx) > 0) {
      val <- pdb$atom$b[atom_idx[1]]
      scores[j] <- round(val, 1)

      if (val > 90) {
        states[j] <- "Very high"
      } else if (val > 70) {
        states[j] <- "Confident"
      } else if (val > 50) {
        states[j] <- "Low"
      } else {
        states[j] <- "Very low"
      }
    } else {
      scores[j] <- NA
      states[j] <- "Unknown"
    }
  }

  plddt_scores[i] <- paste(scores, collapse = ",")
  structural_state[i] <- paste(states, collapse = ",")
}

protein_data$`P-site pLDDT Score` <- plddt_scores
protein_data$`P-site Structural State` <- structural_state

if ("Annotation" %in% names(protein_data)) {
  protein_data <- protein_data %>% relocate(Annotation, .after = last_col())
}

write_tsv(protein_data, output_tsv)
