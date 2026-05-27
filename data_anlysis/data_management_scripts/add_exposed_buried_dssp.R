library(readr)
library(dplyr)
library(bio3d)

# Routes (Do after adding ordered/disordered data)
tsv_file <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.tsv" 
pdb_dir <- "./data/proteins_pdb"
output_tsv <- "./data/mt_nucleoid_PTMs_list_P-sites_processed.tsv"

protein_data <- read_tsv(tsv_file, show_col_types = FALSE)
location_state <- character(nrow(protein_data))

mkdssp_path <- Sys.which("mkdssp")
if (nchar(mkdssp_path) == 0) stop("mkdssp not found on PATH — install DSSP (e.g. conda install -c bioconda dssp) and ensure it is accessible.")

for (i in 1:nrow(protein_data)) {
  uniprot_id <- protein_data$`Uniprot ID`[i]
  p_sites_raw <- protein_data$`P-site positions`[i]

  p_sites <- as.numeric(trimws(unlist(strsplit(as.character(p_sites_raw), ","))))
  p_sites <- p_sites[!is.na(p_sites)]

  pdb_files <- list.files(path = pdb_dir, pattern = uniprot_id, full.names = TRUE)
  pdb_files <- pdb_files[grep("\\.pdb$", pdb_files, ignore.case = TRUE)]

  if (length(pdb_files) == 0) {
    cat("No PDB for:", uniprot_id, "\n")
    location_state[i] <- NA
    next
  }
  
  pdb <- read.pdb(pdb_files[1], verbose = FALSE)
  locs <- character(length(p_sites))
  
  # Start DSSP
  dssp_data <- tryCatch({
    dssp(pdb, exefile = mkdssp_path)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(dssp_data)) {
    cat("DSSP failed for:", uniprot_id, "\n")
    location_state[i] <- paste(rep("Unknown", length(p_sites)), collapse = ",")
    next
  }

  cat(paste0("\nProcessing ", uniprot_id, ": Looking for positions [", paste(p_sites, collapse=","), "]\n"))

  for (j in seq_along(p_sites)) {
    site <- p_sites[j]

    # Find CA for position
    atom_idx <- which(pdb$atom$resno == site & pdb$atom$elety == "CA")
    
    if (length(atom_idx) > 0) {
      res_index_in_pdb <- pdb$atom$resno[atom_idx[1]] 
      
      # Get SASA value from DSSP object
      sasa_val <- dssp_data$acc[as.character(res_index_in_pdb)]
      if (is.null(sasa_val) || is.na(sasa_val)) {
        sasa_val <- dssp_data$acc[res_index_in_pdb]
      }

      cat(paste0("Position ", site, ": found in PDB. SASA value from DSSP = ", sasa_val, "\n"))
      
      if (!is.na(sasa_val)) {
        if (sasa_val > 20) {
          locs[j] <- "Exposed"
        } else {
          locs[j] <- "Buried"
        }
      } else {
        locs[j] <- "Unknown"
      }
      
    } else {
      cat(paste0("Position ", site, " not found in PDB file!\n"))
      locs[j] <- "Unknown"
    }
  }
  
  location_state[i] <- paste(locs, collapse = ",")
}

protein_data$`P-site 3D Location (SASA)` <- location_state

if ("Annotation" %in% names(protein_data)) {
  protein_data <- protein_data %>% relocate(Annotation, .after = last_col())
}

write_tsv(protein_data, output_tsv)
