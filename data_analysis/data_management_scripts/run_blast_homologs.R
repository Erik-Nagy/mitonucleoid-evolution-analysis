library(Biostrings)
library(dplyr)

# Routes
db_fasta <- "./data/proteins_fasta/combined_proteins.fasta"
orthologs_dir <- "./data/orthologs"
homologs_dir <- "./data/homologs"
db_name <- "yeast_db"

if(!dir.exists(homologs_dir)) dir.create(homologs_dir)

# Create BLAST database
makeblastdb_exit <- system(paste('makeblastdb -in', shQuote(db_fasta), '-dbtype prot -out', db_name))
if (makeblastdb_exit != 0) stop("makeblastdb failed — ensure NCBI BLAST+ is installed and on PATH.")

# Load database and extract species
all_proteins <- readAAStringSet(db_fasta)
full_headers <- names(all_proteins)
short_ids <- sub("\\s.*", "", full_headers)
names(all_proteins) <- short_ids 

# Extract names of species
species_list <- character(length(full_headers))

# A) NCBI format: [Saccharomyces cerevisiae]
ncbi_idx <- grepl("\\[.*?\\]", full_headers)
species_list[ncbi_idx] <- gsub("\\[|\\]", "", regmatches(full_headers[ncbi_idx], regexpr("\\[.*?\\]", full_headers[ncbi_idx])))

# B) Uniport format: OS=Saccharomyces cerevisiae
uniprot_idx <- !ncbi_idx & grepl("OS=", full_headers)
species_list[uniprot_idx] <- gsub("OS=", "", regmatches(full_headers[uniprot_idx], regexpr("OS=[a-zA-Z0-9_\\. ]+", full_headers[uniprot_idx])))

# C) Fallback
species_list[!ncbi_idx & !uniprot_idx] <- "Unknown_Species"

metadata <- data.frame(SubjectID = short_ids, Species = species_list, stringsAsFactors = FALSE)

# Get files
query_files <- list.files(orthologs_dir, pattern = "\\.(fasta|fa|faa)$", full.names = TRUE)
if(length(query_files) == 0) stop("No FASTA files in orthologs dir")

# Process sequences
for (query_fasta in query_files) {

  file_name_clean <- tools::file_path_sans_ext(basename(query_fasta))
  base_name <- sub("_orthologs", "", file_name_clean) 
  output_fasta <- file.path(homologs_dir, paste0(base_name, "_homologs.fasta"))
  results_table <- paste0("temp_results_", base_name, ".tab")

  cat("Processing gene:", base_name, "\n")

  query_seqs <- readAAStringSet(query_fasta)
  ref_seq <- query_seqs[1] 

  ref_id <- sub("\\s.*", "", names(ref_seq))
  names(ref_seq) <- paste0(ref_id, " | Saccharomyces cerevisiae (Query)")

  cmd <- paste('blastp -query', shQuote(query_fasta), '-db', db_name, '-out', shQuote(results_table), '-outfmt 6 -evalue 1e-5')
  blastp_exit <- system(cmd)
  if (blastp_exit != 0) stop(paste("blastp failed for gene:", base_name, "— ensure NCBI BLAST+ is installed and on PATH."))

  if (file.exists(results_table) && file.info(results_table)$size > 0) {
    blast_results <- read.table(results_table, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(blast_results) <- c("QueryID", "SubjectID", "PercID", "Length", "Mismatches", "GapOpens", "Qstart", "Qend", "Sstart", "Send", "Evalue", "Bitscore")

    blast_results <- left_join(blast_results, metadata, by = "SubjectID")

    best_matches <- blast_results %>% 
      arrange(Evalue) %>% 
      distinct(Species, .keep_all = TRUE) %>%
      filter(!grepl("Saccharomyces cerevisiae", Species, ignore.case = TRUE)) 

    valid_ids <- best_matches$SubjectID

    if(length(valid_ids) > 0) {
      homolog_seqs <- all_proteins[valid_ids]
      names(homolog_seqs) <- paste0(best_matches$SubjectID, " | ", best_matches$Species)

      final_seqs <- c(ref_seq, homolog_seqs)

      writeXStringSet(final_seqs, output_fasta)
      cat("Homologs count:", length(valid_ids) + 1, "\n")
      cat("Stored in:", output_fasta, "\n\n")
    } else {
      writeXStringSet(ref_seq, output_fasta)
      cat("Other species not found, stored just reference sequence.\n\n")
    }
  } else {
    writeXStringSet(ref_seq, output_fasta)
    cat("BLAST found no results, stored just reference sequence.\n\n")
  }

  if(file.exists(results_table)) file.remove(results_table)
}
