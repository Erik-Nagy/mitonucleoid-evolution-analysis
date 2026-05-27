library(Biostrings)
library(DECIPHER)
library(readr)
library(dplyr)
library(tidyr)

# Routes
tsv_file <- "./data/mt_nucleoid_PTMs_list_P-sites_20260414.tsv"
homologs_dir <- "./data/homologs"
output_csv <- "./data/psite_conservation_results.csv"

# Load proteins
protein_data <- read_tsv(tsv_file)

# Results for each p-site
results_list <- list()

cat("Executing MSA and analysis of P-sites conservation\n")

# Main cycle for each gene
for (i in 1:nrow(protein_data)) {
  
  gene_name <- protein_data$`Standard gene name`[i]
  p_sites_raw <- protein_data$`P-site positions`[i]
  
  p_sites <- as.numeric(trimws(unlist(strsplit(as.character(p_sites_raw), ","))))
  
  fasta_file <- file.path(homologs_dir, paste0(gene_name, "_homologs.fasta"))
  
  # Check FASTA file
  if (!file.exists(fasta_file)) next
  seqs <- readAAStringSet(fasta_file)
  if (length(seqs) < 2) next
  
  cat("Aligning and analyzing:", gene_name, "(Number of sequences:", length(seqs), ")\n")
  
  # A: Aligne sequences with library DECIPHER
  aligned_seqs <- AlignSeqs(seqs, verbose = FALSE)
  
  # B: Conversion to matrix
  aln_matrix <- as.matrix(aligned_seqs)
  
  # Get S. cerevisiae — match by name rather than assuming position 1
  ref_idx <- grep("cerevisiae", rownames(aln_matrix), ignore.case = TRUE)[1]
  if (is.na(ref_idx)) stop(paste("S. cerevisiae reference sequence not found in alignment for gene:", gene_name))
  ref_seq_aligned <- aln_matrix[ref_idx, ]

  # C: Map original positions to aligned positions
  # Find indexes in alignments, that are NOT spaces
  non_gap_indices <- which(ref_seq_aligned != "-")
  
  for (site in p_sites) {
    
    # Check length
    if (is.na(site) || site > length(non_gap_indices)) {
      cat("P-site", site, "is larger than length of sequence", gene_name, "- skipping.\n")
      next
    }
    
    aligned_pos <- non_gap_indices[site]
    column_aa <- aln_matrix[, aligned_pos]

    # Reference aminoacid (S, T, Y)
    ref_aa <- column_aa[ref_idx]

    # Homolog aminoacids
    homolog_aa <- column_aa[-ref_idx]
    num_homologs <- length(homolog_aa)
    
    # D: Calculate consevations
    # 1. Exact match (Same aminoacid as reference)
    exact_matches <- sum(homolog_aa == ref_aa)
    
    # 2. Functional match (Still S, T or Y)
    functional_matches <- sum(homolog_aa %in% c("S", "T", "Y"))
    
    # 3. Mutated
    mutated <- sum(!(homolog_aa %in% c("S", "T", "Y", "-")))
    
    # 4. Gap
    gaps <- sum(homolog_aa == "-")
    
    # Saved to temp data.frame
    res <- data.frame(
      Gene = gene_name,
      Uniprot_ID = protein_data$`Uniprot ID`[i],
      Original_Position = site,
      Aligned_Position = aligned_pos,
      Reference_AA = ref_aa,
      Num_Homologs = num_homologs,
      Exact_Match_Count = exact_matches,
      Exact_Match_Perc = round((exact_matches / num_homologs) * 100, 1),
      Functional_STY_Count = functional_matches,
      Functional_STY_Perc = round((functional_matches / num_homologs) * 100, 1),
      Gap_Perc = round((gaps / num_homologs) * 100, 1),
      stringsAsFactors = FALSE
    )
    
    results_list[[length(results_list) + 1]] <- res
  }
}

# 4. Save to CSV
final_results <- bind_rows(results_list)
write_csv(final_results, output_csv)
