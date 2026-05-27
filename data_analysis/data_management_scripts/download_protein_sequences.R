library(httr)
library(readr)
library(dplyr)

# List of target species
target_species <- c(
  "Saccharomyces paradoxus", "Saccharomyces mikatae", "Saccharomyces kudriavzevii",
  "Saccharomyces arboricola", "Saccharomyces uvarum", "Saccharomyces eubayanus",
  "Saccharomyces jurei", "Saccharomyces cariocanus", "Saccharomyces bayanus",
  "Candida nivariensis", "Nakaseomyces delphensis", "Candida bracarensis",
  "Candida glabrata", "Nakaseomyces bacillosporus", "Candida castellii"
)

# Routes
protein_data <- read_tsv("./data/mt_nucleoid_PTMs_list_P-sites_20260414.tsv")
orthologs_dir <- "./data/orthologs"

if(!dir.exists(orthologs_dir)) dir.create(orthologs_dir)

get_with_retry <- function(url, attempts = 3, wait = 1) {
  for (i in seq_len(attempts)) {
    res <- GET(url)
    if (status_code(res) == 200) return(res)
    if (i < attempts) Sys.sleep(wait)
  }
  warning(paste("All", attempts, "attempts failed for URL:", url))
  res
}

# Get FASTA files
for (i in 1:nrow(protein_data)) {
  id <- protein_data$`Uniprot ID`[i]
  gene_name <- protein_data$`Standard gene name`[i]

  cat("Processing:", gene_name, "\n")
  fasta_file <- file.path(orthologs_dir, paste0(gene_name, "_orthologs.fasta"))
  
  # A: S. cerevisiae (Reference)
  ref_url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".fasta")
  ref_fasta <- get_with_retry(ref_url)

  if(status_code(ref_fasta) == 200) {
    content_text <- content(ref_fasta, "text", encoding = "UTF-8")

    # Edit header
    content_text <- sub("^>", ">Saccharomyces_cerevisiae | ", content_text)

    # Clean text
    clean_text <- trimws(content_text, which = "right")
    cat(clean_text, "\n", file = fasta_file, sep = "")
  } else {
    cat("  ! Error while getting reference", id, "\n")
    next
  }

  # B: Other species (Via UniParc)
  for (species in target_species) {
    query <- URLencode(paste0('(gene:"', gene_name, '") AND (taxonomy_name:"', species, '")'))
    search_url <- paste0("https://rest.uniprot.org/uniparc/search?query=", query, "&format=fasta&size=1")
    
    res <- get_with_retry(search_url)
    if(status_code(res) == 200) {
      content_text <- content(res, "text", encoding = "UTF-8")
      if(nchar(content_text) > 0) {
        
        #Edit header
        safe_species_name <- gsub(" ", "_", species)
        content_text <- sub("^>", paste0(">", safe_species_name, " | "), content_text)
        
        # Clean text
        clean_text <- trimws(content_text, which = "right")
        cat(clean_text, "\n", file = fasta_file, append = TRUE, sep = "")
        cat("  +", species, "- found\n")
      } else {
        cat("  -", species, "- not found\n")
      }
    }
    Sys.sleep(0.15) 
  }
  cat("Stored in:", fasta_file, "\n\n")
}
