library(readr)

args <- commandArgs(trailingOnly = TRUE)
p_sites_table <- read_tsv(args[1])
pdb_dir <- args[2]

uniport_ids <- p_sites_table$`Uniprot ID`

if(!dir.exists(pdb_dir)) dir.create(pdb_dir)

for (id in uniport_ids) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", id, "-F1-model_v6.pdb")
  destfile <- file.path(pdb_dir, paste0(id, ".pdb"))
  if (!file.exists(destfile)) {
    try(download.file(url, destfile, mode = "wb"), silent = TRUE)
  }
  
}
