library(readr)

p_sites_table <- read_tsv("./data/mt_nucleoid_PTMs_list_P-sites_20260414.tsv")

uniport_ids <- p_sites_table$`Uniprot ID`

if(!dir.exists("./data/proteins_pdb")) dir.create("./data/proteins_pdb")

for (id in uniport_ids) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", id, "-F1-model_v6.pdb")
  destfile <- paste0("./data/proteins_pdb/", id, ".pdb")
  if (!file.exists(destfile)) {
    try(download.file(url, destfile, mode = "wb"), silent = TRUE)
  }
  
}
