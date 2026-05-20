library(r3dmol)
library(readr)


p_sites_table <- read_tsv("../data/mt_nucleoid_processed_20260414.tsv")

uniport_ids <- p_sites_table$`Uniprot ID`

for (id in uniport_ids) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", id, "-F1-model_v6.pdb")
  destfile <- paste0("../data/proteins/", id, ".pdb")
  if (!file.exists(destfile)) {
    try(download.file(url, destfile, mode = "wb"), silent = TRUE)
  }
  
}

id <- uniport_ids[37]
cesta_k_pdb <- paste0("../data/proteins/", id, ".pdb")
pdb_text <- paste(readLines(cesta_k_pdb, warn = FALSE), collapse = "\n")

viewer <- r3dmol() %>%
  m_add_model(data = pdb_text, format = "pdb") %>%
  m_set_style(style = m_style_cartoon(color = "spectrum")) %>%
  m_add_style(
    sel = m_sel(resi = c(11,12,14,15,16,21,31,61,63,64,66,71,97,125,138,137,139,169,177,201,269,311,319,330,332,353,354,355,367,378,385,440,450,457,458,481,516,528,572,594,617,623,624,636,649), chain = "A"),
    style = c(
      m_style_stick(),
      m_style_sphere(scale = 0.3)
    )
  ) %>%
  m_zoom_to()

# Explicitné zobrazenie
viewer

