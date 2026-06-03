cran_pkgs <- c("readxl", "readr", "dplyr", "tidyr", "stringr",
               "httr", "bio3d", "writexl", "shiny", "r3dmol")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("Biostrings", "DECIPHER"))
