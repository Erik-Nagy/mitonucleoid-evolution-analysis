library(readr)
library(dplyr)
library(bio3d)

# Max ASA (Å²) per residue — Tien et al. (2013) empirical values
MAX_ASA <- c(
  ALA = 121, ARG = 265, ASN = 187, ASP = 187, CYS = 148,
  GLN = 214, GLU = 214, GLY =  97, HIS = 216, ILE = 195,
  LEU = 191, LYS = 230, MET = 203, PHE = 228, PRO = 154,
  SER = 143, THR = 163, TRP = 264, TYR = 255, VAL = 165
)

args <- commandArgs(trailingOnly = TRUE)
tsv_file   <- args[1]
pdb_dir    <- args[2]
output_tsv <- args[1]

protein_data   <- read_tsv(tsv_file, show_col_types = FALSE)
location_state <- character(nrow(protein_data))
sasa_raw       <- character(nrow(protein_data))
rsa_out        <- character(nrow(protein_data))

mkdssp_path <- Sys.which("mkdssp")
if (nchar(mkdssp_path) == 0) stop("mkdssp not found on PATH — install DSSP (e.g. conda install -c bioconda dssp) and ensure it is accessible.")

# bio3d::dssp() strips CRYST1 via write.pdb, which mkdssp 4.x requires.
# Run mkdssp directly on the original PDB file instead.
run_dssp <- function(pdb_file, mkdssp_path) {
  outfile <- tempfile(fileext = ".dssp")
  on.exit(unlink(outfile, force = TRUE))
  stderr_out <- system2(
    mkdssp_path,
    args = c(shQuote(pdb_file), shQuote(outfile)),
    stdout = FALSE, stderr = TRUE
  )
  exit_code <- attr(stderr_out, "status")
  if (!is.null(exit_code)) {
    cat("mkdssp failed (exit", exit_code, ") for:", basename(pdb_file), "\n")
    if (length(stderr_out) > 0) cat("  stderr:", paste(stderr_out, collapse = "\n  "), "\n")
    return(NULL)
  }
  if (!file.exists(outfile)) {
    cat("mkdssp produced no output file for:", basename(pdb_file), "\n")
    if (length(stderr_out) > 0) cat("  stderr:", paste(stderr_out, collapse = "\n  "), "\n")
    return(NULL)
  }
  lines      <- readLines(outfile)
  header_idx <- which(substring(lines, 1, 3) == "  #")
  if (length(header_idx) == 0) return(NULL)
  data_lines <- lines[(header_idx[1] + 1):length(lines)]
  data_lines <- data_lines[substring(data_lines, 14, 14) != "!"]
  res_num <- as.numeric(substring(data_lines, 6, 10))
  acc     <- as.numeric(substring(data_lines, 35, 38))
  list(acc = setNames(acc, as.character(res_num)))
}

for (i in 1:nrow(protein_data)) {
  uniprot_id  <- protein_data$`Uniprot ID`[i]
  p_sites_raw <- protein_data$`P-site positions`[i]

  p_sites <- as.numeric(trimws(unlist(strsplit(as.character(p_sites_raw), ","))))
  p_sites <- p_sites[!is.na(p_sites)]

  pdb_files <- list.files(path = pdb_dir, pattern = uniprot_id, full.names = TRUE)
  pdb_files <- pdb_files[grep("\\.pdb$", pdb_files, ignore.case = TRUE)]

  if (length(pdb_files) == 0) {
    cat("No PDB for:", uniprot_id, "\n")
    location_state[i] <- NA
    sasa_raw[i]       <- NA
    rsa_out[i]        <- NA
    next
  }

  pdb       <- read.pdb(pdb_files[1], verbose = FALSE)
  locs      <- character(length(p_sites))
  sasa_vals <- character(length(p_sites))
  rsa_vals  <- character(length(p_sites))

  dssp_data <- run_dssp(pdb_files[1], mkdssp_path)

  if (is.null(dssp_data)) {
    cat("DSSP failed for:", uniprot_id, "\n")
    location_state[i] <- paste(rep("Unknown", length(p_sites)), collapse = ",")
    sasa_raw[i]       <- paste(rep("NA", length(p_sites)), collapse = ",")
    rsa_out[i]        <- paste(rep("NA", length(p_sites)), collapse = ",")
    next
  }

  cat(paste0("\nProcessing ", uniprot_id, ": Looking for positions [",
             paste(p_sites, collapse = ","), "]\n"))

  for (j in seq_along(p_sites)) {
    site     <- p_sites[j]
    atom_idx <- which(pdb$atom$resno == site & pdb$atom$elety == "CA")

    if (length(atom_idx) > 0) {
      res_index <- pdb$atom$resno[atom_idx[1]]
      res_name  <- pdb$atom$resid[atom_idx[1]]

      sasa_val <- dssp_data$acc[as.character(res_index)]
      if (is.null(sasa_val) || is.na(sasa_val) || sasa_val > 500) sasa_val <- NA

      max_asa <- MAX_ASA[res_name]
      rsa_val <- if (!is.na(sasa_val) && !is.na(max_asa) && max_asa > 0)
                   round((sasa_val / max_asa) * 100, 2) else NA

      cat(paste0("Position ", site, " (", res_name, "): SASA = ", sasa_val,
                 " Å², RSA = ", rsa_val, "%\n"))

      sasa_vals[j] <- if (!is.na(sasa_val)) as.character(round(sasa_val, 2)) else "NA"
      rsa_vals[j]  <- if (!is.na(rsa_val))  as.character(rsa_val)            else "NA"
      locs[j]      <- if (!is.na(rsa_val)) { if (rsa_val > 20) "Exposed" else "Buried" } else "Unknown"

    } else {
      cat(paste0("Position ", site, " not found in PDB file!\n"))
      sasa_vals[j] <- "NA"
      rsa_vals[j]  <- "NA"
      locs[j]      <- "Unknown"
    }
  }

  location_state[i] <- paste(locs,      collapse = ",")
  sasa_raw[i]       <- paste(sasa_vals, collapse = ",")
  rsa_out[i]        <- paste(rsa_vals,  collapse = ",")
}

protein_data$`P-site SASA (Å²)`          <- sasa_raw
protein_data$`P-site RSA (%)`            <- rsa_out
protein_data$`P-site 3D Location (SASA)` <- location_state

if ("Annotation" %in% names(protein_data)) {
  protein_data <- protein_data %>% relocate(Annotation, .after = last_col())
}

write_tsv(protein_data, output_tsv)
