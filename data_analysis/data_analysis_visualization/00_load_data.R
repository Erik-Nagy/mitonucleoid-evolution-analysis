# 00_load_data.R
# Load and preprocess P-site data.
# Sources this file to get: psite_long, save_fig(), DATA_DIR, FIG_DIR, COL_* palettes

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR <- "../data"
FIG_DIR  <- "figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Global theme ───────────────────────────────────────────────────────────────
theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "gray95"),
      legend.background = element_blank()
    )
)

# ── Helper: save figure as PNG ────────────────────────────────────────────────
save_fig <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
  invisible(p)
}

# ── Colour palettes ───────────────────────────────────────────────────────────
COL_DISORDER <- c("Ordered" = "#4393c3", "Disordered" = "#d6604d")
COL_SASA     <- c("Buried"  = "#762a83", "Exposed"    = "#1b7837")
COL_AA       <- c("pS" = "#e41a1c", "pT" = "#377eb8", "pY" = "#ff7f00")

# ── Load main TSV ──────────────────────────────────────────────────────────────
raw <- read.delim(
  file.path(DATA_DIR, "mt_nucleoid_PTMs_list_P-sites_processed.tsv"),
  sep             = "\t",
  check.names     = FALSE,
  stringsAsFactors = FALSE,
  encoding        = "UTF-8"
)

# Columns containing comma-separated per-site values
CSV_COLS <- c(
  "P-site positions",
  "P-site Exact Conservation (%)",
  "P-site Functional STY (%)",
  "P-site pLDDT Score",
  "P-site Structural State",
  "Metapredict Disorder Score",
  "Metapredict State",
  "P-site SASA (Å²)",
  "P-site 3D Location (SASA)"
)

# Detect SASA column regardless of Å encoding (may be Å, Â², A2, etc.)
sasa_col_name <- grep("P-site SASA", names(raw), value = TRUE)[1]
if (!is.na(sasa_col_name)) CSV_COLS[8] <- sasa_col_name

# ── Explode: one row per P-site ───────────────────────────────────────────────
psite_long <- raw %>%
  mutate(across(all_of(CSV_COLS), ~ strsplit(as.character(.), ",\\s*"))) %>%
  unnest(cols = all_of(CSV_COLS)) %>%
  mutate(across(all_of(CSV_COLS), trimws))

# ── Rename to clean snake_case ────────────────────────────────────────────────
# rename() expects c(new_name = "old_name"); handle SASA column separately
# because its name may differ depending on UTF-8 round-trip of the Å character.
psite_long <- psite_long %>%
  rename(
    systematic       = `Systematic gene name`,
    gene             = `Standard gene name`,
    uniprot_id       = `Uniprot ID`,
    n_psites         = `Number of P-sites`,
    position         = `P-site positions`,
    exact_cons       = `P-site Exact Conservation (%)`,
    functional_sty   = `P-site Functional STY (%)`,
    plddt            = `P-site pLDDT Score`,
    structural_state = `P-site Structural State`,
    disorder_score   = `Metapredict Disorder Score`,
    disorder_state   = `Metapredict State`,
    sasa_location    = `P-site 3D Location (SASA)`,
    annotation       = Annotation
  )

# Rename SASA column regardless of how Å was encoded
sasa_col_actual <- grep("P-site SASA", names(psite_long), value = TRUE)[1]
if (!is.na(sasa_col_actual) && sasa_col_actual != "sasa") {
  psite_long <- psite_long %>% rename(sasa = !!sasa_col_actual)
}

psite_long <- psite_long %>%
  mutate(
    position       = as.integer(position),
    n_psites       = as.integer(n_psites),
    exact_cons     = as.numeric(exact_cons),
    functional_sty = as.numeric(functional_sty),
    plddt          = as.numeric(plddt),
    disorder_score = as.numeric(disorder_score),
    sasa           = as.numeric(sasa),
    structural_state = factor(structural_state,
      levels = c("Very low", "Low", "Confident", "Very high"), ordered = TRUE),
    disorder_state = factor(disorder_state, levels = c("Ordered", "Disordered")),
    sasa_location  = factor(sasa_location,  levels = c("Buried", "Exposed"))
  )

# ── Add residue type (S/T/Y) from conservation CSV ────────────────────────────
psite_long <- psite_long %>%
  left_join(
    read.csv(file.path(DATA_DIR, "psite_conservation_results.csv"),
             stringsAsFactors = FALSE) %>%
      transmute(
        gene         = Gene,
        position     = Original_Position,
        residue_type = factor(
          case_when(trimws(Reference_AA) == "S" ~ "pS",
                    trimws(Reference_AA) == "T" ~ "pT",
                    trimws(Reference_AA) == "Y" ~ "pY",
                    TRUE ~ NA_character_),
          levels = c("pS", "pT", "pY")
        )
      ),
    by = c("gene", "position")
  )

message(sprintf("Loaded %d P-sites from %d proteins.",
                nrow(psite_long), n_distinct(psite_long$gene)))

# ── Load all-STY background data (produced by extract_all_sty_data.py) ───────
all_sty_path <- file.path(DATA_DIR, "all_sty_data.tsv")
if (file.exists(all_sty_path)) {
  all_sty <- read.delim(all_sty_path, stringsAsFactors = FALSE) %>%
    mutate(
      is_psite      = is_psite %in% c("True", "TRUE", TRUE),
      residue_type  = factor(
        ifelse(is_psite,
               paste0("p", residue),  # pS, pT, pY
               residue),              # S, T, Y
        levels = c("S", "pS", "T", "pT", "Y", "pY")
      ),
      disorder_state   = factor(disorder_state, levels = c("Ordered", "Disordered")),
      structural_state = factor(structural_state,
        levels = c("Very low", "Low", "Confident", "Very high"), ordered = TRUE),
      sasa_location    = factor(sasa_location, levels = c("Buried", "Exposed"))
    )
  message(sprintf("Loaded all-STY background: %d residues (%d P-sites, %d background).",
                  nrow(all_sty), sum(all_sty$is_psite), sum(!all_sty$is_psite)))
} else {
  all_sty <- NULL
  message("NOTE: all_sty_data.tsv not found — run 'make process' to generate it (plots 12-14 will be skipped).")
}
