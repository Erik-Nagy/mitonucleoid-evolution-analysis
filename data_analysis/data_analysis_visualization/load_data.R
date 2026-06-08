# load_data.R
# Sources this file to get: psite_long, all_sty, save_fig(), DATA_DIR, FIG_DIR, COL_* palettes
# Requires data/psite_long.tsv — run 'make process' from data_analysis/ first.

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

# ── Load P-site long table ────────────────────────────────────────────────────
psite_long_path <- file.path(DATA_DIR, "psite_long.tsv")
if (!file.exists(psite_long_path))
  stop("psite_long.tsv not found — run 'make process' from data_analysis/ first.")

psite_long <- read_tsv(psite_long_path, show_col_types = FALSE,
                       col_types = cols(.default = "c")) %>%
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
    sasa_location    = `P-site Exposure`,
    rsa_pct          = `P-site RSA (%)`,
    annotation       = Annotation,
    residue_type     = `Residue Type`
  ) %>%
  mutate(
    # Rename SASA column regardless of how Å was encoded in the TSV
    across(matches("P-site SASA"), ~ .x, .names = "sasa")
  ) %>%
  select(-matches("P-site SASA")) %>%
  mutate(
    position         = as.integer(position),
    n_psites         = as.integer(n_psites),
    exact_cons       = as.numeric(exact_cons),
    functional_sty   = as.numeric(functional_sty),
    plddt            = as.numeric(plddt),
    disorder_score   = as.numeric(disorder_score),
    sasa             = as.numeric(sasa),
    rsa_pct          = as.numeric(rsa_pct),
    structural_state = factor(structural_state,
      levels = c("Very low", "Low", "Confident", "Very high"), ordered = TRUE),
    disorder_state   = factor(disorder_state, levels = c("Ordered", "Disordered")),
    sasa_location    = factor(sasa_location,  levels = c("Buried", "Exposed")),
    residue_type     = factor(residue_type,   levels = c("pS", "pT", "pY"))
  )

message(sprintf("Loaded %d P-sites from %d proteins.",
                nrow(psite_long), n_distinct(psite_long$gene)))

# ── Load all-STY background data (produced by extract_all_sty_data.py) ───────
all_sty_path <- file.path(DATA_DIR, "all_sty_data.tsv")
if (file.exists(all_sty_path)) {
  all_sty <- read_tsv(all_sty_path, show_col_types = FALSE) %>%
    mutate(
      is_psite      = is_psite %in% c("True", "TRUE", TRUE),
      residue_type  = factor(
        ifelse(is_psite, paste0("p", residue), residue),
        levels = c("S", "pS", "T", "pT", "Y", "pY")
      ),
      rsa              = as.numeric(rsa),
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
