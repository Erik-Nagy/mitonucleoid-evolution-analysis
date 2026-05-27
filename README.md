# Mitonucleoid Evolution Analysis

A bioinformatics pipeline that annotates phosphorylation sites (P-sites) from 37 *Saccharomyces cerevisiae* mitochondrial nucleoid-associated proteins with evolutionary conservation, structural confidence, solvent accessibility, and intrinsic disorder predictions. The dataset covers P-sites from mass-spectrometry experiments (2026-04-14). A companion Shiny app provides interactive 3D structure visualisation of annotated P-sites.

---

## Repository Structure

```
mitonucleoid-evolution-analysis/
├── data_anlysis/                        # pipeline scripts and data (note: folder name typo)
│   ├── data_management_scripts/         # 8 R scripts + 1 Python script
│   └── data/
│       ├── mt_nucleoid_PTMs_list_P-sites_20260414.xlsx   # source input
│       ├── mt_nucleoid_PTMs_list_P-sites_20260414.tsv    # converted input
│       ├── mt_nucleoid_PTMs_list_P-sites_processed.tsv   # fully enriched output
│       ├── psite_conservation_results.csv                # intermediate conservation table
│       ├── orthologs/                   # per-gene FASTA from UniProt/UniParc
│       ├── homologs/                    # per-gene FASTA from BLAST
│       ├── proteins_pdb/                # AlphaFold2 PDB structures
│       └── proteins_fasta/              # combined proteome FASTA + per-species FAA files
└── p_sites_explorer/                    # Shiny visualisation app
    ├── app.r
    └── data/
        ├── mt_nucleoid_processed_20260414.tsv
        └── proteins/                    # PDB files for viewer
```

---

## Prerequisites

### System tools (must be on PATH)

| Tool | Purpose | Install |
|------|---------|---------|
| R ≥ 4.2 | pipeline scripts | [r-project.org](https://www.r-project.org/) |
| Python ≥ 3.9 | disorder prediction | [python.org](https://www.python.org/) |
| NCBI BLAST+ | `makeblastdb`, `blastp` | `conda install -c bioconda blast` |
| DSSP | solvent accessibility | `conda install -c bioconda dssp` |

### R packages

```r
install.packages(c("readxl", "readr", "dplyr", "tidyr", "stringr", "httr", "bio3d", "shiny", "r3dmol", "BiocManager"))
BiocManager::install(c("Biostrings", "DECIPHER"))
```

### Python packages

```bash
pip install metapredict pandas requests
```

---

## Input Data

- **`data_anlysis/data/mt_nucleoid_PTMs_list_P-sites_20260414.xlsx`** — source table with 37 proteins (columns: Systematic gene name, Standard gene name, UniProt ID, Number of P-sites, P-site positions, Annotation)
- **`data_anlysis/data/proteins_fasta/combined_proteins.fasta`** — combined proteome FASTA for the BLAST database (not tracked in git; must be provided manually before step 4)

---

## Pipeline

All scripts in steps 1–7 and 9 use `data_anlysis/` as the working directory. Step 8 must be run from `data_anlysis/data_management_scripts/`.

| Step | Script | What it does | Output |
|------|--------|--------------|--------|
| 1 | `convert_xls.R` | Convert Excel source to TSV | `data/mt_nucleoid_PTMs_list_P-sites_20260414.tsv` |
| 2 | `download_protein_sequences.R` | Fetch *S. cerevisiae* + 15-species ortholog sequences from UniProt/UniParc | `data/orthologs/*.fasta` |
| 3 | `download_pdb_alphafold2.R` | Download AlphaFold2 v6 PDB structures from EBI | `data/proteins_pdb/*.pdb` |
| 4 | `run_blast_homologs.R` | BLAST orthologs against combined proteome; keep best hit per species | `data/homologs/*.fasta` |
| 5 | `find_conservations.R` | MSA (DECIPHER) + per-site exact-match and functional STY conservation | `data/psite_conservation_results.csv` |
| 6 | `merge_conservation_data.R` | Join conservation results into main TSV | `data/mt_nucleoid_PTMs_list_P-sites_processed.tsv` |
| 7 | `add_pLDDT_data.R` | Extract AlphaFold2 pLDDT (B-factor) at P-site Cα atoms | updates processed TSV |
| 8 | `add_ordered_disordered_data.py` | Metapredict V3 disorder scores per P-site | updates processed TSV |
| 9 | `add_exposed_buried_dssp.R` | DSSP SASA → Exposed / Buried classification | updates processed TSV (final) |

---

## Target Species

Conservation is calculated against homolog sequences from the 15 fungal species below. Species marked with † had no retrievable sequences in NCBI and are excluded from the analysis; 8 species are effectively used.

**_Saccharomyces_ clade:** *S. paradoxus*, *S. mikatae*, *S. kudriavzevii*, *S. arboricola*, *S. uvarum*, *S. eubayanus*, *S. jurei* †, *S. cariocanus* †, *S. bayanus* †

**CTG clade:** *C. bracarensis*, *C. glabrata*, *C. nivariensis* †, *Nakaseomyces delphensis* †, *Nakaseomyces bacillosporus* †, *C. castellii* †

---

## Output

The final file `data/mt_nucleoid_PTMs_list_P-sites_processed.tsv` contains the original 6 columns plus the following per-site annotations (comma-separated when a protein has multiple P-sites):

| Column | Source | Description |
|--------|--------|-------------|
| `P-site Exact Conservation (%)` | MSA | % of homologs with the identical amino acid |
| `P-site Functional STY (%)` | MSA | % of homologs retaining S, T, or Y at that position |
| `P-site pLDDT Score` | AlphaFold2 | Per-residue model confidence (0–100) |
| `P-site Structural State` | AlphaFold2 | Very high (>90) / Confident (70–90) / Low (50–70) / Very low (<50) |
| `Metapredict Disorder Score` | Metapredict V3 | Per-residue disorder probability (0–1) |
| `Metapredict State` | Metapredict V3 | Ordered (<0.5) / Disordered (≥0.5) |
| `P-site SASA (Å²)` | DSSP | Raw solvent-accessible surface area at P-site Cα (Å²; comma-separated) |
| `P-site 3D Location (SASA)` | DSSP | Exposed (SASA > 20 Å²) / Buried (≤20 Å²) |

---

## P-sites Explorer (Shiny App)

Interactive viewer for inspecting P-site positions on 3D protein structures.

**Launch from the repository root:**

```r
shiny::runApp("p_sites_explorer")
```

Or open `p_sites_explorer/app.r` in RStudio and click **Run App**.

**Features:**
- Protein selector (37 UniProt IDs)
- Editable comma-separated P-site list with reset button
- Three colour schemes: Spectrum (N→C gradient), Secondary structure, Gray
- Interactive 3D molecular viewer (r3dmol) with P-sites highlighted as spheres

**Required data files** (not tracked in git):
- `p_sites_explorer/data/mt_nucleoid_processed_20260414.tsv`
- `p_sites_explorer/data/proteins/*.pdb` (one PDB per UniProt ID)

---

## Data Availability

Large binary files (FASTA, PDB, BLAST databases) are excluded from git via `.gitignore`. Steps 2 and 3 of the pipeline download them automatically from public APIs (UniProt REST API, AlphaFold EBI). The combined proteome FASTA (`proteins_fasta/combined_proteins.fasta`) must be assembled manually from the per-species FAA files before running step 4.
