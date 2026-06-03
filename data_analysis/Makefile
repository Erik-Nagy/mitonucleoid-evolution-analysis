# Run from data_analysis/ with Git Bash or MSYS2 (needs touch, mkdir -p, rm)
# Usage:
#   make download        fetch ortholog sequences and AlphaFold2 structures
#   make process         annotate P-sites (BLAST -> conservation -> pLDDT -> disorder -> SASA)
#   make all             run both pipelines in order
#   make clean           remove all generated data files and stamps
#   make clean-process   remove only processed output (keep downloads)

SCRIPTS := data_management_scripts
DATA    := data

INPUT_BASE     ?= mt_nucleoid_PTMs_list_P-sites_20260414
PROCESSED_BASE ?= mt_nucleoid_PTMs_list_P-sites_processed

INPUT_XLSX     := $(DATA)/$(INPUT_BASE).xlsx
INPUT_TSV      := $(DATA)/$(INPUT_BASE).tsv
COMBINED_FA    := $(DATA)/proteins_fasta/combined_proteins.fasta

ORTHOLOGS_DONE := $(DATA)/orthologs/.done
PDB_DONE       := $(DATA)/proteins_pdb/.done
HOMOLOGS_DONE  := $(DATA)/homologs/.done

CONSERVATION   := $(DATA)/psite_conservation_results.csv
PROCESSED_TSV  := $(DATA)/$(PROCESSED_BASE).tsv
PLDDT_DONE     := $(DATA)/.plddt_done
DISORDER_DONE  := $(DATA)/.disorder_done
SASA_DONE      := $(DATA)/.sasa_done
PROCESSED_XLSX := $(DATA)/$(PROCESSED_BASE).xlsx

.PHONY: all download process clean clean-process

all: download process

# ── Download pipeline ──────────────────────────────────────────────────────

download: $(ORTHOLOGS_DONE) $(PDB_DONE)

$(INPUT_TSV): $(INPUT_XLSX)
	Rscript $(SCRIPTS)/convert_xls.R $(INPUT_XLSX) $(INPUT_TSV)

$(ORTHOLOGS_DONE): $(INPUT_TSV)
	mkdir -p $(DATA)/orthologs
	Rscript $(SCRIPTS)/download_protein_sequences.R $(INPUT_TSV) $(DATA)/orthologs
	touch $@

$(PDB_DONE): $(INPUT_TSV)
	mkdir -p $(DATA)/proteins_pdb
	Rscript $(SCRIPTS)/download_pdb_alphafold2.R $(INPUT_TSV) $(DATA)/proteins_pdb
	touch $@

# ── Processing pipeline ────────────────────────────────────────────────────

process: $(PROCESSED_XLSX)

$(HOMOLOGS_DONE): $(ORTHOLOGS_DONE) $(COMBINED_FA)
	mkdir -p $(DATA)/homologs
	Rscript $(SCRIPTS)/run_blast_homologs.R $(COMBINED_FA) $(DATA)/orthologs $(DATA)/homologs
	touch $@

$(CONSERVATION): $(INPUT_TSV) $(HOMOLOGS_DONE)
	Rscript $(SCRIPTS)/find_conservations.R $(INPUT_TSV) $(DATA)/homologs $(CONSERVATION)

$(PROCESSED_TSV): $(INPUT_TSV) $(CONSERVATION)
	Rscript $(SCRIPTS)/merge_conservation_data.R $(INPUT_TSV) $(CONSERVATION) $(PROCESSED_TSV)

$(PLDDT_DONE): $(PROCESSED_TSV) $(PDB_DONE)
	Rscript $(SCRIPTS)/add_pLDDT_data.R $(PROCESSED_TSV) $(DATA)/proteins_pdb
	touch $@

$(DISORDER_DONE): $(PLDDT_DONE)
	python $(SCRIPTS)/add_ordered_disordered_data.py $(PROCESSED_TSV)
	touch $@

$(SASA_DONE): $(DISORDER_DONE) $(PDB_DONE)
	Rscript $(SCRIPTS)/add_exposed_buried_dssp.R $(PROCESSED_TSV) $(DATA)/proteins_pdb
	touch $@

$(PROCESSED_XLSX): $(SASA_DONE)
	Rscript $(SCRIPTS)/convert_tsv_to_xlsx.R $(PROCESSED_TSV) $(PROCESSED_XLSX)

# ── Clean ──────────────────────────────────────────────────────────────────

clean: clean-process
	rm -f $(INPUT_TSV)
	rm -rf $(DATA)/orthologs $(DATA)/proteins_pdb $(DATA)/homologs
	rm -f yeast_db.*

clean-process:
	rm -f $(CONSERVATION) $(PROCESSED_TSV) $(PROCESSED_XLSX)
	rm -f $(PLDDT_DONE) $(DISORDER_DONE) $(SASA_DONE)
