"""
extract_all_sty_data.py
For each of the 37 mitochondrial nucleoid proteins, extract disorder scores
(Metapredict V3), pLDDT scores (AlphaFold2 B-factors), and SASA (mkdssp)
for every S, T, and Y residue in the full sequence, then mark which ones
are known P-sites.

Usage:
    python extract_all_sty_data.py <processed_tsv> <orthologs_dir> <pdb_dir> <output_csv>

Requires: metapredict, pandas
Requires on PATH: mkdssp (DSSP 4.x — same as the existing pipeline)
"""

import sys
import os
import shutil
import subprocess
import tempfile
import pandas as pd
import metapredict as meta


def parse_cerevisiae_sequence(orthologs_dir, gene_name):
    """Return the S. cerevisiae sequence from orthologs FASTA (first entry)."""
    fasta_path = os.path.join(orthologs_dir, f"{gene_name}_orthologs.fasta")
    if not os.path.exists(fasta_path):
        return None
    seq_lines = []
    found_header = False
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if found_header:
                    break  # second header → done with first sequence
                found_header = True
            elif found_header:
                seq_lines.append(line.strip())
    return "".join(seq_lines) if seq_lines else None


def parse_pdb_plddt(pdb_dir, uniprot_id):
    """Return {residue_number: pLDDT} from AlphaFold2 PDB (B-factor on Cα atoms)."""
    pdb_path = os.path.join(pdb_dir, f"{uniprot_id}.pdb")
    if not os.path.exists(pdb_path):
        return {}
    plddt_by_res = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resnum = int(line[22:26])
                bfactor = float(line[60:66])
                plddt_by_res[resnum] = round(bfactor, 1)
    return plddt_by_res


def run_dssp(pdb_path, mkdssp_path):
    """
    Run mkdssp on a PDB file and return {residue_number: SASA_Å²}.
    Returns an empty dict if mkdssp fails.
    DSSP ACC column (positions 35-38 in data lines) is the absolute SASA in Å².
    """
    with tempfile.NamedTemporaryFile(suffix=".dssp", delete=False) as tmp:
        out_path = tmp.name
    try:
        result = subprocess.run(
            [mkdssp_path, pdb_path, out_path],
            capture_output=True, text=True
        )
        if result.returncode != 0 or not os.path.exists(out_path):
            return {}
        sasa_map = {}
        with open(out_path) as fh:
            lines = fh.readlines()
        header_idx = next(
            (i for i, l in enumerate(lines) if l.startswith("  #")), None
        )
        if header_idx is None:
            return {}
        for line in lines[header_idx + 1:]:
            if line[13:14] == "!":  # chain break marker
                continue
            try:
                resnum = int(line[5:10])
                acc = float(line[34:38])
                if acc <= 500:  # sanity cap (same as R script)
                    sasa_map[resnum] = round(acc, 2)
            except (ValueError, IndexError):
                continue
        return sasa_map
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)


def structural_state(plddt):
    if plddt is None:
        return None
    if plddt < 50:
        return "Very low"
    if plddt < 70:
        return "Low"
    if plddt < 90:
        return "Confident"
    return "Very high"


def sasa_location(sasa):
    if sasa is None:
        return None
    return "Exposed" if sasa > 20 else "Buried"


def main():
    if len(sys.argv) != 5:
        sys.exit(
            "Usage: extract_all_sty_data.py "
            "<processed_tsv> <orthologs_dir> <pdb_dir> <output_csv>"
        )

    tsv_path, orthologs_dir, pdb_dir, out_csv = sys.argv[1:]

    mkdssp_path = shutil.which("mkdssp")
    if mkdssp_path is None:
        print("WARNING: mkdssp not found on PATH — SASA columns will be empty.")
        print("         Install DSSP (e.g. conda install -c bioconda dssp) to enable SASA.")

    df = pd.read_csv(tsv_path, sep="\t")

    # Build set of known P-site (gene, position) pairs
    psite_set = set()
    for _, row in df.iterrows():
        gene = str(row["Standard gene name"]).strip()
        for p in str(row["P-site positions"]).split(","):
            try:
                psite_set.add((gene, int(float(p.strip()))))
            except ValueError:
                pass

    records = []
    for _, row in df.iterrows():
        gene = str(row["Standard gene name"]).strip()
        uniprot_id = str(row["Uniprot ID"]).strip()

        print(f"  {gene} ({uniprot_id})", flush=True)

        sequence = parse_cerevisiae_sequence(orthologs_dir, gene)
        if sequence is None:
            print(f"    WARNING: no ortholog FASTA for {gene}, skipping")
            continue

        try:
            disorder = meta.predict_disorder(sequence)
        except Exception as exc:
            print(f"    ERROR: Metapredict failed for {gene}: {exc}")
            continue

        pdb_path = os.path.join(pdb_dir, f"{uniprot_id}.pdb")
        plddt_map = parse_pdb_plddt(pdb_dir, uniprot_id)

        sasa_map = {}
        if mkdssp_path and os.path.exists(pdb_path):
            sasa_map = run_dssp(pdb_path, mkdssp_path)
            if not sasa_map:
                print(f"    WARNING: DSSP failed for {uniprot_id}")

        for i, aa in enumerate(sequence):
            if aa not in ("S", "T", "Y"):
                continue
            pos = i + 1  # 1-indexed
            dis_score = round(float(disorder[i]), 3)
            dis_state = "Disordered" if dis_score >= 0.5 else "Ordered"
            plddt_val = plddt_map.get(pos)
            sasa_val = sasa_map.get(pos)
            is_psite = (gene, pos) in psite_set

            records.append({
                "gene": gene,
                "uniprot_id": uniprot_id,
                "position": pos,
                "residue": aa,
                "is_psite": is_psite,
                "disorder_score": dis_score,
                "disorder_state": dis_state,
                "plddt": plddt_val,
                "structural_state": structural_state(plddt_val),
                "sasa": sasa_val,
                "sasa_location": sasa_location(sasa_val),
            })

    result = pd.DataFrame(records)
    result.to_csv(out_csv, index=False)

    n_psite = int(result["is_psite"].sum())
    n_bg = int((~result["is_psite"]).sum())
    n_sasa = int(result["sasa"].notna().sum())
    print(f"\nSaved {len(result)} S/T/Y residues to {out_csv}")
    print(f"  P-sites   : {n_psite}")
    print(f"  Background: {n_bg}")
    print(f"  With SASA : {n_sasa}")


if __name__ == "__main__":
    main()
