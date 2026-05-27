import pandas as pd
import requests
import metapredict as meta
import time

# Routes (Do after adding pLDDT data)
input_tsv = "../data/mt_nucleoid_PTMs_list_P-sites_processed.tsv" 
output_tsv = "../data/mt_nucleoid_PTMs_list_P-sites_processed.tsv"

# Load table
try:
    df = pd.read_csv(input_tsv, sep='\t')
except FileNotFoundError:
    print(f"{input_tsv} not found")
    exit(1)

meta_scores = []
meta_states = []

# Main cycle
for index, row in df.iterrows():
    uniprot_id = row['Uniprot ID']
    p_sites_raw = row['P-site positions']
    gene_name = row['Standard gene name']

    if pd.isna(uniprot_id) or pd.isna(p_sites_raw):
        meta_scores.append(float('nan'))
        meta_states.append(float('nan'))
        continue

    p_sites = []
    for x in str(p_sites_raw).split(','):
        try:
            p_sites.append(int(float(x.strip())))
        except ValueError:
            pass
    print(f"Analyzing: {gene_name} ({uniprot_id})...")

    # Download seuences
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        lines = response.text.strip().split('\n')
        sequence = "".join(lines[1:])

        try:
            # Analyse sequences via Metapredict V3
            predicted_disorder = meta.predict_disorder(sequence)

            scores_list = []
            states_list = []

            # Match to p-site
            for site in p_sites:
                if site <= len(sequence):
                    score = round(predicted_disorder[site - 1], 3)
                    scores_list.append(str(score))

                    # Threshold 0.5 (standard for Metapredict V3)
                    if score >= 0.5:
                        states_list.append("Disordered")
                    else:
                        states_list.append("Ordered")
                else:
                    scores_list.append("NA")
                    states_list.append("Unknown")

            meta_scores.append(",".join(scores_list))
            meta_states.append(",".join(states_list))

        except Exception as e:
            print(f"Failed to get score {gene_name}: {e}")
            meta_scores.append(float('nan'))
            meta_states.append(float('nan'))

    else:
        print(f"Failed to get sequence from Uniprot: {uniprot_id}")
        meta_scores.append(float('nan'))
        meta_states.append(float('nan'))

    time.sleep(0.2)

df['Metapredict Disorder Score'] = meta_scores
df['Metapredict State'] = meta_states

cols = list(df.columns)
if 'Annotation' in cols:
    cols.remove('Annotation')
    cols.append('Annotation')
    df = df[cols]

df.to_csv(output_tsv, sep='\t', index=False)
