import os
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from Bio import SeqIO
import subprocess
from datetime import datetime

# === Paths ===
EXPR_PATH = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Expression_data_At.xlsx"
UPSTREAM_FASTA = "/home/15712745/personal/Gene_selection/TAIR10_upstream_1000_20101104.txt"
MOTIF_FILE = "/home/15712745/personal/TF_prediction_genomes/TF_bindingsite_motifs/ALL_plant_motifs_JASPAR.meme"
FIMO_OUTPUT_DIR = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder"
OUTPUT_PREFIX = "background_cluster3_" + datetime.today().strftime("%d_%m_%Y")

GENE_ID_OUTPUT = f"/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/{OUTPUT_PREFIX}_gene_ids.txt"
PROMOTER_FASTA_OUTPUT = f"/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/{OUTPUT_PREFIX}_promoters.fasta"
FIMO_RUN_DIR = os.path.join(FIMO_OUTPUT_DIR, f"fimo_{OUTPUT_PREFIX}")

# === Step 1: Load Expression Data ===
df = pd.read_excel(EXPR_PATH)
df = df.replace(",", ".", regex=True)
df.iloc[:, 9:-1] = df.iloc[:, 9:-1].astype(float)
df['avg_expr'] = df.iloc[:, 9:-1].mean(axis=1)

# === Step 2: Match Background Genes ===
foreground_df = df[df['cl'] == 3].copy()
background_pool = df[df['cl'] != 3].copy()

if len(foreground_df) == 0 or len(background_pool) == 0:
    raise ValueError("Foreground or background pool is empty — check clustering column or input file.")

nn = NearestNeighbors(n_neighbors=1)
nn.fit(background_pool[['avg_expr']])
_, indices = nn.kneighbors(foreground_df[['avg_expr']])

matched_background_df = background_pool.iloc[indices.flatten()].copy()
matched_background_df['matched_to'] = foreground_df['ID'].values
background_ids = matched_background_df['ID'].dropna().astype(str).unique()

print(f"✅ Background genes selected: {len(background_ids)}")

# === Step 3: Write Background Gene IDs ===
with open(GENE_ID_OUTPUT, "w") as f:
    for gene_id in background_ids:
        f.write(gene_id + "\n")
print(f"✅ Written gene list to: {GENE_ID_OUTPUT}")

# === Step 4: Filter Promoter Sequences ===
background_ids_set = set(gene_id.split('.')[0] for gene_id in background_ids)
records = list(SeqIO.parse(UPSTREAM_FASTA, "fasta"))
filtered_records = [rec for rec in records if rec.id.split('.')[0] in background_ids_set]

if len(filtered_records) == 0:
    print("⚠️ No matching promoter sequences found. Check FASTA headers and gene ID format.")
else:
    print(f"🔍 Example FASTA headers:\n" + "\n".join(f">{r.id}" for r in filtered_records[:5]))
    print(f"✅ {len(filtered_records)} promoter sequences matched out of {len(background_ids)} genes")

SeqIO.write(filtered_records, PROMOTER_FASTA_OUTPUT, "fasta")
print(f"✅ Filtered promoter FASTA written to: {PROMOTER_FASTA_OUTPUT}")

# === Step 5: Run FIMO ===
os.makedirs(FIMO_RUN_DIR, exist_ok=True)
fimo_cmd = ["fimo", "--oc", FIMO_RUN_DIR, "--verbosity", "1", MOTIF_FILE, PROMOTER_FASTA_OUTPUT]

print("🔍 Running FIMO...")
subprocess.run(fimo_cmd, check=True)
print(f"✅ Done. FIMO output written to: {os.path.join(FIMO_RUN_DIR, 'fimo.tsv')}")
