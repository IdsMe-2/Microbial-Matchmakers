import random
from Bio import SeqIO
import subprocess
import os
import pandas as pd
import sys
from collections import defaultdict
from statsmodels.stats.multitest import multipletests

# --- Config ---
input_fasta = "/home/15712745/personal/Gene_selection/selected_upstream_sequences_cluster3.fasta"
shuffled_fasta_base = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/shuffled_promoters_Arabidopsis_Cluster3_"
fimo_output_shuffled_base = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_shuffled_output_Arabidopsis_Cluster3_"
fimo_output_real = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_real_output_Arabidopsis_Cluster3"
motif_file = "/home/15712745/personal/TF_prediction_genomes/TF_bindingsite_motifs/ALL_plant_motifs_JASPAR.meme"
comparison_output = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_score_pval_comparison_Arabidopsis_Cluster3.csv"
num_shuffles = 100

# --- Functions ---
def shuffle_sequence(sequence):
    seq_list = list(sequence)
    random.shuffle(seq_list)
    return ''.join(seq_list)

def run_fimo(input_fasta, output_dir):
    fimo_cmd = f"fimo --oc {output_dir} {motif_file} {input_fasta}"
    try:
        subprocess.run(fimo_cmd, shell=True, check=True)
        print(f"FIMO run completed: results in {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"FIMO command failed: {e}")

def calculate_score_pvalues(real_fimo_tsv, shuffled_folder_base, num_shuffles, output_csv):
    real_df = pd.read_csv(real_fimo_tsv, sep='\t')
    real_scores = real_df.groupby('motif_id')['score'].sum()

    shuffled_scores_dict = defaultdict(list)

    for i in range(num_shuffles):
        shuffled_fimo_tsv = os.path.join(f"{shuffled_folder_base}{i}", "fimo.tsv")
        if os.path.exists(shuffled_fimo_tsv):
            df = pd.read_csv(shuffled_fimo_tsv, sep='\t')
            score_sums = df.groupby('motif_id')['score'].sum()
            for motif in set(score_sums.index).union(real_scores.index):
                shuffled_scores_dict[motif].append(score_sums.get(motif, 0))
        else:
            print(f"Warning: Missing shuffled FIMO result {i}")

    result_data = []
    for motif in real_scores.index:
        real_val = real_scores[motif]
        shuffled_vals = shuffled_scores_dict[motif]
        count_ge_real = sum(1 for val in shuffled_vals if val >= real_val)
        p_val = (count_ge_real + 1) / (len(shuffled_vals) + 1)
        mean_shuffled = sum(shuffled_vals) / len(shuffled_vals)

        result_data.append({
            'Motif_ID': motif,
            'Real_Score_Sum': real_val,
            'Shuffled_Mean_Score': mean_shuffled,
            'P_Value': p_val
        })

    # Create DataFrame and apply FDR correction
    result_df = pd.DataFrame(result_data)
    _, adj_pvals, _, _ = multipletests(result_df['P_Value'], method='fdr_bh')
    result_df['Adj_P_Value'] = adj_pvals

    result_df.to_csv(output_csv, index=False)
    print(f"Empirical p-values (with FDR correction) saved to {output_csv}")

# --- Main ---

# Run FIMO on real sequences
print("Running FIMO on real sequences...")
os.makedirs(fimo_output_real, exist_ok=True)
run_fimo(input_fasta, fimo_output_real)

# Run FIMO on shuffled sequences
for nr in range(num_shuffles):
    print(f"\nShuffling round {nr}")
    shuffled_fasta = f"{shuffled_fasta_base}{nr}.f"
    fimo_output_shuffled = f"{fimo_output_shuffled_base}{nr}"
    os.makedirs(fimo_output_shuffled, exist_ok=True)

    with open(shuffled_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            shuffled_seq = shuffle_sequence(record.seq)
            output_handle.write(f">{record.id}_shuffled\n{shuffled_seq}\n")
    print(f"Shuffled sequences saved to {shuffled_fasta}")

    run_fimo(shuffled_fasta, fimo_output_shuffled)

# Calculate empirical p-values (score comparison only)
print("\nComparing FIMO outputs using scores...")
real_fimo_tsv = os.path.join(fimo_output_real, "fimo.tsv")
calculate_score_pvalues(real_fimo_tsv, fimo_output_shuffled_base, num_shuffles, comparison_output)

print("Shuffled control complete")
