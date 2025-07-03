"""
Script Name: Shuffled_control_At_100_times.py

Purpose:
This script performs a full shuffled sequence control for motif enrichment in Arabidopsis Cluster 3:
1. Shuffles real promoter sequences 100 times
2. Runs FIMO on real and shuffled sets
3. Computes empirical p-values based on motif score sums
4. Applies multiple testing correction (FDR)
5. Visualizes significantly enriched motifs

Inputs:
- Real promoter FASTA file for Cluster 3
- JASPAR motif file (.meme format)
- FIMO must be installed and accessible

Outputs:
- Empirical p-values and FDR-adjusted CSV
- Filtered CSV of significant motifs
- Barplot and scatterplot visualizations

Thesis Reference:
- Section 2.4.4: "Statistical Controls Using Shuffling"
- Complements Fisher’s exact test (see Motif_distribution_visualization.py)
"""

import os
import random
import subprocess
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# === Config ===
input_fasta = "/home/15712745/personal/Gene_selection/selected_upstream_sequences_cluster3.fasta"
shuffled_fasta_base = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/shuffled_promoters_Arabidopsis_Cluster3_"
fimo_output_shuffled_base = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_shuffled_output_Arabidopsis_Cluster3_"
fimo_output_real = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_real_output_Arabidopsis_Cluster3"
motif_file = "/home/15712745/personal/TF_prediction_genomes/TF_bindingsite_motifs/ALL_plant_motifs_JASPAR.meme"
result_csv = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fimo_score_pval_comparison_Arabidopsis_Cluster3.csv"
filtered_csv = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fdr_significant_motifs_Cluster3.csv"
barplot_path = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fdr_corrected_enriched_motifs_barplot_7_6_2025.png"
scatterplot_path = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Shuffled_control/fdr_corrected_motif_score_scatter_7_6_2025.png"
num_shuffles = 100

# === Functions ===

def shuffle_sequence(sequence):
    """Return a randomly shuffled version of a sequence."""
    seq_list = list(sequence)
    random.shuffle(seq_list)
    return ''.join(seq_list)

def run_fimo(input_fasta, output_dir):
    """Run FIMO on a given input FASTA."""
    cmd = f"fimo --oc {output_dir} {motif_file} {input_fasta}"
    subprocess.run(cmd, shell=True, check=True)

def generate_shuffled_controls():
    """Generate shuffled promoter sequences and run FIMO."""
    os.makedirs(fimo_output_real, exist_ok=True)
    run_fimo(input_fasta, fimo_output_real)

    for i in range(num_shuffles):
        print(f"Shuffle round {i+1}/{num_shuffles}")
        shuffled_fasta = f"{shuffled_fasta_base}{i}.f"
        fimo_dir = f"{fimo_output_shuffled_base}{i}"
        os.makedirs(fimo_dir, exist_ok=True)

        with open(shuffled_fasta, "w") as out:
            for record in SeqIO.parse(input_fasta, "fasta"):
                shuffled = shuffle_sequence(record.seq)
                out.write(f">{record.id}_shuffled\n{shuffled}\n")

        run_fimo(shuffled_fasta, fimo_dir)

def compute_empirical_pvalues():
    """Calculate empirical p-values and adjust with FDR."""
    real_fimo = os.path.join(fimo_output_real, "fimo.tsv")
    real_df = pd.read_csv(real_fimo, sep='\t')
    real_scores = real_df.groupby('motif_id')['score'].sum()

    shuffled_scores = defaultdict(list)
    for i in range(num_shuffles):
        tsv = os.path.join(f"{fimo_output_shuffled_base}{i}", "fimo.tsv")
        if not os.path.exists(tsv):
            print(f"Missing FIMO shuffle file {i}")
            continue
        df = pd.read_csv(tsv, sep='\t')
        scores = df.groupby('motif_id')['score'].sum()
        for motif in set(scores.index).union(real_scores.index):
            shuffled_scores[motif].append(scores.get(motif, 0))

    results = []
    for motif in real_scores.index:
        real_score = real_scores[motif]
        shuf_scores = shuffled_scores[motif]
        p_val = (sum(s >= real_score for s in shuf_scores) + 1) / (len(shuf_scores) + 1)
        mean_shuffled = sum(shuf_scores) / len(shuf_scores)
        results.append({
            'motif_id': motif,
            'Real_Score_Sum': real_score,
            'Shuffled_Mean_Score': mean_shuffled,
            'P_Value': p_val
        })

    df = pd.DataFrame(results)
    df['Adjusted_P'] = multipletests(df['P_Value'], method='fdr_bh')[1]
    df['Score_Difference'] = df['Real_Score_Sum'] - df['Shuffled_Mean_Score']
    df['Significant'] = df['Adjusted_P'] < 0.05
    df.to_csv(result_csv, index=False)
    print(f"Empirical p-values saved to: {result_csv}")
    return df

def visualize_results(df):
    """Create barplot and scatterplot from results."""
    sig_df = df[(df['Score_Difference'] > 0) & (df['Significant'])].copy()
    sig_df.sort_values('Score_Difference', ascending=False).to_csv(filtered_csv, index=False)

    # Barplot
    sns.set(style="whitegrid")
    plt.figure(figsize=(20, 12))
    sns.barplot(
        x='Score_Difference',
        y='motif_id',
        data=sig_df,
        palette='viridis'
    )
    plt.title('FDR-Corrected Significantly Enriched Motifs (Cluster 3)', fontsize=18)
    plt.xlabel('Score Difference (Real - Shuffled)', fontsize=14)
    plt.ylabel('Motif ID', fontsize=14)
    plt.tight_layout()
    plt.savefig(barplot_path, dpi=300)
    plt.close()

    # Scatterplot
    df['Color_Label'] = df['Adjusted_P'].apply(lambda p: 'FDR < 0.05' if p < 0.05 else 'Not Significant')
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        x='Real_Score_Sum',
        y='Shuffled_Mean_Score',
        data=df,
        hue='Color_Label',
        palette={'FDR < 0.05': 'green', 'Not Significant': 'gray'},
        size='Score_Difference',
        sizes=(20, 300),
        alpha=0.7
    )
    plt.title('Real vs. Shuffled Motif Scores (Cluster 3)', fontsize=16)
    plt.xlabel('Real Score Sum')
    plt.ylabel('Shuffled Mean Score')
    plt.axline((0, 0), slope=1, linestyle='--', color='black')
    plt.legend(title='Significance (Adjusted)', loc='upper left')
    plt.tight_layout()
    plt.savefig(scatterplot_path, dpi=300)
    plt.close()

    print("Visualizations saved.")
    print(f"Filtered significant motifs saved to: {filtered_csv}")

# === Execute Pipeline ===
if __name__ == "__main__":
    print("Starting shuffled control analysis for Cluster 3...")
    generate_shuffled_controls()
    result_df = compute_empirical_pvalues()
    visualize_results(result_df)
    print("Full analysis complete.")
