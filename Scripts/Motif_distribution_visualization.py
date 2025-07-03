"""
Script Name: Motif_distribution_visualization.py

Purpose:
This script performs motif enrichment analysis using Fisher's exact test. It compares motif occurrences
in a foreground gene cluster (e.g., Cluster 3 from Arabidopsis SC samples) versus a matched background set.
It visualizes significantly enriched or depleted motifs and corrects for multiple hypothesis testing.

Steps:
1. Count motif occurrences in FIMO outputs from foreground and background gene sets.
2. Build 2x2 contingency tables per motif and run Fisher's exact test.
3. Adjust p-values for multiple testing using Benjamini-Hochberg FDR.
4. Filter significantly enriched motifs and plot the odds ratios.
5. Save result tables and plots for downstream interpretation.

Inputs:
- FIMO output files (`fimo.tsv`) from MEME Suite for both foreground and background gene sets.

Outputs:
- CSV with enrichment statistics for each motif
- PNG barplot of enriched motifs
- TXT file with list of significantly enriched motif IDs

Thesis Reference:
- Section 2.4.3: Motif Enrichment Statistical Analysis (used in Figure 5)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# === File paths ===
foreground_path = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_At_cluster3_ALL_plant_motifs_output/fimo.tsv"
background_path = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_background_cluster3_03_06_2025/fimo.tsv"

output_csv = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_At_cluster3_ALL_plant_motifs_output/At_cluster3_enrichment_results_3_6_2025.csv"
output_plot = "/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/At_cluster3_motif_enrichment_plot_ALL_plant_promoters.png"
output_enriched_ids = "enriched_motifs_list.txt"

# === Load FIMO outputs ===
foreground = pd.read_csv(foreground_path, sep='\t')
background = pd.read_csv(background_path, sep='\t')

print("Unique motifs in foreground:", foreground['motif_id'].nunique())
print("Unique motifs in background:", background['motif_id'].nunique())

# === Count motif occurrences ===
foreground_counts = foreground['motif_id'].value_counts().reset_index()
foreground_counts.columns = ['motif_id', 'count_fg']

background_counts = background['motif_id'].value_counts().reset_index()
background_counts.columns = ['motif_id', 'count_bg']

# === Merge counts and fill missing values with 0 ===
merged = pd.merge(foreground_counts, background_counts, on='motif_id', how='outer').fillna(0)

# === Total hits ===
total_fg = len(foreground)
total_bg = len(background)

# === Fisher's Exact Test per motif ===
results = []
for _, row in merged.iterrows():
    a = int(row['count_fg'])  # motif in foreground
    b = int(row['count_bg'])  # motif in background
    c = total_fg - a          # non-motif in foreground
    d = total_bg - b          # non-motif in background

    oddsratio, pvalue = fisher_exact([[a, c], [b, d]])

    results.append({
        'Motif': row['motif_id'],
        'Foreground_Count': a,
        'Background_Count': b,
        'Odds_Ratio': oddsratio,
        'P_Value': pvalue
    })

enrichment_df = pd.DataFrame(results)

# === Multiple testing correction ===
enrichment_df['Adj_P_Value'] = multipletests(enrichment_df['P_Value'], method='fdr_bh')[1]

# === Filter significant motifs ===
significant = enrichment_df[enrichment_df['Adj_P_Value'] < 0.05]
significant = significant.sort_values('Odds_Ratio', ascending=False)

# === Save enriched motif IDs ===
enriched_motif_ids = significant['Motif'].tolist()
with open(output_enriched_ids, "w") as f:
    for motif in enriched_motif_ids:
        f.write(motif + "\n")

# === Save enrichment table ===
enrichment_df.to_csv(output_csv, index=False)
print(f"Enrichment results saved to: {output_csv}")
print(f"Enriched motif list saved to: {output_enriched_ids}")

# === Visualization ===
plt.figure(figsize=(10, 6))
sns.barplot(data=significant, x='Odds_Ratio', y='Motif', palette='viridis')
plt.axvline(1, color='red', linestyle='--')
plt.title("Significantly Enriched or Depleted TF Motifs in At Cluster 3")
plt.xlabel("Odds Ratio (Foreground vs Background)")
plt.ylabel("")
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
print(f"Plot saved to: {output_plot}")
