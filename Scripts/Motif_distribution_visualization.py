import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Load both datasets
foreground = pd.read_csv("/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_At_cluster3_ALL_plant_motifs_output/fimo.tsv", sep='\t')
background = pd.read_csv("/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_background_cluster3_03_06_2025/fimo.tsv", sep='\t')

print("Unique motifs in foreground:", foreground['motif_id'].nunique())
print("Unique motifs in background:", background['motif_id'].nunique())


# Get motif counts
foreground_counts = foreground['motif_id'].value_counts().reset_index()
foreground_counts.columns = ['motif_id', 'count_fg']

background_counts = background['motif_id'].value_counts().reset_index()
background_counts.columns = ['motif_id', 'count_bg']

# Merge counts
merged = pd.merge(foreground_counts, background_counts, on='motif_id', how='outer').fillna(0)

# Total counts
total_fg = len(foreground)
total_bg = len(background)

# Calculate enrichment for all motifs
enrichment = []
for _, row in merged.iterrows():
    a = int(row['count_fg'])
    b = int(row['count_bg'])
    c = total_fg - a
    d = total_bg - b

    # 2x2 contingency table: [[a, c], [b, d]]
    oddsratio, pvalue = fisher_exact([[a, c], [b, d]])

    enrichment.append({
        'Motif': row['motif_id'],
        'Odds_Ratio': oddsratio,
        'P_Value': pvalue,
        'Foreground_Count': a,
        'Background_Count': b
    })

# Convert to DataFrame
enrichment_df = pd.DataFrame(enrichment)

# Multiple testing correction
enrichment_df['Adj_P_Value'] = multipletests(enrichment_df['P_Value'], method='fdr_bh')[1]

# Filter significant motifs (adj. p < 0.05)
significant = enrichment_df[enrichment_df['Adj_P_Value'] < 0.05].sort_values('Odds_Ratio', ascending=False)

# Save enriched motif IDs to a file
enriched_motif_ids = significant['Motif'].tolist()
with open("enriched_motifs_list.txt", "w") as f:
    for motif in enriched_motif_ids:
        f.write(motif + "\n")

# Save the enrichment_df to a CSV file
enrichment_df.to_csv("/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_At_cluster3_ALL_plant_motifs_output/At_cluster3_enrichment_results_3_6_2025.csv", index=False)

# Visualization
plt.figure(figsize=(10, 6))
sns.barplot(data=significant, x='Odds_Ratio', y='Motif', palette='viridis')
plt.axvline(1, color='red', linestyle='--')
plt.title("Significantly Enriched or Depleted TF Motifs in At Cluster 3")
plt.xlabel("Odds Ratio (Foreground vs Background)")
plt.ylabel("")
plt.savefig("/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/At_cluster3_motif_enrichment_plot_ALL_plant_promoters.png", dpi=300, bbox_inches='tight')

print("Plot saved as At_Cluster3_motif_enrichment_plot_ALL_plant_promoters.png")
print("Enrichment results saved as At_Cluster3_enrichment_results.csv")
