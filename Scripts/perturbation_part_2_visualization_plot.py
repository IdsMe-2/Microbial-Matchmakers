"""
Script Name: perturbation_part_2_visualization_plot.py

Purpose:
This script visualizes the effect of removing individual transcription factors (TFs) from the SC-specific Arabidopsis gene regulatory network (GRN).
It plots the number of weakly connected components remaining after each TF removal, highlighting the top 20 most disruptive TFs.

Steps:
1. Load perturbation analysis results from part 1.
2. Rank TFs by the number of components caused after their removal.
3. Highlight top 20 most disruptive TFs in the barplot.
4. Save and display the plot.

Input:
- CSV file from `perturbation_analysis_part_1.py` containing per-TF network disruption metrics

Output:
- PNG barplot showing network fragmentation caused by each TF

Thesis Reference:
- Section 2.6.4: "Network Validation and Robustness Analysis" (Figure 9)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Configuration ===
input_csv = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/At_Network_perturbation_results_cutoff2_9_6_2025.csv"
output_plot = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/TF_disruption_num_components_SC_specific_9_6_2025.png"

# === Load data ===
df = pd.read_csv(input_csv)

# === Sort and annotate ===
df_sorted = df.sort_values("num_components", ascending=False)
top_disruptors = df_sorted.head(20)["Removed_TF"].tolist()
df_sorted["Highlight"] = df_sorted["Removed_TF"].apply(lambda x: "Top 20" if x in top_disruptors else "Other")

# === Plot ===
sns.set(style="whitegrid")
plt.figure(figsize=(17, 6))
sns.barplot(
    data=df_sorted,
    x="Removed_TF",
    y="num_components",
    hue="Highlight",
    dodge=False,
    palette={"Top 20": "crimson", "Other": "lightgrey"}
)

plt.xticks(rotation=90, fontsize=10)
plt.yticks(fontsize=12)
plt.title("Impact of TF Removal on GRN Fragmentation (SC-specific Network)", fontsize=16)
plt.ylabel("Number of Weakly Connected Components", fontsize=13)
plt.xlabel("Removed Transcription Factor", fontsize=13)
plt.legend(title="TF Category", loc='upper right')
plt.tight_layout()

# === Save and show ===
plt.savefig(output_plot, dpi=300)
plt.show()
print(f"Plot saved to: {output_plot}")
