"""
Script Name: perturbation_analysis_part_1.py

Purpose:
This script performs a topological robustness analysis on the SC-specific Arabidopsis GRN. It iteratively removes
the top transcription factors (TFs) by degree and measures the impact on:
  - Number of components
  - Size of the largest component
  - Average shortest path length

Steps:
1. Load GRN inferred by GRNBoost2 and filter by importance threshold.
2. Rank TFs by their node degree in the directed graph.
3. Iteratively remove each top TF and recalculate graph metrics.
4. Save the resulting perturbation statistics to CSV.

Inputs:
- GRNBoost2 output TSV file with columns: TF, target, importance

Outputs:
- CSV file with perturbation metrics per TF (degree, #components, avg path length, largest CC)

Thesis Reference:
- Section 2.6.4: "Network Validation and Robustness Analysis" (used in Figure 9)
"""

import pandas as pd
import networkx as nx
import numpy as np
from tqdm import tqdm

# === Configuration ===
network_file = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output_AtSC_vs_LjSC_final_8_6_2025.tsv"
output_file = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/At_Network_perturbation_results_cutoff2_9_6_2025.csv"

importance_threshold = 2.0
top_n_tfs = 100

# === Step 1: Load and filter GRN ===
print("Loading GRNBoost2 output...")
df = pd.read_csv(network_file, sep="\t")
df = df[df["importance"] > importance_threshold]
print(f"Retained {len(df)} edges with importance > {importance_threshold}")

# === Step 2: Construct directed network ===
G_real = nx.from_pandas_edgelist(df, source="TF", target="target", create_using=nx.DiGraph())

# Ensure all TFs are included (some may only regulate others)
for tf in df["TF"].unique():
    if tf not in G_real:
        G_real.add_node(tf)

print(f"Network: {G_real.number_of_nodes()} nodes, {G_real.number_of_edges()} edges")

# === Step 3: Compute degree centrality and select top TFs ===
print("📈 Ranking TFs by degree...")
degree_dict = dict(G_real.degree())
top_tfs = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:top_n_tfs]
print(f"Selected top {top_n_tfs} TFs for perturbation.")

# === Step 4: Define metric computation ===
def compute_network_stats(G):
    stats = {}

    # Count number of weakly connected components
    stats["num_components"] = nx.number_weakly_connected_components(G)

    # Analyze the largest component
    largest_cc = max(nx.weakly_connected_components(G), key=len)
    subgraph = G.subgraph(largest_cc)

    try:
        stats["avg_path_length"] = nx.average_shortest_path_length(subgraph)
    except Exception:
        stats["avg_path_length"] = np.nan

    stats["largest_component_size"] = len(largest_cc)
    return stats

# === Step 5: Run perturbation analysis ===
print("Running TF removal simulations...")
perturbation_results = []

for tf, degree in tqdm(top_tfs):
    if tf not in G_real:
        print(f"Skipping TF {tf} — not found in graph.")
        continue

    G_copy = G_real.copy()
    G_copy.remove_node(tf)

    stats = compute_network_stats(G_copy)
    stats["Removed_TF"] = tf
    stats["Original_Degree"] = degree
    perturbation_results.append(stats)

# === Step 6: Save output ===
df_out = pd.DataFrame(perturbation_results)
df_out = df_out[["Removed_TF", "Original_Degree", "num_components", "avg_path_length", "largest_component_size"]]
df_out.to_csv(output_file, index=False)

print(f"\nPerturbation analysis complete. Results saved to:\n{output_file}")
