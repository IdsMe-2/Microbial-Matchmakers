import pandas as pd
import networkx as nx
import numpy as np
from tqdm import tqdm

# === CONFIGURATION ===
network_file = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output_AtSC_vs_LjSC_final_8_6_2025.tsv"
output_file = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/At_Network_perturbation_results_cutoff2_9_6_2025.csv"
importance_threshold = 2.0
top_n_tfs = 100

# === Load and Filter GRN ===
print("🔄 Loading GRNBoost2 network...")
df = pd.read_csv(network_file, sep="\t")

# Filter by importance threshold
df = df[df["importance"] > importance_threshold]
print(f"✅ Filtered edges with importance > {importance_threshold}: {len(df)} edges remaining.")

# Build directed graph
G_real = nx.from_pandas_edgelist(df, source="TF", target="target", create_using=nx.DiGraph())

# Ensure all TFs are nodes (in case some have no incoming edges)
for tf in df["TF"].unique():
    if tf not in G_real:
        G_real.add_node(tf)

print(f"Final network has {G_real.number_of_nodes()} nodes and {G_real.number_of_edges()} edges.")

# === Compute real centrality scores (to rank TFs) ===
print("Computing TF degrees...")
real_centrality = {
    "degree": dict(G_real.degree())
}

# === Select top N TFs based on degree ===
top_tfs = sorted(real_centrality["degree"].items(), key=lambda x: x[1], reverse=True)[:top_n_tfs]
print(f"Top {top_n_tfs} TFs selected for perturbation.")

# === Function to compute network stats ===
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

# === Run perturbation analysis ===
print("Running perturbation analysis...")
perturbation_results = []

for tf, _ in tqdm(top_tfs):
    if tf not in G_real:
        print(f"⚠️  Skipping TF {tf} — not found in graph.")
        continue

    G_copy = G_real.copy()
    G_copy.remove_node(tf)

    stats = compute_network_stats(G_copy)
    stats["Removed_TF"] = tf
    stats["Original_Degree"] = real_centrality["degree"].get(tf, 0)
    perturbation_results.append(stats)

# === Save Results ===
df_perturb = pd.DataFrame(perturbation_results)
df_perturb = df_perturb[["Removed_TF", "Original_Degree", "num_components", "avg_path_length", "largest_component_size"]]
df_perturb.to_csv(output_file, index=False)

print(f"\n✅ Perturbation analysis complete. Results saved to:\n{output_file}")
