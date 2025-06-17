import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv("/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/At_Network_perturbation_results_cutoff2_9_6_2025.csv")

# Sort by number of components (higher = more disruption)
df_sorted = df.sort_values("num_components", ascending=False)

# Highlight top 20 most disruptive TFs
top_disruptors = df_sorted.head(20)["Removed_TF"].tolist()
df_sorted["Highlight"] = df_sorted["Removed_TF"].apply(lambda x: "Top 20" if x in top_disruptors else "Other")

# Plot
plt.figure(figsize=(17, 6))
sns.barplot(data=df_sorted, x="Removed_TF", y="num_components", hue="Highlight", dodge=False,
            palette={"Top 20": "crimson", "Other": "lightgrey"})
plt.xticks(rotation=90)
plt.title("Number of Network Components After TF Removal")
plt.ylabel("Number of Weakly Connected Components")
plt.xlabel("Removed Transcription Factor")
plt.legend(title="TF Category")
plt.tight_layout()
plt.savefig("/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Network_validation_and_robustness_of_network/TF_disruption_num_components_SC_specific_9_6_2025.png", dpi=300)
plt.show()
