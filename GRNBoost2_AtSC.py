import os
import pandas as pd
import numpy as np
from arboreto.algo import grnboost2
from dask.distributed import Client, LocalCluster
import logging

# Monkey patch for arboreto compatibility
pd.DataFrame.as_matrix = lambda self: self.to_numpy()

# ========== CONFIG LOGGING ==========
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # ========== CONFIG ==========
    excel_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Expression_data_At.xlsx"
    tf_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Ath_TF_list.txt"
    output_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output_AtSC_vs_LjSC_final_8_6_2025.tsv"

    # ========== LOAD EXPRESSION ==========
    logging.info("Loading expression matrix...")
    df = pd.read_excel(excel_path)
    expr_cols = [col for col in df.columns if "C_AtSC" in col or "C_LjSC" in col]

    df[expr_cols] = df[expr_cols].replace(",", ".", regex=True).astype(float)
    df['Gene_ID'] = df['ID'].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df.set_index('Gene_ID', inplace=True)

    expression_matrix = df[expr_cols].T
    expression_matrix.columns = expression_matrix.columns.str.replace(r"\.\d+$", "", regex=True)
    expression_matrix.columns = expression_matrix.columns.astype(str)
    expression_matrix = expression_matrix.loc[:, ~expression_matrix.columns.duplicated()]

    # Drop non-variable and NaN genes
    expression_matrix = expression_matrix.loc[:, expression_matrix.std() > 0]
    expression_matrix = expression_matrix.dropna(axis=1)

    if expression_matrix.empty:
        logging.error("Expression matrix is empty after filtering.")
        return

    logging.info(f"Expression matrix shape: {expression_matrix.shape}")
    logging.info(f"First genes: {list(expression_matrix.columns[:5])}")

    # ========== LOAD TFs ==========
    logging.info("Loading TF list...")
    tf_df = pd.read_csv(tf_path, sep="\t")
    tf_df['Gene_ID'] = tf_df['Gene_ID'].astype(str).str.replace(r"\.\d+$", "", regex=True)

    matched_tfs = [tf for tf in tf_df['Gene_ID'] if tf in expression_matrix.columns]
    tf_variance_filtered = [
    tf for tf in matched_tfs 
    if tf in expression_matrix.columns and expression_matrix[tf].std().item() > 0
]


    if not tf_variance_filtered:
        logging.error("No usable TFs found after filtering for expression and variability.")
        return

    logging.info(f"TFs matched and variable: {len(tf_variance_filtered)}")

    # ========== RUN GRNBoost2 ==========
    logging.info("Starting Dask client...")
    cluster = LocalCluster(n_workers=4, threads_per_worker=1)
    client = Client(cluster)

    try:
        logging.info("Running GRNBoost2...")
        network = grnboost2(expression_data=expression_matrix, tf_names=tf_variance_filtered)
        if network.empty:
            logging.warning("GRNBoost2 returned an empty network.")
    except Exception as e:
        logging.error(f"GRNBoost2 execution failed: {e}")
        return

    # ========== SAVE ==========
    try:
        logging.info(f"Saving output to {output_path}...")
        network.to_csv(output_path, sep="\t", index=False)
        logging.info("✅ GRNBoost2 run completed.")
    except Exception as e:
        logging.error(f"Error saving output: {e}")

if __name__ == "__main__":
    main()
