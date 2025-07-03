"""
Script Name: GRNBoost2_AtSC.py

Purpose:
This script infers a gene regulatory network (GRN) using GRNBoost2 from expression data of Arabidopsis thaliana roots inoculated with synthetic microbial communities (SC).
It uses tree-based regression to predict regulatory edges between transcription factors (TFs) and their putative target genes.

Inputs:
- Excel expression matrix (TPM-normalized RNA-seq values) from Wippel et al. (2021)
  → Includes samples for Arabidopsis-SC (At-SC) and Lotus-SC (Lj-SC) conditions
- TF list in TSV format with Arabidopsis gene IDs (from PlantTFDB or JASPAR)
  
Output:
- A GRNBoost2-inferred network saved as a TSV file, with columns: regulator, target, importance score

Thesis Reference:
- Section 2.6: "Gene Regulatory Network Inference"
- Builds the SC-specific GRN used in Figures 6–8

Requirements:
- Python packages: pandas, numpy, dask, arboreto
- Dask installed and configured for parallelism
"""

import os
import pandas as pd
import numpy as np
from arboreto.algo import grnboost2
from dask.distributed import Client, LocalCluster
import logging

# Monkey patch for deprecated method in arboreto (if using older versions)
pd.DataFrame.as_matrix = lambda self: self.to_numpy()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # === Configuration ===
    excel_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Expression_data_At.xlsx"
    tf_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Ath_TF_list.txt"
    output_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output_AtSC_vs_LjSC_final_8_6_2025.tsv"

    # === Load expression matrix ===
    logging.info("Loading expression matrix...")
    df = pd.read_excel(excel_path)
    expr_cols = [col for col in df.columns if "C_AtSC" in col or "C_LjSC" in col]

    # Normalize and clean
    df[expr_cols] = df[expr_cols].replace(",", ".", regex=True).astype(float)
    df['Gene_ID'] = df['ID'].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df.set_index('Gene_ID', inplace=True)

    # Transpose for GRNBoost2 format: genes = columns, samples = rows
    expression_matrix = df[expr_cols].T
    expression_matrix.columns = expression_matrix.columns.str.replace(r"\.\d+$", "", regex=True)
    expression_matrix = expression_matrix.loc[:, ~expression_matrix.columns.duplicated()]
    expression_matrix = expression_matrix.loc[:, expression_matrix.std() > 0].dropna(axis=1)

    if expression_matrix.empty:
        logging.error("Expression matrix is empty after filtering.")
        return

    logging.info(f"Expression matrix shape: {expression_matrix.shape}")
    logging.info(f"Example genes: {list(expression_matrix.columns[:5])}")

    # === Load TF list ===
    logging.info("Loading TF list...")
    tf_df = pd.read_csv(tf_path, sep="\t")
    tf_df['Gene_ID'] = tf_df['Gene_ID'].astype(str).str.replace(r"\.\d+$", "", regex=True)

    matched_tfs = [tf for tf in tf_df['Gene_ID'] if tf in expression_matrix.columns]
    tf_variance_filtered = [tf for tf in matched_tfs if expression_matrix[tf].std() > 0]

    if not tf_variance_filtered:
        logging.error("No usable TFs found after filtering for expression and variability.")
        return

    logging.info(f"Variable TFs matched: {len(tf_variance_filtered)}")

    # === Run GRNBoost2 ===
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

    # === Save Output ===
    try:
        logging.info(f"Saving GRN to {output_path}...")
        network.to_csv(output_path, sep="\t", index=False)
        logging.info("GRNBoost2 run completed successfully.")
    except Exception as e:
        logging.error(f"Error saving GRN output: {e}")

if __name__ == "__main__":
    main()
