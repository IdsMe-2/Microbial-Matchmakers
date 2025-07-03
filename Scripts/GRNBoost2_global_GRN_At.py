"""
Script Name: GRNBoost2_global_GRN_At.py

Purpose:
This script infers a *global* gene regulatory network (GRN) using GRNBoost2 from full expression data of Arabidopsis thaliana roots under all treatment conditions. 
Unlike the SC-specific GRN, this network captures broader transcriptional regulatory architecture across At and Lj samples.

Inputs:
- Excel expression matrix containing all 16 samples (SC and control)
- TSV file of transcription factors (Arabidopsis TFs from PlantTFDB/JASPAR)

Output:
- TSV file with inferred GRN edges: columns = ['regulator', 'target', 'importance']

Thesis Reference:
- Section 2.6: Global GRN construction
- Used in comparison with SC-specific GRN in Figures 6, 9, and 11
"""

import os
import pandas as pd
from arboreto.algo import grnboost2
from dask.distributed import Client
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # === Configuration ===
    excel_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Expression_data_At.xlsx"
    tf_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Ath_TF_list.txt"
    output_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output.tsv"

    # === Check for input files ===
    if not os.path.exists(excel_path):
        logging.error(f"Expression file not found: {excel_path}")
        return
    if not os.path.exists(tf_path):
        logging.error(f"TF list file not found: {tf_path}")
        return

    # === Load expression matrix ===
    logging.info("Loading expression matrix...")
    try:
        df = pd.read_excel(excel_path)
        df.iloc[:, 9:-1] = df.iloc[:, 9:-1].replace(",", ".", regex=True).astype(float)
        df['Gene_ID'] = df['ID'].astype(str).str.replace(r"\.\d+$", "", regex=True)
        df.set_index('Gene_ID', inplace=True)

        # Transpose for GRNBoost2 (samples as rows, genes as columns)
        expression_matrix = df.iloc[:, 9:-1].T
        expression_matrix.columns.name = None
        expression_matrix.columns = expression_matrix.columns.str.replace(r'\.\d+$', '', regex=True)

        if expression_matrix.empty:
            logging.error("Expression matrix is empty after processing.")
            return
    except Exception as e:
        logging.error(f"Error loading or processing expression matrix: {e}")
        return

    # === Load transcription factors ===
    logging.info("Loading transcription factor list...")
    try:
        tf_df = pd.read_csv(tf_path, sep="\t")
        if tf_df.empty:
            logging.error("TF list file is empty.")
            return

        tf_names = [tf for tf in tf_df['Gene_ID'] if tf in expression_matrix.columns]
        if not tf_names:
            logging.warning("No TFs found in expression matrix.")
            return

        logging.info(f"TFs matched in expression data: {len(tf_names)}")
    except Exception as e:
        logging.error(f"Error loading TF list: {e}")
        return

    # === Start Dask and run GRNBoost2 ===
    logging.info("Starting Dask client...")
    try:
        client = Client()
    except Exception as e:
        logging.error(f"Failed to start Dask client: {e}")
        return

    logging.info("Running GRNBoost2 (global)...")
    try:
        network = grnboost2(expression_data=expression_matrix, tf_names=tf_names)
        if network.empty:
            logging.warning("GRNBoost2 returned an empty network.")
    except Exception as e:
        logging.error(f"GRNBoost2 failed: {e}")
        return

    # === Save output ===
    try:
        logging.info(f"Saving GRN to {output_path}...")
        network.to_csv(output_path, sep="\t", index=False)
        logging.info("GRNBoost2 global network saved.")
    except Exception as e:
        logging.error(f"Failed to save output: {e}")
        return

if __name__ == "__main__":
    main()
