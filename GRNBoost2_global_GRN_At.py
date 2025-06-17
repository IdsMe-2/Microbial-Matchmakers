import os
import pandas as pd
from arboreto.algo import grnboost2
from dask.distributed import Client
import logging

# ========== CONFIG LOGGING ==========
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # ========== CONFIG ==========
    excel_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Expression_data_At.xlsx"
    tf_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/Ath_TF_list.txt"
    output_path = "/home/15712745/personal/TF_prediction_genomes/Gene_regulatory_network/grnboost2_output.tsv"

    # ========== CHECK FILES EXIST ==========
    if not os.path.exists(excel_path):
        logging.error(f"Expression file not found: {excel_path}")
        return
    if not os.path.exists(tf_path):
        logging.error(f"TF list file not found: {tf_path}")
        return

    # ========== LOAD DATA ==========
    logging.info("Loading expression matrix...")
    try:
        df = pd.read_excel(excel_path)
    except Exception as e:
        logging.error(f"Error loading expression file: {e}")
        return

    try:
        # Clean and convert expression values
        df.iloc[:, 9:-1] = df.iloc[:, 9:-1].replace(",", ".", regex=True)
        df.iloc[:, 9:-1] = df.iloc[:, 9:-1].astype(float)

        df['Gene_ID'] = df['ID']
        df.set_index('Gene_ID', inplace=True)

        expression_matrix = df.iloc[:, 9:-1].T
        expression_matrix.columns.name = None
        expression_matrix.columns = expression_matrix.columns.str.replace(r'\.\d+$', '', regex=True)

        if expression_matrix.empty:
            logging.error("Expression matrix is empty after processing.")
            return
    except Exception as e:
        logging.error(f"Error processing expression matrix: {e}")
        return

    # ========== LOAD TF LIST ==========
    logging.info("Loading TF list...")
    try:
        tf_df = pd.read_csv(tf_path, sep="\t")
        if tf_df.empty:
            logging.error("TF list file is empty.")
            return

        tf_names = [tf for tf in tf_df['Gene_ID'] if tf in expression_matrix.columns]
        if not tf_names:
            logging.warning("No TFs found in the expression matrix.")
            return

        logging.info(f"Number of TFs matched in data: {len(tf_names)}")
    except Exception as e:
        logging.error(f"Error loading TF list: {e}")
        return

    # ========== RUN GRNBOOST2 ==========
    logging.info("Starting Dask client...")
    try:
        client = Client()
    except Exception as e:
        logging.error(f"Failed to start Dask client: {e}")
        return

    logging.info("Running GRNBoost2...")
    try:
        network = grnboost2(expression_data=expression_matrix, tf_names=tf_names)
        if network.empty:
            logging.warning("GRNBoost2 output is empty. Check input data or TF list.")
    except Exception as e:
        logging.error(f"GRNBoost2 execution failed: {e}")
        return

    # ========== SAVE RESULTS ==========
    try:
        logging.info(f"Saving results to {output_path}...")
        network.to_csv(output_path, sep="\t", index=False)
    except Exception as e:
        logging.error(f"Error saving GRNBoost2 output: {e}")
        return

    logging.info("GRNBoost2 completed successfully.")

if __name__ == "__main__":
    main()
