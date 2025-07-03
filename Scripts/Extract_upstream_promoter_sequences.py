"""
Script Name: extract_upstream_promoter_sequences.py

Purpose:
This script extracts 1000 bp upstream promoter sequences for a user-defined set of Arabidopsis thaliana gene IDs 
from a reference FASTA file (e.g., TAIR10). It is used to generate custom promoter sets for motif discovery and enrichment analysis.

Inputs:
- A text file containing Arabidopsis gene IDs (with or without transcript versions, e.g. AT1G01010 or AT1G01010.1)
- A FASTA file containing upstream promoter sequences from TAIR
  (e.g., 1000 bp upstream sequences in FASTA format with headers starting with '>AT...')
  
Output:
- A FASTA file with promoter sequences corresponding to the input gene list

Thesis Reference:
- Described in Section 2.2 "Promoter Sequence Extraction"
- Used for generating input to FIMO and STREME motif analysis (Sections 2.3–2.4)
"""

import re

# === File Paths ===
gene_ids_file = "/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/Lotus_cluster6_background_arabidopsishomolog.txt"
tair_file = "/home/15712745/personal/Gene_selection/TAIR10_upstream_1000_20101104.txt"
output_file = "/home/15712745/personal/Gene_selection/selected_upstream_sequences_lotus_cluster6_Background_15_5_2025.fasta"

# === Step 1: Load Gene IDs ===
def get_gene_ids(txt_file):
    """
    Read and normalize gene IDs from a TXT file.
    Removes transcript version suffixes like ".1" if present.
    """
    gene_ids = set()
    with open(txt_file, "r", encoding="utf-8") as file:
        for line in file:
            match = re.match(r'(AT[1-5]G\d{5})(?:\.\d+)?', line.strip())
            if match:
                gene_ids.add(match.group(1))  # Extract only main gene ID, e.g. "AT1G01010"
    
    print(f"Extracted {len(gene_ids)} gene IDs (first 10 shown):")
    print(list(gene_ids)[:10])
    return gene_ids

# === Step 2: Extract Matching Promoter Sequences ===
def extract_sequences(tair_file, gene_ids):
    """
    Parse FASTA file and extract sequences for gene IDs in the input list.
    Assumes headers start with '>AT1G...' style IDs.
    """
    extracted_data = []
    found_genes = set()
    current_gene = None
    current_sequence = []

    with open(tair_file, "r", encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                if current_gene and current_gene in gene_ids:
                    extracted_data.append(f">{current_gene}\n" + "".join(current_sequence) + "\n")
                    found_genes.add(current_gene)
                
                # Parse new gene ID from header
                match = re.search(r'>(AT[1-5]G\d{5})', line)
                current_gene = match.group(1) if match else None
                current_sequence = []
            elif current_gene:
                current_sequence.append(line.strip())

        # Capture final entry if matched
        if current_gene and current_gene in gene_ids:
            extracted_data.append(f">{current_gene}\n" + "".join(current_sequence) + "\n")
            found_genes.add(current_gene)

    print(f"Found {len(found_genes)} matching promoter sequences (first 10 shown):")
    print(list(found_genes)[:10])
    return extracted_data

# === Step 3: Write Output ===
if __name__ == "__main__":
    gene_ids = get_gene_ids(gene_ids_file)
    selected_sequences = extract_sequences(tair_file, gene_ids)

    with open(output_file, "w") as output:
        output.writelines(selected_sequences)

    print(f"Extraction complete! {len(selected_sequences) // 2} genes saved to {output_file}")
