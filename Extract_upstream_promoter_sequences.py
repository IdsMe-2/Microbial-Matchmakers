import re

# File paths
gene_ids_file = "/home/15712745/personal/TF_prediction_genomes/MEME/Visualization_MEME/Lotus_cluster6_background_arabidopsishomolog.txt"
tair_file = "/home/15712745/personal/Gene_selection/TAIR10_upstream_1000_20101104.txt"
output_file = "/home/15712745/personal/Gene_selection/selected_upstream_sequences_lotus_cluster6_Background_15_5_2025.fasta"

# Function to extract gene IDs from the TXT file
def get_gene_ids(txt_file):
    gene_ids = set()
    
    with open(txt_file, "r", encoding="utf-8") as file:
        for line in file:
            match = re.match(r'(AT[1-5]G\d{5})(?:\.\d+)?', line.strip())  # Match ID, remove version
            if match:
                gene_ids.add(match.group(1))  # Store without version number
    
    print(f"Extracted {len(gene_ids)} gene IDs from TXT file (first 10 shown):")
    print(list(gene_ids)[:10])  # Debug: Print first 10 IDs
    
    return gene_ids

# Function to extract sequences from the TAIR file
def extract_sequences(tair_file, gene_ids):
    extracted_data = []
    found_genes = set()
    
    with open(tair_file, "r", encoding="utf-8") as file:
        current_gene = None
        current_sequence = []
        
        for line in file:
            if line.startswith(">"):  # New gene entry
                if current_gene and current_gene in gene_ids:
                    extracted_data.append(f">{current_gene}\n" + "".join(current_sequence) + "\n")
                    found_genes.add(current_gene)
                
                # Extract gene ID from header
                match = re.search(r'>(AT[1-5]G\d{5})', line)
                if match:
                    current_gene = match.group(1)  # Extract only gene ID
                    current_sequence = []
                else:
                    current_gene = None
            elif current_gene:
                current_sequence.append(line.strip())

        # Save last gene if it matches
        if current_gene and current_gene in gene_ids:
            extracted_data.append(f">{current_gene}\n" + "".join(current_sequence) + "\n")
            found_genes.add(current_gene)

    print(f"Found {len(found_genes)} matching genes in TAIR file (first 10 shown):")
    print(list(found_genes)[:10])  # Debug: Print first 10 matched IDs
    
    return extracted_data

# Main execution
if __name__ == "__main__":
    gene_ids = get_gene_ids(gene_ids_file)
    selected_sequences = extract_sequences(tair_file, gene_ids)

    # Write to output file
    with open(output_file, "w") as output:
        output.writelines(selected_sequences)

    print(f"Extraction complete! {len(selected_sequences) // 2} genes saved to {output_file}")
