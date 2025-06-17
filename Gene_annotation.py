import csv
from collections import defaultdict

# Input files
fimo_output_file = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_Lj_At_homolog_cluster6_ALL_plant_motifs_output/fimo.tsv"
upstream_sequences_file = "/home/15712745/personal/Gene_selection/TAIR10_upstream_1000_20101104.txt"
gff_file = "/home/15712745/personal/Gene_selection/TAIR10_GFF3_genes.gff"

# Output file
output_file = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Gene_annotation_folder/annotated_fimo_results_Lotus_Cluster6.csv"

def normalize_gene_id(gene_id):
    """Remove version suffix from gene ID (e.g. AT1G01010.1 → AT1G01010)."""
    return gene_id.split('.')[0] if '.' in gene_id else gene_id

# Step 1: Parse FIMO hits
fimo_hits = []
with open(fimo_output_file, 'r') as fimo_file:
    reader = csv.reader(fimo_file, delimiter='\t')
    next(reader)  # Skip header
    for row in reader:
        if not row or row[0].startswith('#'):
            continue
        if len(row) < 6:
            continue
        motif_id, raw_gene_id, start, stop, strand = row[0], row[2], row[3], row[4], row[5]
        norm_gene_id = normalize_gene_id(raw_gene_id)
        fimo_hits.append((motif_id, norm_gene_id, start, stop, strand))

# Step 2: Load upstream sequences
upstream_sequences = {}
with open(upstream_sequences_file, 'r') as fasta:
    current_gene = None
    sequence = []
    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            if current_gene:
                upstream_sequences[normalize_gene_id(current_gene)] = ''.join(sequence)
            current_gene = line.split()[0][1:]
            sequence = []
        else:
            sequence.append(line)
    if current_gene:
        upstream_sequences[normalize_gene_id(current_gene)] = ''.join(sequence)

# Step 3: Load GFF annotations
gene_annotations = {}
with open(gff_file, 'r') as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 9 or cols[2] != "gene":
            continue
        try:
            attrs = cols[8]
            gene_id_raw = attrs.split("ID=")[1].split(";")[0]
            norm_gene_id = normalize_gene_id(gene_id_raw)
            gene_annotations[norm_gene_id] = {
                "chromosome": cols[0],
                "start": cols[3],
                "stop": cols[4],
                "strand": cols[6]
            }
        except Exception as e:
            print(f"Error parsing GFF line: {line}\nError: {e}")

# Step 4: Output annotated results
with open(output_file, 'w', newline='') as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Motif_ID", "Gene_ID", "Start", "Stop", "Strand", "Upstream_Sequence", "Chromosome", "Gene_Start", "Gene_Stop"])

    for motif_id, gene_id, start, stop, strand in fimo_hits:
        upstream_seq = upstream_sequences.get(gene_id, "N/A")
        annotation = gene_annotations.get(gene_id, {})
        writer.writerow([
            motif_id, gene_id, start, stop, strand,
            upstream_seq,
            annotation.get("chromosome", "N/A"),
            annotation.get("start", "N/A"),
            annotation.get("stop", "N/A")
        ])

print(f"Annotated FIMO results saved to: {output_file}")
