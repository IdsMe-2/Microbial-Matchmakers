"""
Script Name: Gene_annotation.py

Purpose:
This script annotates FIMO motif‐scan hits with gene and genomic context:
  1. Parses FIMO output for motif occurrences in promoter regions.
  2. Retrieves the corresponding upstream promoter sequence for each hit.
  3. Looks up gene coordinates (chromosome, gene start/stop, strand) from a GFF3 annotation.
  4. Writes a combined CSV with motif ID, gene ID, hit position, sequence context, and genomic coordinates.

Inputs:
- FIMO TSV output (fimo.tsv) from MEME Suite scan of promoter FASTA.
- FASTA file of upstream sequences (e.g., 1 kb promoters).
- GFF3 file of gene models (TAIR10).

Output:
- CSV file with one row per FIMO hit, including promoter sequence and gene location.

Thesis Reference:
- Section 2.5 “Gene Annotation” (maps motif occurrences to genes for functional interpretation)
"""

import csv
from collections import defaultdict

# === File paths ===
FIMO_TSV         = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/fimo_Lj_At_homolog_cluster6_ALL_plant_motifs_output/fimo.tsv"
PROMOTER_FASTA   = "/home/15712745/personal/Gene_selection/TAIR10_upstream_1000_20101104.txt"
GFF3_FILE        = "/home/15712745/personal/Gene_selection/TAIR10_GFF3_genes.gff"
OUTPUT_CSV       = "/home/15712745/personal/TF_prediction_genomes/MEME/FIMO_folder/Gene_annotation_folder/annotated_fimo_results_Lotus_Cluster6.csv"

def normalize_gene_id(gene_id):
    """
    Strip version suffix from gene ID.
    E.g. 'AT1G01010.1' → 'AT1G01010'
    """
    return gene_id.split('.')[0]

# === Step 1: Parse FIMO hits ===
fimo_hits = []
with open(FIMO_TSV, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader, None)  # skip header
    for row in reader:
        if not row or row[0].startswith('#') or len(row) < 6:
            continue
        motif_id      = row[0]
        raw_gene_id   = row[2]
        start, stop   = row[3], row[4]
        strand        = row[5]
        gene_id       = normalize_gene_id(raw_gene_id)
        fimo_hits.append((motif_id, gene_id, start, stop, strand))

print(f"Parsed {len(fimo_hits)} FIMO hits.")

# === Step 2: Load promoter sequences into a dict ===
upstream_sequences = {}
with open(PROMOTER_FASTA, 'r', encoding='utf-8') as fasta:
    current_id = None
    seq_chunks = []
    for line in fasta:
        line = line.rstrip()
        if line.startswith('>'):
            # Save previous
            if current_id and seq_chunks:
                upstream_sequences[current_id] = ''.join(seq_chunks)
            # New record
            current_id = normalize_gene_id(line[1:].split()[0])
            seq_chunks = []
        else:
            seq_chunks.append(line)
    # Save last
    if current_id and seq_chunks:
        upstream_sequences[current_id] = ''.join(seq_chunks)

print(f"Loaded {len(upstream_sequences)} promoter sequences.")

# === Step 3: Load GFF3 gene annotations ===
gene_annotations = {}
with open(GFF3_FILE, 'r', encoding='utf-8') as gff:
    for line in gff:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9 or cols[2] != 'gene':
            continue
        chrom, source, feature, start, end, score, strand, phase, attrs = cols
        # Parse ID from attributes field
        try:
            # Assumes attribute string contains "ID=ATxGxxxxx;..."
            gene_id_raw = attrs.split('ID=')[1].split(';')[0]
        except IndexError:
            continue
        gene_id = normalize_gene_id(gene_id_raw)
        gene_annotations[gene_id] = {
            'chromosome': chrom,
            'gene_start': start,
            'gene_stop':  end,
            'strand':     strand
        }

print(f"Loaded annotations for {len(gene_annotations)} genes.")

# === Step 4: Write annotated output ===
with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as out:
    writer = csv.writer(out)
    writer.writerow([
        'Motif_ID','Gene_ID','Hit_Start','Hit_Stop','Hit_Strand',
        'Promoter_Sequence','Chromosome','Gene_Start','Gene_Stop','Gene_Strand'
    ])

    for motif_id, gene_id, start, stop, strand in fimo_hits:
        seq = upstream_sequences.get(gene_id, 'N/A')
        annot = gene_annotations.get(gene_id, {})
        writer.writerow([
            motif_id, gene_id, start, stop, strand,
            seq,
            annot.get('chromosome','N/A'),
            annot.get('gene_start','N/A'),
            annot.get('gene_stop','N/A'),
            annot.get('strand','N/A')
        ])

print(f"Annotated FIMO results saved to: {OUTPUT_CSV}")
