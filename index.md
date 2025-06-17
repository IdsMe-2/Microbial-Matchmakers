
# Materials and Methods

This section describes the bioinformatics pipeline used to identify transcription factor (TF) targets involved in host microbiome assembly in *Arabidopsis thaliana* and *Lotus japonicus*.

---

## 1. Gene Selection

RNA-seq data from Wippel et al. (2021) was used to identify host-specific clusters:
- Arabidopsis: Cluster 3
- Lotus: Cluster 6

Genes were selected based on co-expression under synthetic microbial communities (SCs). Only TFs from these clusters were used for further analysis.

---

## 2. Promoter Sequence Extraction

Promoters (1000 bp upstream of TSS) for Arabidopsis genes were extracted from TAIR10. Lotus homologs were mapped to Arabidopsis gene IDs when needed.

---

## 3. Motif Discovery

Two tools from MEME Suite were used:
- **FIMO**: scanned for known motifs (using JASPAR 2022)
- **STREME**: de novo motif discovery

---

## 4. Motif Validation and Enrichment

### 4.1 GO Enrichment
Functional GO enrichment was performed with PANTHER using genes from selected clusters.

### 4.2 Motif Enrichment
Motif occurrences in target genes vs background genes were compared using Fisher’s Exact Test with FDR correction.

### 4.3 Shuffled Control
Promoter sequences were shuffled (100x) and FIMO was rerun to generate an empirical p-value distribution.

---

## 5. Gene Annotation

Motif occurrences were mapped to Arabidopsis GFF3 annotations. Motif location and orientation were linked to gene IDs.

---

## 6. Gene Regulatory Network (GRN) Inference

**GRNBoost2** (via Arboreto) inferred TF-target relationships using bulk RNA-seq data. Two networks were built:
- Global GRN (all 16 samples)
- SC-specific GRN (only AtSC and LjSC samples)

---

## 7. Network Visualization

TFs and their target motifs were visualized using **Cytoscape**. Key nodes were identified using centrality metrics.

---

## 8. Robustness Testing

TFs were removed one-by-one from the GRNs. Changes in network fragmentation were recorded using NetworkX.

---

This pipeline enables the discovery of context-specific regulatory patterns in host microbiome interactions and is applicable to other plant systems.
