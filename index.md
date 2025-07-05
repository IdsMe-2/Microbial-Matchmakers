# Microbial Matchmakers – Materials and Methods

This project explores how plants regulate microbiome assembly by identifying key transcription factors (TFs) and their gene targets in *Arabidopsis thaliana* and *Lotus japonicus*. Below is the full documentation of the bioinformatics pipeline, supporting scripts, and figures.

---

## Overview of Methods

### 1. Gene Selection

Gene expression clusters responsive to synthetic microbial communities (SCs) were extracted from Wippel et al. (2021). Clusters were filtered for transcription factors (TFs) using PlantTFDB annotations.
**Script used**: `background_cluster3.py`

### 2. Promoter Sequence Extraction

Arabidopsis promoter sequences were extracted (1000 bp upstream of TSS) using TAIR10. For Lotus, homologous Arabidopsis gene IDs were used due to limited genomic annotation.
**Scripts used**: `extract_upstream_promoter_sequences.py`, `background_cluster3.py`

### 3. Motif Discovery

* **FIMO**: scanned for known motifs from the JASPAR database
* **STREME**: discovered novel enriched motifs in promoter sequences
  **Scripts used**: `Gene_annotation.py`, `Shuffled_control_At_100_times.py`

### 4. Motif Validation and Enrichment

* GO enrichment using PANTHER
* Motif enrichment statistics using Fisher's exact test and FDR correction
* 100x shuffled motif controls
  **Scripts used**: `Motif_distribution_visualization.py`, `Shuffled_control_At_100_times.py`

### 5. Gene Regulatory Network (GRN) Inference

Used **GRNBoost2** (via Arboreto) on transcriptome data from Arabidopsis and Lotus to infer regulatory relationships between TFs and target genes.
**Scripts used**: `GRNBoost2_AtSC.py`, `GRNBoost2_global_GRN_At.py`

### 6. Network Visualization and Robustness Testing

Networks were visualized in Cytoscape. Perturbation analysis simulated TF removals to assess network stability.
**Scripts used**: `perturbation_analysis_part_1.py`, `perturbation_part_2_visualization_plot.py`

---

## Scripts

| Script Name                                                                                      | Description                                                                         |
| ------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| [`Gene_annotation.py`](scripts/Gene_annotation.py)                                               | Maps FIMO motifs to gene features using GFF3 annotations                            |
| [`GRNBoost2_global_GRN_At.py`](scripts/GRNBoost2_global_GRN_At.py)                               | Infers a global GRN using expression data across all samples                        |
| [`GRNBoost2_AtSC.py`](scripts/GRNBoost2_AtSC.py)                                                 | Infers a GRN using only At-SC samples to highlight context-specific regulation      |
| [`Motif_distribution_visualization.py`](scripts/Motif_distribution_visualization.py)             | Visualizes positional bias of TF motifs in promoters and computes enrichment        |
| [`Shuffled_control_At_100_times.py`](scripts/Shuffled_control_At_100_times.py)                   | Performs FIMO scans on 100 randomized promoter sets to calculate empirical p-values |
| [`Extract_upstream_promoter_sequences.py`](scripts/Extract_upstream_promoter_sequences.py)       | Extracts 1kb upstream sequences from Arabidopsis and Lotus genomes                  |
| [`Background_Arabidopsis_vs_Lotus.py`](scripts/Background_Arabidopsis_vs_Lotus.py)               | Creates background gene sets for motif enrichment analysis                          |
| [`perturbation_analysis_part_1.py`](scripts/perturbation_analysis_part_1.py)                     | Removes individual TFs from GRN and computes network fragmentation                  |
| [`perturbation_part_2_visualization_plot.py`](scripts/perturbation_part_2_visualization_plot.py) | Plots network robustness metrics from perturbation results                          |

---

## Figures

Each figure visualizes a key result or method step in the pipeline.

| Figure                                                                                                                                                                                                   | Description                                                    |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| !\[]\(figures/FigureS1 flowchart of microbial matchmakers.png)                                                                                                                                           | **Figure S1**: Full bioinformatics pipeline overview           |
| !\[]\(figures/Figure1 Transcriptional responses specific to syncoms in roots of lotus and arabidopsis, from Wippel et al 2021 with permission.png)                                                       | **Figure 1**: PCA and expression heatmaps from Wippel et al.   |
| !\[]\(figures/Figure2 comparative analysis of enriched TF binding motifs in arabidopsis cluster 3 and lotus cluster 6 via arabidopsis homologs.png)                                                      | **Figure 2**: Top TF motifs detected by FIMO                   |
| !\[]\(figures/Figure3 significantly enriched or depleted transcription factor motifs in arabdopsis cluster 3 promoters.png)                                                                              | **Figure 3**: Statistical enrichment of specific motifs        |
| !\[]\(figures/Figure4 real vs shuffled motif scores with FDR correction lotus cluster 6.png)                                                                                                             | **Figure 4**: Motif significance in Lotus via shuffled control |
| !\[]\(figures/Figure5 real vs shuffled motif scores with FDR correction arabidopsis cluster 3.png)                                                                                                       | **Figure 5**: Same as above for Arabidopsis                    |
| !\[]\(figures/Figure6 visual comparison of TFs between the global GRN and the SC specific GRN for arabidopsis.png)                                                                                       | **Figure 6**: Differential TF roles in global vs. SC GRNs      |
| !\[]\(figures/Figure7 gene regulatory network inferred from arabidopsis roots inoculated with synthetic communities visualized in cytoscape.png)                                                         | **Figure 7**: Cytoscape-rendered GRN of Arabidopsis            |
| !\[]\(figures/Figure8 gene regulatory subnetwork centered on the top three transcription factors with the highest number of predicted targets in the arabidopsis global GRN inferred with GRNBoost2.png) | **Figure 8**: Zoom-in on top TFs and their targets             |
| !\[]\(figures/Figure9 comparison of network disruption after individual removal of the top 100 most connected TFs from the SC specific GRN.png)                                                          | **Figure 9**: Effect of TF removals on network fragmentation   |
| !\[]\(figures/Figure10 impact of transcription factor removal on network connectivity.png)                                                                                                               | **Figure 10**: Quantitative summary of perturbation outcomes   |

---

## Data & Resources

* Arabidopsis genome: [TAIR10](https://www.arabidopsis.org/)
* Motif database: [JASPAR 2022](https://jaspar.genereg.net/)
* Expression data: From Wippel et al. (2021), Nature Microbiology
* Pipeline source repo: [MPIPZ\_Kathrin\_Persistence\_RNASeq (GitHub)](https://github.com/YulongNiu/MPIPZ_Kathrin_Persistence_RNASeq)

---

## Citation

If reusing this pipeline, please cite the original data sources and this GitHub repository.
