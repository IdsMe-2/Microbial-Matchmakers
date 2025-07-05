# Microbial Matchmakers – Regulatory Analysis of Host Microbiome Assembly

This repository contains the materials, methods, and analysis scripts used in a bioinformatics project to identify transcription factors (TFs) involved in microbiome assembly in *Arabidopsis thaliana* and *Lotus japonicus*. The goal is to reconstruct gene regulatory networks (GRNs) that explain how plants modulate their root microbiome in response to host-specific synthetic communities (SynComs).

---

## Project Overview

This project integrates:
- RNA-seq clustering and gene selection
- Promoter sequence extraction
- Motif discovery and statistical enrichment testing
- Shuffled-sequence statistical controls
- Gene regulatory network inference (GRNBoost2)
- Perturbation/robustness analysis of TFs

The results are presented in the thesis:  
**_"Microbial Matchmakers: Using Machine Learning to Identify Plant Transcription Factor Targets Involved in Host Microbiome Assembly"_**  
**Ids Hartman**, University of Amsterdam, June 2025

---

## Repository Structure


```txt
microbial-matchmakers-methods/
├── index.md                      # Materials & Methods (web version)
├── README.md                     # This file
├── scripts/                      # All custom Python scripts
├── figures/                      # PNGs used in report/thesis
├── data/                         # Sample inputs or Crunchomics references
├── requirements.txt              # Python environment dependencies

```
---

## Pipeline Overview


RNA-seq data → Cluster genes → Extract promoters → Motif scan with FIMO/STREME → GRN inference with GRNBoost2 → Network robustness by TF deletion


---

## Scripts and Methods (Linked to Thesis Sections)

| Script | Purpose | Thesis Section |
|--------|---------|----------------|
| `extract_upstream_promoter_sequences.py` | Extract 1kb upstream TAIR promoters | §2.2 |
| `background_cluster3.py` | Select background genes matched by expression | §2.4.1 |
| `Gene_annotation.py` | Map FIMO hits to gene location & metadata | §2.5 |
| `Motif_distribution_visualization.py` | Fisher’s exact test + FDR on motif counts | §2.4.3 |
| `Shuffled_control_analysis_Cluster3.py` | Empirical p-values from shuffled controls | §2.4.4 |
| `GRNBoost2_AtSC.py` | Build SC-specific Arabidopsis GRN | §2.6 |
| `GRNBoost2_global_GRN_At.py` | Build full-condition GRN | §2.6 |
| `perturbation_analysis_part_1.py` | Remove top TFs & compute network disruption | §2.6.4 |
| `perturbation_part_2_visualization_plot.py` | Plot impact of TF deletions | §2.6.4 |

---

## Example Outputs

- `At_cluster3_enrichment_results.csv` – motif enrichment table (Fisher + FDR)
- `fdr_significant_motifs_Cluster3.csv` – motifs with empirical p < 0.05
- `grnboost2_output.tsv` – inferred gene–TF edges with importance scores
- `TF_disruption_num_components_SC_specific.png` – TF perturbation impact plot

---

## How to Run

1. Clone and enter the repo:
   ```bash
   git clone https://github.com/IdsMe-2/Microbial-Matchmakers.git
   cd Microbial-Matchmakers

2. Create a Python environment:

   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. Run a script:

   ```bash
   python scripts/Motif_distribution_visualization.py
   ```

> Note: FIMO (from the MEME Suite) must be installed and accessible via your shell to run motif scanning steps.

---

## Python Dependencies

Install from `requirements.txt`:

```txt
pandas==2.2.2
numpy==1.26.4
matplotlib==3.8.4
seaborn==0.13.2
scipy==1.13.1
statsmodels==0.14.1
tqdm==4.66.4
networkx==3.3
biopython==1.83
dask==2024.4.1
arboreto==0.1.6
```

---

## Data Sources

| Resource                                                                                      | Description                               |
| --------------------------------------------------------------------------------------------- | ----------------------------------------- |
| [TAIR10](https://www.arabidopsis.org/)                                                        | Arabidopsis genome and upstream sequences |
| [JASPAR 2022](https://jaspar.genereg.net/)                                                    | TF binding motifs (Plant MEME-format)     |
| Wippel et al., 2021                                                                           | Raw RNA-seq data (Crunchomics)            |
| [Expression clustering GitHub](https://github.com/YulongNiu/MPIPZ_Kathrin_Persistence_RNASeq) | Gene clustering basis                     |

> Contact the author or your lab for access to restricted files on Crunchomics.

---

## Citation

If you use this repository or code, please cite:

> Wippel et al., *Nature Microbiology* (2021)

> Hartman, I. (2025) *Microbial Matchmakers: Using Machine Learning to Identify Plant Transcription Factor Targets Involved in Host Microbiome Assembly*

> GitHub: [https://github.com/IdsMe-2/Microbial-Matchmakers](https://github.com/IdsMe-2/Microbial-Matchmakers)

---

## Author

**Ids Hartman**
University of Amsterdam – Green Life Sciences
Focus: Biosystems Data Analysis & Plant Hormone Biology
June 2025

