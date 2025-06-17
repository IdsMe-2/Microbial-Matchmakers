# Microbial Matchmakers – Regulatory Analysis of Host Microbiome Assembly

This repository contains the materials & methods, and supporting scripts used in a bioinformatics project to identify transcription factors involved in microbiome assembly in *Arabidopsis thaliana* and *Lotus japonicus*. The goal is to reconstruct gene regulatory networks that explain how plants modulate their microbiome in response to host-specific synthetic communities.

---

## Repository Structure

📦 microbial-matchmakers-methods/

┣ 📄 index.md ← Webpage content (Materials & Methods)

┣ 📄 README.md ← This file

┣ 📁 scripts/ ← All custom Python scripts used for analysis

┣ 📁 figures/ ← Figures used in the report and web documentation


## Contents

### Scripts
Includes all analysis steps:
- Gene selection & filtering
- Promoter extraction
- Motif discovery (FIMO/STREME)
- Motif enrichment validation
- Gene regulatory network construction (GRNBoost2)
- Network robustness testing

### Figures
Visuals from each stage of the pipeline, such as:
- TF motif enrichment
- GRN visualizations
- Perturbation analysis
- Bioinformatics flowchart

---

## Data Sources

- Arabidopsis genome: [TAIR10](https://www.arabidopsis.org/)
- Motif database: [JASPAR 2022](https://jaspar.genereg.net/)
- Raw expression data: Wippel et al., 2021
- Expression clustering source: [GitHub – Kathrin/Yulong Niu](https://github.com/YulongNiu/MPIPZ_Kathrin_Persistence_RNASeq)

---

## Citation

If using this repository or pipeline, please cite:
- Wippel et al., Nature Microbiology (2021)
- This GitHub repository (https://github.com/IdsMe-2/Microbial-Matchmakers.git) and it's corresponding thesis: *Microbial Matchmakers: Using Machine Learning to Identify Plant Transcription Factor Targets Involved in Host Microbiome Assembly*

---

## Author

**Ids Hartman**  
University of Amsterdam – Green Life Sciences (Biosystems Data Analysis & Plant Hormone Biology)  
June 2025  

