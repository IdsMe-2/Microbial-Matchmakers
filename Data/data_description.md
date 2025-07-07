# Data Description – Microbial Matchmakers

This file describes the contents of the `data/` folder, providing clarification on each input table or sequence file used in the thesis project.

---

## 1. `Expression_data_At.xlsx`

**Description:**

* Contains normalized gene expression values (FPKM) for *Arabidopsis thaliana*.
* Source: Wippel et al., 2021, *Nature Microbiology*.
* Includes **16 samples**, covering **4 experimental conditions**:

| Sample Prefix | Condition Description                       |
| ------------- | ------------------------------------------- |
| `C_mock_*`    | Sterile control (no SynCom)                 |
| `C_AtSC_*`    | Inoculated with Arabidopsis SynCom          |
| `C_LjSC_*`    | Inoculated with Lotus SynCom                |
| `C_fSC_*`     | Mixed SynCom (combined Arabidopsis + Lotus) |

**Columns:**

* Metadata columns (e.g. `ID`, gene annotation)
* `cl` – gene cluster label (as defined in the original study)
* Expression values for each of the 16 samples

**Note:**

* While the **heatmap** in the thesis/report reflects only the first 3 conditions (Mock, AtSC, LjSC), the **fSC** (mixed SynCom) condition is included in the expression matrix and used in GRN inference and motif discovery analyses.

---

## 2. `selected_upstream_sequences_cluster3.fasta`

**Description:**

* Contains 1000 bp upstream promoter sequences for genes in Cluster 3 of *Arabidopsis*.
* FASTA headers correspond to TAIR10 gene IDs.
* These sequences were scanned with known TF motifs using **FIMO**.

---

## 3. `selected_upstream_sequences_lotus_cluster6_Background_15_5_2025.fasta`

**Description:**

* Promoter sequences (Arabidopsis homologs) corresponding to *Lotus* Cluster 6 background genes.
* Used for comparative motif analysis between species.

---

## 4. Other Notes

* All expression matrices and promoter sets are reproducible using the scripts provided in the repository (see `extract_upstream_promoter_sequences.py`, `background_cluster3.py`, etc.).
* Data files reflect the structure used during final analyses for motif enrichment, GRN construction, and network perturbation.

---

For further details, see the [README.md](../README.md) or the full thesis report.
