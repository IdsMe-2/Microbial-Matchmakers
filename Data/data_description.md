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

## 2. `Lj_expression_data.xlsx`

**Description:**

* Same structure as `Expression_data_At.xlsx`, but for *Lotus japonicus* gene expression data.
* Also contains 16 samples under the 4 conditions described above.
* Used for cross-species comparison and identification of Arabidopsis homologs for Lotus clusters.

---

## 3. `cluster_table_KW.xlsx`

**Description:**

* Cross-species cluster table linking *Arabidopsis* genes and their *Lotus* homologs.
* Includes cluster assignments and GO term annotations.

**Columns:**

* `ID`: Arabidopsis gene ID (with transcript version)
* `cl`: Arabidopsis cluster number
* `ljID`: corresponding Lotus homolog ID
* `ljcl1`, `ljcl2`: cluster assignments for the Lotus gene in two different clustering schemes
* `GOID`: semicolon-separated Gene Ontology IDs
* `Term`: associated biological processes
* `TF_activity_AT`: whether the Arabidopsis gene is a transcription factor

---

## 4. `selected_upstream_sequences_cluster3.fasta`

**Description:**

* Contains 1000 bp upstream promoter sequences for genes in Cluster 3 of *Arabidopsis*.
* FASTA headers correspond to TAIR10 gene IDs.
* These sequences were scanned with known TF motifs using **FIMO**.

---

## 5. `selected_upstream_sequences_lotus_cluster6_Background_15_5_2025.fasta`

**Description:**

* Promoter sequences (Arabidopsis homologs) corresponding to *Lotus* Cluster 6 background genes.
* Used for comparative motif analysis between species.

---

## 6. `ALL_plant_motifs_JASPAR.meme`

**Description:**

* MEME-format motif database for plant transcription factors.
* Downloaded from the JASPAR 2022 database.
* Used in all FIMO-based motif scans in the project.

---

## 7. `TAIR10_GFF3_genes.zip`

**Description:**

* Compressed file containing the GFF3 genome annotation for *Arabidopsis thaliana*.
* Used for mapping FIMO hits to gene features in `Gene_annotation.py`.

---

## 8. `TAIR10_upstream_1000_20101104.zip`

**Description:**

* Compressed file of Arabidopsis promoter sequences (1 kb upstream of TSS).
* Source: TAIR10 upstream annotations.
* Used in `extract_upstream_promoter_sequences.py` to extract relevant sequences.

---

For further details, see the [README.md](../README.md) or the full thesis report.
