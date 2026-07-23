### 10x Genomics Multi-Omic & Spatial Pipelines Repository
Welcome to the central repository for processing, analyzing, and integrating single-cell and spatial omics datasets. This repository hosts standardized, reproducible, and modular computational workflows designed for 10x Genomics technologies.

Different assays require different Cell Ranger pipelines before downstream analysis in R.

| Assay                                                                | Cell Ranger Pipeline    | Primary Output                                     | Typical R Package       |
| -------------------------------------------------------------------- | ----------------------- | -------------------------------------------------- | ----------------------- |
| Single-cell RNA-seq (Gene Expression)                                | `cellranger count`      | Gene × Cell UMI matrix                             | Seurat                  |
| Multi-library Gene Expression (RNA + ADT, CRISPR, Cell Multiplexing) | `cellranger multi`      | Integrated Feature × Cell matrix                   | Seurat                  |
| Single-cell V(D)J (TCR/BCR)                                          | `cellranger vdj`        | Clonotypes, CDR3 sequences                         | scRepertoire, immunarch |
| Single-cell ATAC-seq                                                 | `cellranger-atac count` | Peak × Cell accessibility matrix                   | Signac                  |
| Multiome (scRNA + scATAC)                                            | `cellranger-arc count`  | Gene expression + chromatin accessibility matrices | Seurat + Signac         |
| Multiple samples                                                     | `cellranger aggr`       | Aggregated feature matrices                        | Seurat / Signac         |

---

## Typical Output Files

### Gene Expression (`cellranger count` / `cellranger multi`)

```text
filtered_feature_bc_matrix/
raw_feature_bc_matrix/
molecule_info.h5
cloupe.cloupe
web_summary.html
```

**Key outputs**

| File                         | Description                                                                                                                                 |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `filtered_feature_bc_matrix` | Filtered gene expression count matrix.                                                                                                      |
| `raw_feature_bc_matrix`      | Unfiltered count matrix containing all detected barcodes.                                                                                   |
| `molecule_info.h5`           | Molecule-level information for advanced analyses (e.g., RNA velocity, downsampling).                                                        |
| `web_summary.html`           | Sequencing, mapping, and cell-calling quality metrics.                                                                                      |
| `cloupe.cloup`               | Interactive visualization file for Loupe Browser.                                                                                           |

---

### Single-cell ATAC-seq (`cellranger-atac count`)

```text
filtered_peak_bc_matrix.h5
fragments.tsv.gz
fragments.tsv.gz.tbi
singlecell.csv
peaks.bed
web_summary.html
```

**Key outputs**

| File                         | Description                                                                                                                                 |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `filtered_peak_bc_matrix.h5` | Sparse Peak × Cell accessibility matrix used to create a Signac object.                                                                     |
| `fragments.tsv.gz`           | Fragment-level genomic coordinates for every sequenced fragment. Required for QC, coverage plots, footprinting, and gene activity analyses. |
| `fragments.tsv.gz.tbi`       | Tabix index enabling rapid genomic access to the fragment file.                                                                             |
| `singlecell.csv`             | Per-cell quality metrics including fragments, FRIP, TSS enrichment, nucleosome signal, and barcode statistics.                              |
| `peaks.bed`                  | Consensus peak coordinates identified across all cells.                                                                                     |
| `web_summary.html`           | Sequencing and ATAC-specific quality metrics.                                                                                               |

---

### Multiome (`cellranger-arc count`)

Outputs include both RNA and ATAC modalities.

```text
filtered_feature_bc_matrix/
filtered_peak_bc_matrix.h5
fragments.tsv.gz
per_barcode_metrics.csv
web_summary.html
```

These files enable integrated analysis of gene expression and chromatin accessibility using **Seurat**, **Signac**, or **ArchR**.

---

## Recommended Downstream R Packages

| Data Type  | Recommended Packages         |
| ---------- | ---------------------------- |
| scRNA-seq  | Seurat, SingleCellExperiment |
| CITE-seq   | Seurat, dsb, scuttle         |
| scATAC-seq | Signac, ArchR, chromVAR      |
| Multiome   | Seurat + Signac, ArchR       |
| TCR/BCR    | scRepertoire, immunarch      |

---

### Repository Architecture
#### (1) scRNA-seq - single-cell RNA transcriptomics pipeline
* README.md
* scripts - quality control and doublet removal, normalization, PCA dimensionality reduction, UMAP visualization, clustering, and marker gene identification, etc.
#### (2) scATAC-seq - single-cell chromatin accessibility pipeline
* README.md
* scripts - ATAC-specific QC, TF-IDF normalization, SVD/LSI dimensionality reduction, UMAP clustering, differential accessibility testing, and chromVAR motif enrichment, etc.

#### (3) CITE-seq - surface epitope + transcriptomic + TCR integration
* README.md        
* scripts - WNN and clonal tracking scripts

#### (4) Visium (Spatial Transcriptomics)           
* [Coming Soon] Standardized workflow utilizing Space Ranger and Seurat/SeuratDisk or Squidpy to map whole-transcriptome expression overlayed onto H&E histology sections.

#### (5) Xenium (In Situ Subcellular Mapping)          
* [Coming Soon] High-resolution targeted transcripts analysis pipeline optimized for cell segmentation boundaries, subcellular molecule visualization, and spatial neighborhood indexing.


