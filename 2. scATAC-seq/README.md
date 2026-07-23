## scATAC-seq Data Analysis Pipeline
This repository contains an end-to-end bioinformatics pipeline for processing, analyzing, and interpreting single-cell chromatin accessibility data (scATAC-seq). The workflow guides you from raw fragment files to peak calling, cell clustering, motif enrichment, and transcription factor (TF) activity inference.

### Key workflow steps in scATAC-seq analysis include:
#### 1. Quality Control & Cell Filtering
* Compute single-cell quality metrics—including TSS enrichment scores, nucleosome signal ratios, fragment counts in peaks, and genomic blacklist ratios—to filter out damaged cells, low-depth barcodes, and potential doublets.

#### 2. Normalization & Dimensionality Reduction
* Normalize cell-by-peak matrices using Term Frequency-Inverse Document Frequency (TF-IDF) normalization, followed by Singular Value Decomposition (SVD) to generate Latent Semantic Indexing (LSI) low-dimensional embeddings.

#### 3. Clustering & Non-Linear Visualization
* Construct Shared Nearest Neighbor (SNN) graphs on selected biological LSI dimensions and project cellular landscapes using UMAP to identify distinct regulatory cell states.

#### 4. Gene Activity Quantification
* Convert sparse chromatin accessibility (peaks) into pseudo-gene expression levels based on accessibility within gene bodies and promoters, enabling initial cell-type identification using known RNA markers.

#### 5. Cell Type Annotation & Label Transfer
* Map scATAC-seq clusters to reference single-cell RNA-seq (scRNA-seq) datasets via anchor transfer, or annotate manually using canonical cell-type markers (e.g., CD3D for T cells, MS4A1 for B cells, CD14 for Monocytes).

#### 6. Differential Accessibility Analysis & Track Visualization
* Perform marker peak detection to identify cell-type-specific open chromatin regions, promoter accessibility profiles, and candidate cis-regulatory elements (CREs) that are significantly differentially accessible across clusters. Inspect accessibility profiles and generate locus-specific track plots at key genes, promoters, and distal enhancers across distinct cell types.

#### 7. Motif & Transcription Factor Enrichment
* Integrate transcription factor binding databases (e.g., JASPAR) to identify overrepresented TF binding motifs within differentially accessible peaks and uncover regulatory drivers.

#### 8. ChromVAR Footprinting & Single-Cell TF Activity
* Infer transcription factor activity on a per-cell basis using chromVAR z-scores and evaluate fine-scale TF binding footprints at target genomic sites.
---
### Repository Architecture
#### scATAC-seq - single-cell chromatin accessibility pipeline
* README.md        
* **standard_scATACseq_analysis**: This pipeline processes sparse, high-dimensional single-cell chromatin accessibility data from raw fragment files into biologically interpretable cell clusters:
  * **Workflow**: quality control filtering, normalization & scaling, linear dimensionality reduction, graph construction and clustering.
  * **Input data**: 10k_pbmc_ATACv2_nextgem_Chromium_Controller (filtered_peak_bc_matrix.h5, fragments.tsv, singlecell.csv)
  * **Analysis pipeline**: scATAC-seq_analysis.md
  * **Figures**: scATAC-seq_analysis_files
---
### Resources
1. https://www.youtube.com/watch?v=yEKZJVjc5DY&list=PLJefJsd1yfhagnkss5B1YCsHaH0GWQfFT&index=14
2. https://stuartlab.org/signac/articles/pbmc_vignette#non-linear-dimension-reduction-and-clustering
