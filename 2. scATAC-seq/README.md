## scATAC-seq Data Analysis Pipeline
This repository contains an end-to-end bioinformatics pipeline for processing, analyzing, and interpreting single-cell chromatin accessibility data (scATAC-seq). The workflow guides you from raw fragment files to peak calling, cell clustering, motif enrichment, and transcription factor (TF) activity inference.

### Four major steps in scATAC-seq analysis include:
#### (1) Quality Control.
* Compute TSS enrichment scores and nucleosome signal metrics to filter out dead cells, low-quality barcodes, and doublets.
  
#### (2) Normalization & Dimensionality Reduction.
* Run TF-IDF (Term Frequency-Inverse Document Frequency) normalization followed by Singular Value Decomposition (SVD) / LSI.
  
#### (3) Clustering & Embedding.
* Build a Shared Nearest Neighbor (SNN) graph and visualize cellular landscapes via UMAP.
  
#### (4) Differential Accessibility.
* Identify cell-type-specific peaks, open chromatin regions, and cis-regulatory elements.
  
#### (5) Motif Enrichment & TF Activity.
* Compute per-cell motif activity scores using ChromVAR to identify driver transcription factors.

### Resources
