### 10x Genomics Multi-Omic & Spatial Pipelines Repository
Welcome to the central repository for processing, analyzing, and integrating single-cell and spatial omics datasets. This repository hosts standardized, reproducible, and modular computational workflows designed for 10x Genomics technologies.

### Repository Architecture

```text
├── scRNA-seq/           # Single-cell RNA transcriptomics pipeline
│   ├── README.md        # Detailed scRNA instructions
│   └── scripts/         # QC, clustering, and marker scripts
├── scATAC-seq/          # Single-cell chromatin accessibility pipeline
│   │   ├── README.md        # Detailed scATAC instructions
│   │   └── scripts/         # TF-IDF, LSI, and peak calling scripts
├── CITE-seq/            # Surface Epitope + Transcriptomic integration (TCR/BCR)
│   │   ├── README.md        # Detailed CITE/VDJ instructions
│   │   └── scripts/         # WNN and clonal tracking scripts
├── Visium/              # [Coming Soon] Spatial Transcriptomics
├──  Xenium/              # [Coming Soon] In Situ Subcellular Mapping
