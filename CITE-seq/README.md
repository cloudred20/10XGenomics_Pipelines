### Five major steps in CITE-seq analysis include:

#### (1) Preprocessing & QC.
* Filter low-quality cells based on RNA counts, mitochondrial percentage, and ADT non-specific binding.
  
#### (2) ADT Normalization
* CLR (Centered Log-Ratio) normalization or DSB (Denoised and Scaled by Background) normalization for surface proteins.
  
#### (3) RNA Normalization & Clustering.
* SCTransform followed by dimensionality reduction (PCA, UMAP).
  
#### (4) TCR V(D)J Integration.
* Parse Cell Ranger `filtered_contig_annotations.csv` files and attach clonotype definitions to the single-cell expression object.
  
#### (5) Downstream Analysis.
* Differential expression, clonal diversity metrics, and cluster-specific TCR expansion visualization
