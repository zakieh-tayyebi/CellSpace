# Small-scale human hematopoiesis dataset

## Pre-processing scATAC-seq data

### [Variable tiles](pre-processing/var_tiles/)

1.  Downloading FASTQ files,
2.  Trimming adapter sequences,
3.  Aligning reads to the reference genome,
4.  Filtering and sorting the alignments,
5.  Filtering fragments and cells with ArchR (quality-control and doublet removal), performing iterative LSI (as the ArchR itLSI embedding), and preparing CellSpace inputs (var. tiles).

### [All peaks](pre-processing/all_peaks/)

1.  Downloading peak-by-cell count matrix from GEO,
2.  Filtering cells (filtered by ArchR previously) and peaks,
3.  Computing LSI embedding (all peaks), and preparing CellSpace inputs (all peaks).

## Training a CellSpace model

-   [**Variable tiles**](CellSpace-train-var_tiles.sh)
-   [**All peaks**](CellSpace-train-all_peaks.sh)

## Downstream analyses

### ArchR itLSI (var. tiles)

-   [Clustering cells and computing UMAP embedding.](pre-processing/var_tiles/5-ArchR.R)

### CellSpace (var. tiles)

-   [Clustering cells, embedding transcription factor (TF) motifs, computing TF activity scores, and computing UMAP embedding for cells and important TFs.](CellSpace-downstream.R)
-   [Pseudotime analysis by Palantir.](palantir.py)
-   [De novo motif discovery, computing UMAP embedding for cells and de novo motifs, and finding similar known motifs.](denovo-motifs.R)

### Benchmarking

-   The results of all methods are available in ['embeddings/'](embeddings/). This repository includes the scripts used to perform CellSpace (var. tiles), CellSpace (all peaks), ArchR itLSI (var. tiles), and LSI (all peaks). All other methods (scBasset, SIMBA, PeakVI, and chromVAR) were performed, as described in [our paper](https://doi.org/10.1101/2022.05.02.490310), on the filtered peak-by-cell count matrix that was used for the 'all peaks' version.
-   [Each embedding was used to compute a UMAP embedding and cluster the cells, excluding the 'unknown' cell type.](benchmarking-Cluster_and_UMAP.R)
-   [The embedding and clustering results from each method were evaluated by biological conservation and batch correction metrics.](benchmarking-compute_metrics.py)

## Visualizing the results

All the figures related to this dataset in [our paper](https://doi.org/10.1101/2022.05.02.490310) were [visualized](plot.R) with the data available in this repository. These figures are available in ['plots/'](plots/).
