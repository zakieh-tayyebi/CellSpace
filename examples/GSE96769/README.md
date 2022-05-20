# Small-scale human hematopoiesis dataset ([GSE96769](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96769))

## [Pre-processing](pre-processing)

Paired-end FASTQ files for 2,779 samples were downloaded using [SRA toolkit v2.11](pre-processing/0-fastq_dump.sh). Adapters were detected and trimmed by [TrimGalore v0.6.6](pre-processing/1-TrimGalore.sh) (FastQC and Cutadapt wrapper). Reads were aligned to hg19 by [Bowtie2 v2.4.2](pre-processing/2-bowtie2.sh) with "--maxins 2000 --very-sensitive --no-unal --no-mixed”. Barcodes were added to the aligned BAM files as ‘CB’ tags using [a custom script](pre-processing/addBarcodeTag.cpp). BAM files were filtered (MAPQ>30), coordinate-sorted, indexed, and merged by [Samtools v1.11](pre-processing/3-samtools.sh).

The final BAM file was inputted to [ArchR v1.0.1](pre-processing/4-ArchR.R), barcodes were filtered (TSS enrichment score > 4 and 3.5 < log10(number of fragments) < 5), and doublets were detected and removed with default settings. 1,849 cells were retained for downstream analyses. Using ArchR’s implementation of iterative LSI, the top 50K most variable tiles (genome-wide 500 bp bins) were identified after 10 iterations. This was restricted to the standard chromosomes chr1, …, chr22, chrX. The first 30 LSI components were used to build a NN graph and its SNN graph and to cluster the cells (method="Seurat") to identify ArchR clusters. The ArchR UMAP embedding was computed from the NN graph.

The selected tiles and their corresponding cell-by-tile counts matrix [were extracted](pre-processing/5-CellSpace_inputs.R) for training a CellSpace model.

## [Training a CellSpace model](CellSpace-train.sh)

A CellSpace model was trained with “--k 8 --sampleLen 150 --dim 30 --ngrams 3 --exmpPerPeak 20 --epoch 20” on the cell-by-tile count matrix and genomic sequences extracted from the corresponding tiles. The resulting CellSpace embedding for cells and k-mers (k = 8) was used for downstream analyses.

## [CellSpace visualization, clustering, and motif embeddings](CellSpace-downstream.R)

