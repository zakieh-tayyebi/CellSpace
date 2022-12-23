setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

source("pre-processing/all_peaks/LSI-functions.R") # ArchR functions
library(GenomicRanges)
library(Biostrings)
library(Matrix)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19

# downloaded from GEO (GSE96769):
#   GSE96769_PeakFile_20160207.bed.gz
#   GSE96769_scATACseq_counts.txt.gz

peaks <- fread(
  file = "pre-processing/all_peaks/GSE96769_PeakFile_20160207.bed.gz",
  header = F, stringsAsFactors = F
)[, c(1:3, 7)]
colnames(peaks) <- c("chr", "start", "end", "annot")
peaks.gr <- makeGRangesFromDataFrame(
  df = peaks,
  seqinfo = seqinfo(genome),
  starts.in.df.are.0based = T
)

cellname <- strsplit(
  x = readLines("pre-processing/all_peaks/GSE96769_scATACseq_counts.txt.gz", n = 1),
  split = ";|(#\t)"
)[[1]][-1]
counts <- read.table(
  file = "pre-processing/all_peaks/GSE96769_scATACseq_counts.txt.gz",
  skip = 1, header = F, sep = "\t",
  col.names = c("peak", "cell", "count")
)
counts <- sparseMatrix(
  i = counts$cell, j = counts$peak, x = counts$count,
  dimnames = list(cellname, as.character(peaks.gr))
) %>% t()

# cells filtered by ArchR (pre-processing/var_tiles/5-ArchR.R):
sample.info <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
) %>% subset(!discard)
ci <- match(sample.info$Title, colnames(counts))
counts <- counts[, ci]

# filter peaks
pi <- !(peaks$chr %in% c("chrX", "chrY", "chrM")) &
  !grepl("promoter", peaks$annot) &
  (rowSums(counts > 0) >= 5)
peaks.gr <- peaks.gr[pi]
counts <- counts[pi, ]

saveRDS(peaks.gr, file = "pre-processing/all_peaks/GSE96769-filtered_peaks.rds")
saveRDS(counts, file = "pre-processing/all_peaks/GSE96769-filtered_counts.rds")

# compute LSI for filtered peak-by-cell count matrix:
LSI <- computeLSI(
  mat = counts,
  LSIMethod = 2,
  scaleTo = 10000,
  nDimensions = 30,
  outlierQuantiles = c(0.02, 0.98)
)

# LSI embedding for benchmarking:
write.csv(
  x = scaleDims(LSI$matSVD),
  row.names = sample.info$Run[match(colnames(counts), sample.info$Title)],
  file = "embeddings/LSI-all_peaks.csv"
)

# input for CellSpace (CellSpace-train-all_peaks.sh):
getSeq(peaks.gr, x = genome) %>%
  writeXStringSet(filepath = "CellSpace-inputs/GSE96769-filtered_peaks.fa")
writeMM(t(counts), file = "CellSpace-inputs/GSE96769-cell_by_peak.mtx")

