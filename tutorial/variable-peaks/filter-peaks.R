library(GEOquery)
library(GenomicRanges)
library(Matrix)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

# downloaded data from GEO:
getGEOSuppFiles(GEO = "GSE96769", makeDirectory = F)

# create GRanges object for peaks:
peaks <- fread(
  file = "GSE96769_PeakFile_20160207.bed.gz",
  header = F, stringsAsFactors = F
)[, c(1:3, 7)]
colnames(peaks) <- c("chr", "start", "end", "annot")
peaks.gr <- makeGRangesFromDataFrame(
  df = peaks,
  seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),
  starts.in.df.are.0based = T
)

# create dgCMatrix object for counts:
cellname <- strsplit(
  x = readLines("GSE96769_scATACseq_counts.txt.gz", n = 1),
  split = ";|(#\t)"
)[[1]][-1]
counts <- read.table(
  file = "GSE96769_scATACseq_counts.txt.gz",
  skip = 1, header = F, sep = "\t",
  col.names = c("peak", "cell", "count")
)
counts <- sparseMatrix(
  i = counts$cell, j = counts$peak, x = counts$count,
  dimnames = list(cellname, as.character(peaks.gr))
) %>% t()

# cells filtered by ArchR (../variable-tiles/5-ArchR.R):
sample.info <- read.table(
  file = "sample-info.tsv", sep = "\t",
  header = T, check.names = F, row.names = 1
) %>% subset(!discard)
ci <- match(sample.info$Title, colnames(counts))
counts <- counts[, ci]
colnames(counts) <- sample.info$Run

# filter peaks
num.cells <- rowSums(counts > 0)
pi <- !(peaks$chr %in% c("chrX", "chrY", "chrM")) &
  !grepl("promoter", peaks$annot) &
  (num.cells >= 5)
peaks.gr <- peaks.gr[pi]
counts <- counts[pi, ]

