library(ArchR)
library(parallel)

addArchRThreads(threads = 16)
addArchRGenome("hg19")
set.seed(1)

# sample info of 2,971 cells from GSE74310 & GSE96769:
sample.info <- read.table(
  file = "../data/sample-info.tsv", sep = "\t",
  header = T, check.names = F, row.names = 1
)

# filter low-quality ATAC-seq fragments, low-quality cells, and doublets:
arrow.files <- createArrowFiles(
  inputFiles = c(buenrostro2018 = "aligned/BC_MQ30_posSorted_merged.BAM"),
  minTSS = 4, minFrags = 5000, maxFrags = 2 * 10^5,
  excludeChr = c("chrX", "chrY", "chrM"),
  bcTag = "CB", validBarcodes = sample.info$Run
)
doubScores <- addDoubletScores(arrow.files)
archr.obj <- filterDoublets(ArchRProject(arrow.files))

# 2,154 cells retained for all analyses:
sample.info$discard <- !(sample.info$Run %in% gsub("^.+#", "", archr.obj$cellNames))
write.table(sample.info, file = "../data/sample-info.tsv", sep = "\t", row.names = T)

