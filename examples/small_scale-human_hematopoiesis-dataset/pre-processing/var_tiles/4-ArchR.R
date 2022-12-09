setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

library(ArchR)
library(GenomicRanges)
library(Biostrings)
library(Matrix)
library(dplyr)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg19)

addArchRThreads(threads = 16)
addArchRGenome("hg19")
set.seed(1)

# sample info of 2,971 cells from GSE74310 & GSE96769:
sample.info <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
)

# pre-processing:
arrow.files <- createArrowFiles(
  inputFiles = c(small_hematopoietic = "pre-processing/var_tiles/aligned/BC_MQ30_posSorted_merged.BAM"),
  minTSS = 4, minFrags = 5000, maxFrags = 2 * 10^5,
  excludeChr = c("chrX", "chrY", "chrM"),
  bcTag = "CB", validBarcodes = sample.info$Run,
  addGeneScoreMat = F
)
doubScores <- addDoubletScores(arrow.files)
archr.obj <- ArchRProject(arrow.files) %>% filterDoublets()

archr.obj$Run <- gsub("^.+#", "", archr.obj$cellNames)
write.csv(
  x = archr.obj@cellColData[, c(16, 2:15)], row.names = F,
  file = "ArchR-results/cellColData.csv"
)

# 2,154 cells retained for all analyses:
sample.info$discard <- !(sample.info$Run %in% archr.obj$Run)
write.table(
  x = sample.info, sep = "\t", row.names = F,
  file = "sample-info.tsv"
)

# iterative LSI (variable tiles):
archr.obj <- addIterativeLSI(
  ArchRProj = archr.obj,
  iterations = 5,
  varFeatures = 50000,
  dimsToUse = 1:30
)
saveArchRProject(archr.obj)

# ArchR embedding for benchmarking:
LSI.scaled <- getReducedDims(archr.obj)
write.csv(
  x = LSI.scaled, file = "embeddings/itLSI_ArchR-var_tiles.csv",
  row.names = gsub("^.+#", "", rownames(LSI.scaled))
)

# UMAP and clustering:
clusters <- Seurat::FindNeighbors(
  object = LSI.scaled,
  distance.matrix = F,
  nn.method = "annoy",
  annoy.metric = "euclidean",
  k.param = 20
)$snn %>% Seurat::FindClusters(resolution = 0.8)
umap <- Seurat::RunUMAP(
  object = LSI.scaled,
  metric = "euclidean",
  n.neighbors = 20,
  min.dist = 0.4,
  spread = 0.9
)@cell.embeddings %>% data.frame()
write.csv(
  x = cbind(umap, Cluster = as.integer(clusters$res.0.8)),
  file = "ArchR-results/UMAP_and_Clusters.csv",
  row.names = gsub("^.+#", "", rownames(LSI.scaled))
)

# extract cell by variable tile accessibility matrix:
get.tile.matrix <- function(archr.obj, genome){
  tile.mtx <- getMatrixFromProject(archr.obj, useMatrix = "TileMatrix", binarize = T)
  var.tiles <- archr.obj@reducedDims$IterativeLSI$LSIFeatures[, -3]
  var.tiles.gr <- GRanges(
    seqinfo = seqinfo(genome),
    seqnames = var.tiles$seqnames, strand = "+",
    ranges = IRanges(
      start = var.tiles$start,
      width = archr.obj@reducedDims$IterativeLSI$tileSize,
      names = paste0("tile", var.tiles$idx)
    )
  )
  return(list(
    var.tiles = var.tiles, var.tiles.gr = var.tiles.gr, genome = genome,
    var.tile.mtx = assays(tile.mtx)$TileMatrix[match(var.tiles, tile.mtx@elementMetadata), ]
  ))
}
archr.res <- get.tile.matrix(archr.obj, genome = BSgenome.Hsapiens.UCSC.hg19)

# input for CellSpace (CellSpace-train-var_tiles.sh):
sample.info <- subset(sample.info, !discard)
ci <- match(sample.info$Run, archr.obj$Run)
writeMM(
  obj = t(archr.res$var.tile.mtx[, ci]),
  file = "CellSpace-inputs/GSE96769-cell_by_tile.mtx"
)
archr.res$var.tiles.gr %>% getSeq(x = archr.res$genome) %>%
  writeXStringSet(filepath = "CellSpace-inputs/GSE96769-top50K_variable_tiles.fa")

