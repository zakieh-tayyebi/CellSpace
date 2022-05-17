library(ArchR)
library(dplyr)

addArchRThreads(threads = 8)
addArchRGenome("hg19")
set.seed(1)

arrow.files <- createArrowFiles(
  inputFiles = "aligned/BC_MQ30_posSorted_merged.BAM", sampleNames = "GSE96769",
  minFrags = 10 ^ 3.5, addGeneScoreMat = F, bcTag = "CB"
)
doubScores <- addDoubletScores(arrow.files)
archr.obj <- ArchRProject(arrow.files) %>% filterDoublets()

nn.k <- 40
archr.obj <- addIterativeLSI(
  archr.obj, iterations = 10, varFeatures = 50000,
  dimsToUse = 1:30, UMAPParams = list(n_neighbors = nn.k)
)
archr.obj <- addClusters(archr.obj, k.param = nn.k, resolution = c(0.4, 0.8, 1.0, 1.5)) %>%
  addUMAP(nNeighbors = nn.k)

# SRA meta-data table from GEO
SRA.tb <- read.table("SraRunTable.txt", sep = ",", header = T, stringsAsFactors = F, check.names = F)
archr.obj$Celltype <- SRA.tb$cell_type[match(archr.obj$cellNames, paste0("GSE96769#", SRA.tb$Run))]

saveArchRProject(archr.obj)
write.csv(archr.obj@cellColData, file = "ArchROutput/cellColData.csv", quote = F)
write.csv(getReducedDims(archr.obj), file = "ArchROutput/IterativeLSI/LSI.csv")
write.csv(getEmbedding(archr.obj), row.names = T, quote = F, file = "ArchROutput/Embeddings/UMAP.csv")
