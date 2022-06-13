library(ArchR)
library(dplyr)
library(stringr)

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

sample.info <- read.csv("sample-info.csv", header = T)[, 1:2]
SRA.tb <- read.table("SraRunTable.txt", sep = ",", header = T, stringsAsFactors = F, check.names = F)
SRA.tb$Title <- sample.info$Title[match(SRA.tb$`Sample Name`, sample.info$Accession)]
SRA.tb$donor <- str_extract(SRA.tb$Title, "BM\\d+")
SRA.tb$donor[grep("singles-MEP-141017", SRA.tb$Title)] <- "BM4983"
archr.obj@cellColData <- cbind(
  archr.obj@cellColData,
  SRA.tb[match(archr.obj$cellNames, paste0("GSE96769#", SRA.tb$Run)), c("cell_type", "donor")]
)

saveArchRProject(archr.obj)
write.csv(archr.obj@cellColData, file = "ArchROutput/cellColData.csv", quote = F)
write.csv(getReducedDims(archr.obj), file = "ArchROutput/IterativeLSI/LSI.csv")
write.csv(getEmbedding(archr.obj), row.names = T, quote = F, file = "ArchROutput/Embeddings/UMAP.csv")
