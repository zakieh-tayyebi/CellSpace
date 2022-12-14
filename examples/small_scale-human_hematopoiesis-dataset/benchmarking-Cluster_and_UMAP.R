setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

library(dplyr)

sample.info <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
) %>% subset(!discard & Cell_type != "unknown")
num.clusters <- length(unique(sample.info$Cell_type))

k <- 20
resolution <- seq(0.1, 3, 0.05)
min.dist <- c(
  `LSI-all_peaks.csv` = 0.3,
  `CellSpace-all_peaks.csv` = 0.1,
  `itLSI_ArchR-var_tiles.csv` = 0.4,
  `CellSpace-var_tiles.csv` = 0.2,
  `chromVAR-motifs.csv` = 0.1,
  `chromVAR-kmers.csv` = 0.1,
  `SIMBA-peaks.csv` = 0.4,
  `SIMBA-peaks_kmers_motifs.csv` = 0.3,
  `scBasset.csv` = 0.2,
  `PeakVI.csv` = 0.4,
  `scBasset-batch_corrected.csv` = 0.1,
  `PeakVI-batch_corrected.csv` = 0.1,
  `SIMBA-batch_corrected.csv` = 0.3
)
spread <- c(
  `LSI-all_peaks.csv` = 0.9,
  `CellSpace-all_peaks.csv` = 1,
  `itLSI_ArchR-var_tiles.csv` = 0.9,
  `CellSpace-var_tiles.csv` = 1,
  `chromVAR-motifs.csv` = 1,
  `chromVAR-kmers.csv` = 1,
  `SIMBA-peaks.csv` = 0.8,
  `SIMBA-peaks_kmers_motifs.csv` = 0.9,
  `scBasset.csv` = 0.5,
  `PeakVI.csv` = 1,
  `scBasset-batch_corrected.csv` = 1,
  `PeakVI-batch_corrected.csv` = 1,
  `SIMBA-batch_corrected.csv` = 1
)

emb.fn <- list.files(path = "embeddings/", pattern = "csv")
clusters <- lapply(emb.fn, function(fn){
  metric <- ifelse(grepl("CellSpace", fn), "cosine", "euclidean")
  emb <- read.csv(paste0("embeddings/", fn), header = T, row.names = 1)
  emb <- as.matrix(emb[sample.info$Run, ])

  umap <- Seurat::RunUMAP(
    object = emb,
    metric = metric,
    n.neighbors = k,
    min.dist = min.dist[fn],
    spread = spread[fn],
  )@cell.embeddings %>% data.frame()
  write.csv(umap, file = paste0("benchmarking/UMAP/", fn))

  cl <- Seurat::FindNeighbors(
    object = emb,
    distance.matrix = F,
    nn.method = "annoy",
    annoy.metric = metric,
    k.param = k
  )$snn %>% Seurat::FindClusters(resolution = resolution)
  cl.num <- sapply(1:ncol(cl), function(i){ length(levels(cl[, i])) })
  if(any(cl.num == num.clusters))
    return(cl[, which(cl.num == num.clusters)[1]])
}) %>% do.call(what = cbind)
colnames(clusters) <- gsub(".csv", "", emb.fn)
write.csv(clusters, "benchmarking/Clusters.csv")

