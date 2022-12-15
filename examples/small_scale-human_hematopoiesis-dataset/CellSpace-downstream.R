setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

library(CellSpace)
library(JASPAR2020)
library(chromVARmotifs)
library(TFBSTools)
library(Biostrings)
library(parallel)
library(dplyr)

sample.info <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
) %>% subset(!discard) # cells filtered by ArchR (pre-processing/var_tiles/5-ArchR.R)


# CellSpace (all peaks) #####
# read CellSpace embedding matrix (CellSpace-train-all_peaks.sh):
cso <- CellSpace(
  project = "CellSpace all_peaks",
  emb.file = "CellSpace-results/all_peaks/CellSpace-embedding.tsv",
  cell.names = sample.info$Run
)

# embedding for benchmarking:
write.csv(cso@cell.emb, "embeddings/CellSpace-all_peaks.csv")


# CellSpace (variable tiles) #####
# read CellSpace embedding matrix (CellSpace-train-var_tiles.sh):
cso <- CellSpace(
  project = "CellSpace var_tiles",
  emb.file = "CellSpace-results/var_tiles/CellSpace-embedding.tsv",
  cell.names = sample.info$Run
)

# embedding for benchmarking:
write.csv(cso@cell.emb, "embeddings/CellSpace-var_tiles.csv")

# clustering:
cso <- find_neighbors(cso, n.neighbors = 20) %>%
  find_clusters(resolution = 1.2)
write.csv(
  data.frame(
    Cluster = as.integer(cso$Clusters.res_1.2),
    row.names = sample.info$Run
  ), file = "CellSpace-results/CellSpace-Clusters.csv",
)

# CellSpace embedding and activity scores for JASPAR motifs:
jaspar <- getMatrixSet(JASPAR2020@db, opts = list(species = "Homo sapiens", collection = "CORE"))
names(jaspar) <- TFBSTools::name(jaspar)
cso <- add_motif_db(cso, PWMs = jaspar, db.name = "JASPAR2020")
write.csv(
  x = cso@motif.emb$JASPAR2020,
  file = "CellSpace-results/var_tiles/TF_motif-embedding/JASPAR2020.csv"
)

# CellSpace embedding and activity scores for chromVAR motifs:
data(human_pwms_v1)
cso <- add_motif_db(cso, PWMs = human_pwms_v1, db.name = "chromVAR_human_v1")
write.csv(
  x = cso@motif.emb$chromVAR_human_v1,
  file = "CellSpace-results/var_tiles/TF_motif-embedding/chromVARmotifs.csv"
)

# CellSpace embedding and activity scores for CisBP motifs:
cisbp <- readRDS("CisBP-Homo_sapiens_2022_09_04_8-24_pm/CisBP-PWMs.rds")
cso <- add_motif_db(cso, PWMs = cisbp, db.name = "CisBP")
write.csv(
  x = cso@motif.emb$CisBP,
  file = "CellSpace-results/var_tiles/TF_motif-embedding/CisBP.csv"
)

# UMAP of cells and TFs in the same space:
motifs <- read.csv("plots/TF-motifs.csv")
motifs <- split(motifs, motifs$db)
umap <- Seurat::RunUMAP(
  object = rbind(
    cso@cell.emb,
    cso@motif.emb$CisBP[motifs$CisBP$motif, ],
    cso@motif.emb$chromVAR_human_v1[motifs$chromVAR_human_v1$motif, ]
  ),
  metric = "cosine", n.neighbors = 50,
  min.dist = 0.2, spread = 1,
  seed.use = 1
)@cell.embeddings %>% data.frame()
write.csv(umap, file = "CellSpace-results/var_tiles/UMAP-cells_and_TFs.csv")

