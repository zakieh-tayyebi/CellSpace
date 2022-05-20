library(CellSpace)
library(chromVARmotifs)
library(Biostrings)
library(dplyr)

cell.md <- read.csv("ArchROutput/cellColData.csv", header = T, row.names = 1)
cell.md$ArchR.cluster <- cell.md$Clusters.res.0.8

cso <- CellSpace(
  project = "GSE96769",
  emb.file = "CellSpace-results/GSE96769-embedding.tsv", # cell and k-mer embeddings
  meta.data = cell.md[, c(16, 21)], # FAC-sorted celltypes and ArchR clusters
  cell.names = rownames(cell.md)
)
cso <- find_neighbors(cso, n.neighbors = 40) # NN and SNN graphs for cells
cso <- find_clusters(cso, resolution = c(1, 1.5, 2)) # cell clusters
write.csv(cso@meta.data, "CellSpace-results/cell-clusters.csv")

cso <- run_UMAP(cso, n.neighbors = 40, min.dist = 0.1, spread = 1) # cell UMAP
write.csv(cso@reductions$cells_UMAP, file = "CellSpace-results/UMAP-cells.csv")

data(human_pwms_v1)
hg.motif.emb <- lapply(human_pwms_v1, function(motif.pwm){
  motif_embedding(cso, PWM = motif.pwm) # TF motif embedding
}) %>% do.call(what = rbind) %>% na.omit()
TF.names <- gsub("(^.+LINE\\d+_)|(_.+$)", "", rownames(hg.motif.emb))

mi <- gsub("(^ENSG\\d+_)|(_\\w(_.+)*$)", "", rownames(hg.motif.emb)) %in%
  c("LINE3170_PAX5", "LINE251_EBF1", "LINE216_TCF12", "LINE202_LYL1",
    "LINE2202_MEIS3", "LINE2729_IRF1", "LINE611_PRDM1", "LINE498_CEBPB",
    "LINE3073_ESRRA", "LINE534_MAFF", "LINE397_ATF4", "LINE558_CEBPD",
    "LINE2186_HOXA9", "LINE16121_HOXA4", "LINE2173_ALX4", "LINE2081_GATA1",
    "LINE2094_GATA3", "LINE151_ID3", "LINE430_CEBPG", "LINE568_CEBPA",
    "LINE848_BCL11A", "LINE303_TCF4") # selected TF motifs
cso <- run_UMAP(
  object = cso, name = "cells_and_TFs_UMAP",
  emb = rbind(cso@cell.emb, hg.motif.emb[mi, ]),
  n.neighbors = 40, min.dist = 0.1, spread = 1
) # CellSpace cell and TF UMAP
write.csv(cso@reductions$cells_and_TFs_UMAP, file = "CellSpace-results/UMAP-cells_and_TFs.csv")

TF.score <- cosine_similarity(x = cso@cell.emb, y = hg.motif.emb[mi, ]) %>% # cell-TF similarity
  scale() %>% t() # z-score
write.csv(TF.score, file = "CellSpace-results/TF-activity-scores.csv")
