# #####
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

plot.groups <- function(
  vis, groups, group.name, add.labels = T, pal = NULL,
  point.size = 3, label.size = 3, nx = NULL, ny = NULL
){
  gp <- ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = groups), vis, size = point.size, alpha = 0.9) +
    labs(color = group.name) + theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  if(add.labels){
    G <- sort(unique(groups))
    vis.gp <- data.frame(matrix(NA, nrow = length(G), ncol = ncol(vis),
                                dimnames = list(G, colnames(vis))))
    for(g in G)
      for(col in colnames(vis))
        vis.gp[g, col] <- median(vis[groups == g, col])
    gp <- gp + geom_label(aes(color = G, label = G), vis.gp, show.legend = F,
                          size = label.size, nudge_x = nx, nudge_y = ny)
  }
  gp + scale_color_manual(values = pal)
}

plot.continuous.value <- function(vis, val, val.name, scale.color, point.size = 3){
  ggplot(vis, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = val), size = point.size, alpha = 0.9) +
    scale.color + labs(color = val.name) +
    theme_classic() + theme(axis.ticks = element_blank(), axis.text = element_blank())
}


# #####
load("plots/palette.rda")

cell.md <- read.csv(
  file = "CellSpace-results/cell-metadata.csv",
  row.names = 1, header = T
)
cell.md$cell_type <- gsub("^GMP.+$", "GMP", cell.md$cell_type)

ArchR.UMAP <- read.csv( # ArchR UMAP of cells
  file = "ArchROutput/Embeddings/UMAP.csv",
  row.names = 1, header = T, col.names = c("cell", "UMAP_1", "UMAP_2")
)

cell.tf.UMAP <- read.csv( # CellSpace UMAP of cells and TFs in the same space
  file = "CellSpace-results/UMAP-cells_and_TFs.csv",
  row.names = 1, header = T
)
cell.only.UMAP <- read.csv( # CellSpace UMAP of cells (without TFs)
  file = "CellSpace-results/UMAP-cells.csv",
  row.names = 1, header = T
)

TF.score <- read.csv(
  file = "CellSpace-results/TF-activity-scores.csv",
  row.names = 1, header = T
)

pseudotime <- read.csv(
  file = "CellSpace-results/palantir-results/pseudotime.csv",
  header = T, row.names = 1
)
branch.probs <- read.csv(
  file = "CellSpace-results/palantir-results/branch_probs.csv",
  row.names = 1, header = T, check.names = F
)


# #####
motif.names <- rownames(cell.tf.UMAP)[grep("ENSG", rownames(cell.tf.UMAP))]
TF.names <- gsub("(^.+LINE\\d+_)|(_.+$)", "", motif.names)
nx <- c(-0.4, 0.4, -0.8, 0.8)[c(1, 2, 2, 4, 2, 2, 3, 3, 2, 2, 2, 2,
                                2, 1, 2, 3, 2, 2, 1, 1, 1, 1)]
ny <- c(-0.4, 0.4, -0.8, 0.8)[c(2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 1,
                                1, 3, 2, 2, 1, 1, 2, 3, 1, 1)]
ggsave(
  filename = "plots/Fig2a.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = cell.tf.UMAP[rownames(cell.md), ], add.labels = F,
      groups = cell.md$cell_type, group.name = "cell type", pal = ct.pal
    ) + labs(title = "CellSpace UMAP", subtitle = "UMAP computed for cells and TFs") +
    geom_point(data = cell.tf.UMAP[motif.names, ], color = "black", shape = 17, size = 3) +
    geom_label_repel(aes(label = TF.names), cell.tf.UMAP[motif.names, ], size = 7, fontface = "bold",
                     min.segment.length = 0, max.overlaps = Inf,
                     segment.color = "black", segment.size = unit(0.2, "mm"),
                     nudge_x = nx, nudge_y = ny)
)

ggsave(
  filename = "plots/Fig2c.pdf", width = 12, height = 10,
  plot =
    plot.continuous.value(
      vis = cell.tf.UMAP[rownames(cell.md), ],
      val = pseudotime$X0, val.name = "Palantir\npseudotime",
      scale.color = scale_color_gradient2(low = "#fcff9c", high = "#6e0000",
                                          mid = "#ff4900", midpoint = 0.5)) +
    labs(title = "CellSpace UMAP", subtitle = "UMAP computed for cells and TFs")
)

ggsave(
  filename = "plots/Fig2d.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = cell.tf.UMAP[rownames(cell.md), ], label.size = 5,
      groups = factor(cell.md$Clusters.res_1.5), group.name = "CellSpace\ncluster", pal = cl.pal[1:8]
    ) + labs(title = "CellSpace UMAP", subtitle = "UMAP computed for cells and TFs")
)

ggsave(
  filename = "plots/Fig2e.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = ArchR.UMAP, label.size = 5,
      groups = cell.md$cell_type, group.name = "cell type", pal = ct.pal
    ) + ggtitle("ArchR UMAP")
)

ggsave(
  filename = "plots/Fig2f.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = ArchR.UMAP, label.size = 5,
      groups = factor(cell.md$ArchR.cluster - 1), group.name = "ArchR\ncluster", pal = cl.pal2[1:10]
    ) + ggtitle("ArchR UMAP")
)

plot_grid(
  plotlist = lapply(1:ncol(branch.probs), function(i){
    end.point <- colnames(branch.probs)[i]
    nx <- c(1, 0, -0.7, 1, 0)[i]
    ny <- c(1, -1, -1, 0, 1)[i]
    start <- "GSE96769#SRR5353482"
    plot.continuous.value(
      vis = cell.tf.UMAP[rownames(cell.md), ],
      val = branch.probs[, end.point], val.name = "branch\nprobability",
      scale.color = scale_color_distiller(palette = "Spectral")) +
      geom_point(data = cell.tf.UMAP[start, ], size = 3, shape = 8) +
      geom_point(data = cell.tf.UMAP[end.point, ], size = 3, shape = 8) +
      geom_text_repel(aes(label = "start"), cell.tf.UMAP[start, ],
                      min.segment.length = 0, segment.linetype = 2,
                      nudge_x = -0.5, nudge_y = -0.5) +
      geom_text_repel(aes(label = "terminal\nstate"), cell.tf.UMAP[end.point, ],
                      min.segment.length = 0, segment.linetype = 2,
                      nudge_x = nx, nudge_y = ny)
  }), ncol = 2
) %>% ggsave(filename = "plots/Supp-Fig1a.pdf", width = 12 * 2, height = 10 * 3)

ggsave(
  filename = "plots/Supp-Fig1b.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = cell.only.UMAP, label.size = 5,
      groups = cell.md$cell_type, group.name = "cell type", pal = ct.pal
    ) + labs(title = "CellSpace UMAP", subtitle = "UMAP computed for cells only")
)

ggsave(
  filename = "plots/Supp-Fig1c.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = cell.tf.UMAP[rownames(cell.md), ], add.labels = F,
      groups = cell.md$donor, group.name = "donor", pal = dn.pal
    ) + labs(title = "CellSpace UMAP", subtitle = "UMAP computed for cells and TFs")
)

ggsave(
  filename = "plots/Supp-Fig1d.pdf", width = 12, height = 10,
  plot =
    plot.groups(
      vis = ArchR.UMAP, add.labels = F,
      groups = cell.md$donor, group.name = "donor", pal = dn.pal
    ) + ggtitle("ArchR UMAP")
)

ca <- columnAnnotation(
  df = data.frame(
    `CellSpace\ncluster` = factor(cell.md$Clusters.res_1.5),
    `cell type` = factor(cell.md$cell_type, levels = names(ct.pal)),
    check.names = F
  ),
  col = list(`CellSpace\ncluster` = cl.pal[1:8], `cell type` = ct.pal)
)
gp <- gpar(fontsize = 10)
TF.names <- gsub("(^.+LINE\\d+_)|(_.+$)", "", rownames(TF.score))
pdf("plots/Fig3e.pdf", width = 10, height = 5)
Heatmap(
  as.matrix(TF.score), name = "Similarity\nz-score",
  col = colorRamp2(c(-1, 0, 1), c("#007e8c", "white", "#a3005c")),
  cluster_columns = T, cluster_rows = T,
  top_annotation = ca, column_split = cell.md$cell_type, column_title_rot = 90, show_column_names = F,
  row_title = "", row_labels = TF.names,
  row_names_gp = gp, row_names_max_width = max_text_width(TF.names, gp = gp),
) %>% draw(heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 2, 10, 2), "mm"))
dev.off()

