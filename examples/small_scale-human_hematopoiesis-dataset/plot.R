# #####
setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

library(ggplot2)
library(grid)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

set.seed(1)

plot.groups <- function(
    vis, mapping = aes(x = UMAP_1, y = UMAP_2),
    groups, add.labels = T, pal = NULL,
    point.size = 1.5, label.size = 5,
    et = element_text(size = 15),
    sx = 0, sy = 0, nx = NULL, ny = NULL
){
  eb <- element_blank()
  gp <- ggplot(mapping = mapping) + theme_classic() +
    geom_point(aes(color = groups), vis, size = point.size, alpha = 0.9) +
    theme(axis.ticks = eb, axis.text = eb,
          axis.title = et, title = et,
          legend.title = et, legend.text = et)
  if(add.labels){
    G <- sort(unique(groups))
    NX <- rep(sx, length(G)); names(NX) <- G
    if(!is.null(nx)) NX[names(nx)] <- nx
    NY <- rep(sy, length(G)); names(NY) <- G
    if(!is.null(ny)) NY[names(ny)] <- ny
    vis.gp <- data.frame(matrix(NA, nrow = length(G), ncol = ncol(vis),
                                dimnames = list(G, colnames(vis))))
    for(g in G)
      for(col in colnames(vis))
        vis.gp[g, col] <- median(vis[groups == g, col])
    gp <- gp + geom_label(aes(color = G, label = G), vis.gp, show.legend = F,
                          size = label.size, nudge_x = NX, nudge_y = NY)
  }
  gp +
    scale_color_manual(values = pal) +
    guides(colour = guide_legend(override.aes = list(size = 5), title.hjust = 0.5))
}

plot.continuous.value <- function(
    vis, mapping = aes(x = UMAP_1, y = UMAP_2),
    scale.color,
    val, point.size = 1.5,
    et = element_text(size = 15)
){
  eb <- element_blank()
  ggplot(mapping = mapping) +
    geom_point(aes(color = val), data = vis, size = point.size, alpha = 0.9) +
    scale.color + theme_classic() +
    theme(axis.ticks = eb, axis.text = eb,
          axis.title = et, title = et,
          legend.title = et, legend.text = et) +
    guides(colour = guide_legend(override.aes = list(size = 5), title.hjust = 0.5))
}

pal <- readRDS("plots/palette.rds")
cell.md <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
) %>% subset(!discard)
rownames(cell.md) <- cell.md$Run

emb.names <- c(
  `CellSpace-var_tiles` = "CellSpace (var. tiles)",
  `CellSpace-all_peaks` = "CellSpace (peaks)",
  `itLSI_ArchR-var_tiles` = "ArchR itLSI (var. tiles)",
  `LSI-all_peaks` = "LSI (peaks)",
  `chromVAR-motifs` = "chromVAR (motifs)",
  `chromVAR-kmers` = "chromVAR (kmers)",
  `scBasset` = "scBasset (peaks)",
  `scBasset-batch_corrected` = "scBasset (batch-corrected)",
  `PeakVI` = "PeakVI (peaks)",
  `PeakVI-batch_corrected` = "PeakVI (batch-corrected)",
  `SIMBA-peaks` = "SIMBA (peaks)",
  `SIMBA-peaks_kmers_motifs` = "SIMBA (peaks+kmers+motifs)",
  `SIMBA-batch_corrected` = "SIMBA (batch-corrected)"
)


# #####
umap <- read.csv(
  file = "CellSpace-results/var_tiles/UMAP-cells_and_TFs.csv",
  row.names = 1, header = T
)
cl <- read.csv(
  file = "CellSpace-results/var_tiles/CellSpace-Clusters.csv",
  row.names = 1, colClasses = "factor"
)
pl <- data.frame(
  pseudotime = read.csv("CellSpace-results/var_tiles/palantir/pseudotime.csv", row.names = 1)[, 1],
  branch_probs = read.csv("CellSpace-results/var_tiles/palantir/branch_probs.csv", row.names = 1)
)
motifs <- read.csv("plots/TF-motifs.csv")

nx <- rep(0, nrow(motifs)); names(nx) <- motifs$TF; ny <- nx
nx['EBF1'] <- -0.5; nx['PAX5'] <- 0.2; ny['EBF1'] <- 0.2;  ny['PAX5'] <- -0.5
ggsave(
  filename = "plots/Fig2a.pdf", width = 8, height = 5,
  plot.groups(
    vis = umap[cell.md$Run, ], groups = cell.md$Cell_type,
    pal = pal$Cell_type, add.labels = F
  ) +
  geom_point(data = umap[motifs$motif, ], color = "black", shape = 17, size = 1.5) +
  geom_label_repel(
    mapping = aes(label = motifs$TF), data = umap[motifs$motif, ],
    size = 3, fontface = "bold", min.segment.length = 0, label.padding = unit(0.5, "mm"),
    max.overlaps = Inf, segment.color = "black", segment.size = unit(0.2, "mm"),
    nudge_x = nx, nudge_y = ny
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Cell type")
)

ggsave(
  filename = "plots/Fig2c.pdf", width = 8, height = 5,
  plot.continuous.value(
    vis = umap[cell.md$Run, ], val = pl[cell.md$Run, "pseudotime"],
    scale.color = scale_color_gradient2(
      low = "#fcff9c", high = "#6e0000",
      mid = "#ff4900", midpoint = 0.5
    )
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Palantir\npseudotime")
)

plot_grid(
  plotlist = lapply(2:ncol(pl), function(i){
    end.point <- gsub("branch_probs\\.", "", colnames(pl)[i])
    nx <- c(1, -1, 0, 0, 1, 0)[i - 1]
    ny <- c(0, 0, -1, -1, 0, -1)[i - 1]
    start <- "SRR5353471"
    plot.continuous.value(
      vis = umap[cell.md$Run, ], val = pl[cell.md$Run, i],
      scale.color = scale_color_distiller(palette = "Spectral")
    ) +
    geom_point(data = umap[start, ], size = 3, shape = 8) +
    geom_point(data = umap[end.point, ], size = 3, shape = 8) +
    geom_text_repel(
      mapping = aes(label = "start"), umap[start, ],
      min.segment.length = 0, segment.linetype = 2, nudge_x = 1
    ) +
    geom_text_repel(
      mapping = aes(label = "terminal\nstate"), umap[end.point, 1:2],
      min.segment.length = 0, segment.linetype = 2,
      nudge_x = nx, nudge_y = ny
    ) +
    labs(x = "UMAP1", y = "UMAP2", color = "branch\nprobability")
  }), ncol = 2
) %>% ggsave(filename = "plots/SupFig1a.pdf", width = 16, height = 15)

plot_grid(
  plot.groups(
    vis = umap[cell.md$Run, ], add.labels = F,
    groups = cell.md$Donor, pal = pal$Donor
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Donor"),
  plot.groups(
    vis = umap[cell.md$Run, ], sy = 0.5,
    groups = cl[cell.md$Run, "Cluster"], pal = pal$Cluster[1:9]
  ) + labs(x = "UMAP1", y = "UMAP2", color = "CellSpace\nCluster"),
  ncol = 2
) %>% ggsave(filename = "plots/SupFig1b.pdf", width = 16, height = 5)


