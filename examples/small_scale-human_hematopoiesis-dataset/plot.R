# #####
setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/")

library(ggplot2)
library(grid)
library(ggrepel)
library(ggdendro)
library(ggseqlogo)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

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
          legend.title = et, legend.text = et)
}

cell.md <- read.table(
  file = "sample-info.tsv",
  sep = "\t", header = T, check.names = F
) %>% subset(!discard)
rownames(cell.md) <- cell.md$Run

pal <- readRDS("plots/palette.rds")
cell.md$Cell_type <- factor(cell.md$Cell_type, levels = names(pal$Cell_type))


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
cs.scores <- read.csv(
  file = "CellSpace-results/var_tiles/TF_motif-embedding/selected-motif-scores.csv",
  row.names = 1, check.names = F
)
TFs <- rownames(cs.scores)

nx <- rep(0, length(TFs)); names(nx) <- TFs; ny <- nx
nx['EBF1'] <- -0.5; nx['PAX5'] <- 0.2; ny['EBF1'] <- 0.2;  ny['PAX5'] <- -0.5
ggsave(
  filename = "plots/Fig2a.pdf", width = 8, height = 5,
  plot.groups(
    vis = umap[cell.md$Run, ], groups = cell.md$Cell_type,
    pal = pal$Cell_type, add.labels = F
  ) +
  geom_point(data = umap[TFs, ], color = "black", shape = 17, size = 1.5) +
  geom_label_repel(
    mapping = aes(label = TFs), data = umap[TFs, ],
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

pdf("plots/Fig3c.pdf", width = 10, height = 7)
gp <- gpar(fontsize = 10)
ta <- columnAnnotation(
  df = cbind(
    cell.md[, c("Cell_type", "Donor")],
    Cluster = cl[cell.md$Run, "Cluster"]
  ), col = pal
)
Heatmap(
  as.matrix(cs.scores[TFs, cell.md$Run]), name = "Similarity\nz-score",
  col = colorRamp2(c(-1, 0, 1), c("#007e8c", "white", "#a3005c")),
  top_annotation = ta, column_split = cell.md$Cell_type,
  column_title_rot = 90, show_column_names = F,
  row_title_rot = 0, cluster_rows = F,
  row_names_gp = gp, row_names_max_width = max_text_width(TFs, gp = gp),
  use_raster = T
) %>% draw(
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  padding = unit(c(2, 2, 10, 2), "mm")
)
dev.off()

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


# #####
archr.res <- read.csv(
  file = "ArchR-results/UMAP_and_Clusters.csv",
  row.names = 1, header = T
)
archr.res$Cluster <- factor(archr.res$Cluster)

ggsave(
  filename = "plots/Fig2d.pdf", width = 8, height = 5,
  plot.groups(
    vis = archr.res[cell.md$Run, 1:2], groups = cell.md$Cell_type,
    sy = 1, nx = c(MPP = -2),
    pal = pal$Cell_type, add.labels = T
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Cell type")
)

plot_grid(
  plot.groups(
    vis = archr.res[cell.md$Run, 1:2], add.labels = F,
    groups = cell.md$Donor, pal = pal$Donor
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Donor"),
  plot.groups(
    vis = archr.res[cell.md$Run, 1:2], sy = 0.5,
    groups = archr.res[cell.md$Run, "Cluster"], pal = pal$Cluster[1:14]
  ) + labs(x = "UMAP1", y = "UMAP2", color = "ArchR\nCluster"),
  ncol = 2
) %>% ggsave(filename = "plots/SupFig1c.pdf", width = 16, height = 5)


# #####
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

mt.pal <- c(
  HG = "#F9E600", ASW = "#FF7617", NMI = "#FF29D8", ARI = "#C50000",
  bASW = "#00ED32", bNMI = "#0002A1", GC = "#00DEEC", kBET = "#8818FF"
)
sc <- scale_color_manual(values = mt.pal)
sp <- rep(c(1, 2), c(4, 4)); names(sp) <- names(mt.pal)
ss <- scale_shape_manual(values = sp)
column <- rep(c("Biological Conservation", "Batch Correction"), each = 4)
names(column) <- names(mt.pal)

eval_metrics <- read.csv(
  file = "benchmarking/metrics-summary.csv",
  row.names = 1, check.names = F
)
eval_metrics <- lapply(colnames(eval_metrics), function(metric){
  data.frame(
    embedding = emb.names[rownames(eval_metrics)],
    metric = metric,
    value = eval_metrics[, metric]
  )
}) %>% do.call(what = rbind)
eval_metrics$metric <- factor(eval_metrics$metric, levels = names(mt.pal))

mapping <- aes(
  x = factor(embedding, levels = as.character(emb.names)),
  y = value, group = metric, color = metric, shape = metric
)
th <- theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
            axis.title = element_blank())
ggsave(
  filename = "plots/SupFig2d.pdf", width = 4, height = 6,
  plot =
    ggplot(eval_metrics, mapping = mapping) +
    geom_point() + geom_line() +
    sc + ss + ylim(c(0, 1)) +
    facet_wrap(~ column[metric], ncol = 1) +
    theme_bw() + th
)

cell.md.exUnk <- subset(cell.md, Cell_type != "unknown")
clusters <- read.csv(
  file = "benchmarking/Clusters.csv",
  row.names = 1, check.names = F, colClasses = "factor"
)[cell.md.exUnk$Run, ]
umap <- sapply(names(emb.names), function(emb){
  read.csv(paste0("benchmarking/UMAP/", emb, ".csv"), row.names = 1)[cell.md.exUnk$Run, ]
}, simplify = F)

n <- length(emb.names) + 1
lg.df <- cbind(umap[[1]], cell.md.exUnk)
ps <- 1
et <- element_text(size = 15)
thm <- theme(title = et, axis.title = et, legend.title = et, legend.text = et)
cg <- guides(colour = guide_legend(override.aes = list(size=10), title.hjust = 0.5, ncol=2))

ct.pl <- lapply(names(emb.names), function(emb){
  plot.groups(
    vis = umap[[emb]], groups = cell.md.exUnk$Cell_type, pal = pal$Cell_type[-10],
    add.labels = F, point.size = ps
  ) + labs(x = "UMAP1", y = "UMAP2", title = emb.names[emb]) +
    thm + theme(legend.position = "none")
})
ct.pl[[n]] <- get_legend(
  plot.groups(
    vis = lg.df[, 1:2], groups = lg.df$Cell_type, pal = pal$Cell_type[-10]
  ) + labs(color = "Cell type") + thm + cg
)
ggsave(
  plot = plot_grid(plotlist = ct.pl, ncol = 2),
  width = 15, height = 35, filename = "plots/SupFig2a.pdf"
)

dn.pl <- lapply(names(emb.names), function(emb){
  plot.groups(
    vis = umap[[emb]], groups = cell.md.exUnk$Donor, pal = pal$Donor,
    add.labels = F, point.size = ps
  ) + labs(x = "UMAP1", y = "UMAP2", title = emb.names[emb]) +
    thm + theme(legend.position = "none")
})
dn.pl[[n]] <- get_legend(
  plot.groups(
    vis = lg.df[, 1:2], groups = lg.df$Donor, pal = pal$Donor
  ) + labs(color = "Donor") + thm + cg
)
ggsave(
  plot = plot_grid(plotlist = dn.pl, ncol = 2),
  width = 15, height = 35, filename = "plots/SupFig2b.pdf"
)

cl.pl <- lapply(names(emb.names), function(emb){
  plot.groups(
    vis = umap[[emb]], groups = clusters[, emb], pal = pal$Cluster[1:9],
    add.labels = F, point.size = ps
  ) + labs(x = "UMAP1", y = "UMAP2", title = emb.names[emb]) +
    thm + theme(legend.position = "none")
})
cl.pl[[n]] <- get_legend(
  plot.groups(
    vis = lg.df[, 1:2], groups = clusters[, 1], pal = pal$Cluster[1:9]
  ) + labs(color = "Cluster") + thm + cg
)
ggsave(
  plot = plot_grid(plotlist = cl.pl, ncol = 2),
  width = 15, height = 35, filename = "plots/SupFig2c.pdf"
)


# #####
load("CellSpace-results/var_tiles/denovo-motifs/top-10mers-clustering.rda")
load("CellSpace-results/var_tiles/denovo-motifs/denovo-PWMs.rda")
umap.denovo <- read.csv(
  file = "CellSpace-results/var_tiles/denovo-motifs/UMAP-cells_and_motifs.csv",
  row.names = 1
)
C <- table(cl.10mers)
C <- as.integer(names(C[C >= 4]))

dd <- dendro_data(dend.10mers)
dd$labels$cluster <- factor(cl.10mers[dd$labels$label], levels = C)
align <- lapply(C, function(cl){
  X <- subset(dd$labels, cluster == cl)$x
  data.frame(
    x = c(mean(X), sort(X)),
    y = c(-3.5, rep(-1.5, length(X))),
    label = c(
      paste("motif", match(cl, C)),
      as.character(align.cl.10mers[[cl]])
    ),
    cluster = rep(factor(cl, levels = C), length(X) + 1)
  )
}) %>% do.call(what = rbind)
ggsave(
  plot = ggplot(mapping = aes(x = x, y = y)) +
    geom_segment(aes(xend = xend, yend = yend), dd$segments) +
    geom_text(
      aes(label = label, color = cluster), rbind(dd$labels, align),
      size = 2, show.legend = F, hjust = -0.1, family = "mono"
    ) +
    scale_y_reverse(expand = c(0.5, 0.2)) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    coord_flip() + theme_void(),
  filename = "plots/SupFig1d.pdf", width = 5, height = 30
)

dummy <- lapply(C, function(cl){
  m <- paste0("motif", match(cl, C))
  pwm <- denovo.PWMs[[cl]]
  ggsave(
    plot = ggseqlogo(pwm, method = 'prob') + ggtitle(m),
    width = ncol(pwm) / 2 + 0.5, height = 2,
    filename =  paste0("plots/Fig2e-logos/denovo-", m, ".pdf")
  )
}) %>% do.call(what = c)

cisbp.pwms <- readRDS("CisBP-Homo_sapiens_2022_09_04_8-24_pm/CisBP-PWMs.rds")
known.motif <- c(
  M08841_2.00 = "M08841_2.00 (CEBPA, CEBPB, CEBPE)",
  M07997_2.00 = "M07997_2.00 (PAX2, PAX5, ...)",
  M09795_2.00 = "M09795_2.00 (TCF3, TCF4, TCF12)",
  M02210_2.00 = "M02210_2.00 (HOXA4, HOXB4, TLX1, ...)",
  M03336_2.00 = "M03336_2.00 (IRF4, IRF8)",
  M08126_2.00 = "M08126_2.00 (GATA1, GATA2, GATA3, ...)",
  M09617_2.00 = "M09617_2.00 (RARA, RORA, RORC, ...)",
  M03018_2.00 = "M03018_2.00 (FOXP3, FOXI1, FOXJ1, ...)",
  M10596_2.00 = "M10596_2.00 (GATA1, GATA2, GATA3, ...)",
  M08207_2.00 = "M08207_2.00 (GATA1, GATA2, GATA3, ...)",
  M09430_2.00 = "M09430_2.00 (TBX4, TBX5, TBX15, ...)",
  M08192_2.00 = "M08192_2.00 (PRDM1)",
  M05572_2.00 = "M05572_2.00 (ESRRA, RARA, RORA, RORC, ...)"
)
rc <- rep(F, length(known.motif)); names(rc) <- names(known.motif)
rc[c("M07997_2.00", "M02210_2.00", "M10596_2.00", "M09795_2.00")] <- T
dummy <- lapply(names(known.motif), function(id){
  pwm <- (exp(as.matrix(cisbp.pwms[[id]])) / 4) %>%
    apply(2, function(x){ x / sum(x) })
  if(rc[id]){ # reverse complement the known motif
    pwm.rc <- pwm[, ncol(pwm):1]
    rownames(pwm.rc) <- c(A = "T", C = "G", G = "C", `T` = "A")[rownames(pwm.rc)]
    pwm <- pwm.rc
  }
  ggsave(
    plot = ggseqlogo(pwm, method = 'prob') + ggtitle(known.motif[id]),
    width = ncol(pwm) / 2 + 0.5, height = 2,
    filename =  paste0("plots/Fig2e-logos/known-", id, ".pdf")
  )
}) %>% do.call(what = c)

mn <- paste("motif", 1:29)
ggsave(
  filename = "plots/Fig2e-UMAP.pdf", width = 8, height = 5,
  plot.groups(
    vis = umap.denovo[cell.md$Run, ], groups = cell.md$Cell_type,
    pal = pal$Cell_type, add.labels = F
  ) + geom_point(data = umap.denovo[mn, ], color = "black", shape = 17, size = 1) +
  geom_label_repel(
    mapping = aes(label = mn), data = umap.denovo[mn, ],
    size = 3, fontface = "bold", min.segment.length = 0, label.padding = unit(0.4, "mm"),
    max.overlaps = Inf, segment.color = "black", segment.size = unit(0.2, "mm")
  ) + labs(x = "UMAP1", y = "UMAP2", color = "Cell type")
)

