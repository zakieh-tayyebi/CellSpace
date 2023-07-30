library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

plot.groups <- function(
    object, vis, groups,
    cell.idx = rownames(object@cell.emb),
    mapping = aes(x = UMAP_1, y = UMAP_2),
    add.labels = T, pal = NULL,
    point.size = 1.5, label.size = 5,
    et = element_text(size = 15),
    sx = 0, sy = 0, nx = NULL, ny = NULL
){
  vis <- data.frame(cso@reductions[[vis]])[cell.idx, ]
  group.name <- groups
  groups <- cso@meta.data[cell.idx, groups]
  eb <- element_blank()
  gp <- ggplot(mapping = mapping) +
    geom_point(aes(color = groups), vis, size = point.size, alpha = 0.9) +
    labs(color = group.name) + theme_classic() +
    theme(
      axis.ticks = eb, axis.text = eb,
      axis.title = et, title = et,
      legend.title = et, legend.text = et
    )
  if(add.labels){
    G <- sort(unique(groups))
    NX <- rep(sx, length(G)); names(NX) <- G
    if(!is.null(nx)) NX[names(nx)] <- nx
    NY <- rep(sy, length(G)); names(NY) <- G
    if(!is.null(ny)) NY[names(ny)] <- ny
    vis.gp <- data.frame(matrix(
      NA, nrow = length(G), ncol = ncol(vis),
      dimnames = list(G, colnames(vis))
    ))
    for(g in G)
      for(col in colnames(vis))
        vis.gp[g, col] <- median(vis[groups == g, col])
    gp <- gp + geom_label(aes(color = G, label = G), vis.gp, show.legend = F,
                          size = label.size, nudge_x = NX, nudge_y = NY)
  }
  gp + scale_color_manual(values = pal)
}

plot.scores <- function(motif.score, cell.annot, pal, ...){
  Heatmap(
    motif.score, name = "Similarity\nz-score",
    col = colorRamp2(c(-1, 0, 1), c("#007e8c", "white", "#a3005c")),
    top_annotation = columnAnnotation(df = cell.annot, col = pal),
    column_title_rot = 90, show_column_names = F,
    row_title_rot = 0, cluster_rows = F,
    use_raster = T, ...
  )
}

