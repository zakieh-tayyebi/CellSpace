library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)
library(Matrix)
library(dplyr)

get.tile.matrix <- function(archr.obj, genome){
  tile.mtx <- getMatrixFromProject(archr.obj, useMatrix = "TileMatrix", binarize = T)
  var.tiles <- archr.obj@reducedDims$IterativeLSI$LSIFeatures[, -3]

  tile.width <- archr.obj@reducedDims$IterativeLSI$tileSize
  var.tiles.gr <- GRanges(
    seqinfo = seqinfo(genome), seqnames = var.tiles$seqnames, strand = "+",
    ranges = IRanges(start = var.tiles$start, end = var.tiles$start + tile.width - 1,
                     names = paste0("tile", var.tiles$idx))
  )

  return(list(
    var.tiles = var.tiles, var.tiles.gr = var.tiles.gr, genome = genome,
    var.tile.mtx = assays(tile.mtx)$TileMatrix[match(var.tiles, tile.mtx@elementMetadata), ]
  ))
}

archr.obj <- loadArchRProject("ArchROutput/")
mtx <- get.tile.matrix(archr.obj, genome = BSgenome.Hsapiens.UCSC.hg19)

dir.create("CellSpace-inputs/")
mtx$var.tiles.gr %>% getSeq(x = mtx$genome) %>%
  writeXStringSet(filepath = "CellSpace-inputs/top50K-variable-tiles.fa")
writeMM(t(mtx$var.tile.mtx), file = "CellSpace-inputs/cell-by-tile-bin.mtx")
write(archr.obj$cellNames, file = "CellSpace-inputs/cell-names.txt", ncolumns = 1)
