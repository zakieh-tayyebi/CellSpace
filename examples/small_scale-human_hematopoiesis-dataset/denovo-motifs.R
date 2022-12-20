setwd("~/CellSpace/examples/small_scale-human_hematopoiesis-dataset/CellSpace-results/var_tiles/")

library(CellSpace)
library(Seurat)
library(Biostrings)
library(PWMEnrich)
library(msa)
library(kmer)
library(dplyr)
library(parallel)

# CellSpace embedding for cells and DNA 8-mers:
cl <- read.csv("CellSpace-Clusters.csv", row.names = 1, colClasses = "factor")
cso <- CellSpace(
  emb.file = "CellSpace-embedding.tsv", # from CellSpace-train-all_peaks.sh
  cell.names = rownames(cl), meta.data = cl
)

# map all possible DNA 10-mers to the CellSpace embedding:
#  out of each k-mer and its reverse complement, only one is included
dna.10mers <- readLines("denovo-motifs/all-DNA-10mers.txt")
emb.10mers <- mclapply(dna.10mers, function(seq){
  CellSpace::DNA_sequence_embedding(object = cso, seq = seq)
}, mc.cores = 8) %>%
  do.call(what = rbind) %>% data.frame(row.names = dna.10mers)

# identify K=50 nearest 10-mers to each cell:
source("https://raw.githubusercontent.com/satijalab/seurat/master/R/clustering.R")
nn.10mers <- NNHelper(
  data = emb.10mers, query = cso@cell.emb,
  method = "annoy", metric = "cosine", k = 50
)

top.10mers <- lapply(levels(cso$Cluster), function(c){
  cells <- cso$Cluster == c
  min.cn <- ceiling(sum(cells) * 0.2)

  # find 10-mers that were in the top nearest neighbors to
  #   at least 20% of cells in the cluster:
  nn.10mers <- dna.10mers[c(nn.10mers@nn.idx[cells, ])] %>% table()
  top.kmers <- names(nn.10mers)[nn.10mers >= min.cn]

  # add reverse complements:
  top.kmers <- Biostrings::DNAStringSet(top.kmers) %>% Biostrings::reverseComplement() %>%
    as.character() %>% toupper() %>%
    c(top.kmers) %>% unique() %>% sort() %>% Biostrings::DNAStringSet()
  names(top.kmers) <- as.character(top.kmers)

    return(top.kmers)
}) %>% do.call(what = c) %>% unique()

# cluster top 10-mers:
dend.10mers <- kmer::cluster(
  x = strsplit(as.character(top.10mers), split = ""),
  k = 5, residues = "DNA"
) %>% as.hclust()
cl.10mers <- cutree(tree = dend.10mers, h = 0.5)

# discard clusters with <4 k-mers:
C <- table(cl.10mers)
C <- as.integer(names(C[C >= 4]))

# align k-mers of each cluster:
align.cl.10mers <- lapply(sort(unique(cl.10mers)), function(cl){
  if(sum(cl.10mers == cl) < 4) return(NULL)
  nucleotides <- c("A", "C", "G", "T")
  cl.kmers <- names(cl.10mers)[cl.10mers == cl]
  msa::msa(Biostrings::DNAStringSet(cl.kmers))
})

# compute the position weight matrix for each cluster (a de novo motif):
denovo.PWMs <- lapply(align.cl.10mers, function(align){
  if(is.null(align)) return(NULL)
  nucleotides <- c("A", "C", "G", "T")
  cs.pwm <- apply(as.matrix(align), 2, function(pos){
    tb <- table(pos)[nucleotides]
    if(sum(tb, na.rm = T) < 2) return(rep(NA, 4))
    return(tb / sum(tb, na.rm = T))
  })
  na.pos <- which(apply(cs.pwm, 2, function(x){ all(is.na(x)) }))
  keep.pos <- setdiff(1:ncol(cs.pwm), na.pos)
  cs.pwm <- cs.pwm[, min(keep.pos):max(keep.pos)]
  cs.pwm[is.na(cs.pwm)] <- 0
  rownames(cs.pwm) <- nucleotides
  return(cs.pwm)
})

save(dend.10mers, cl.10mers, file = "denovo-motifs/top-10mers-clustering.rda")
save(align.cl.10mers, denovo.PWMs, file = "denovo-motifs/denovo-PWMs.rda")

# map de novo motifs to the CellSpace embedding:
motif.emb <- lapply(C, function(cl){
  # 10-mers in the corresponding cluster:
  kmers <- gsub("-", "", as.character(align.cl.10mers[[cl]]))
  # find each k-mer or its reverse complement in the 10-mers embedded earlier:
  kmers <- ifelse(
    test = kmers %in% rownames(emb.10mers),
    yes = kmers,
    no = Biostrings::DNAStringSet(kmers) %>% Biostrings::reverseComplement() %>% as.character()
  )
  # average of embedding vectors:
  colMeans(emb.10mers[kmers, ])
}) %>% do.call(what = rbind) %>%
  data.frame(row.names = paste("motif", 1:length(C)))

# UMAP of cells and de novo motifs:
umap.denovo <- Seurat::RunUMAP(
  object = rbind(cso@cell.emb, as.matrix(motif.emb)),
  metric = "cosine", n.neighbors = 20,
  min.dist = 0.2, spread = 1
)@cell.embeddings %>% write.csv(file = "denovo-motifs/UMAP-cells_and_motifs.csv")

# similarity of known CisBP motifs to CellSpace de novo motifs:
cisbp.pwms <- readRDS("../../CisBP-Homo_sapiens_2022_09_04_8-24_pm/CisBP-PWMs.rds")
pairs <- c(M08126_2.00 = 67, M09617_2.00 = 94, M08207_2.00 = 136)
PFM <- function(pwm){
  pfm <- apply(pwm, 2, function(x){ as.integer(x * 100) })
  class(pfm) <- c(class(pfm), "integer")
  rownames(pfm) <- rownames(pwm)
  return(pfm)
}
dummy <- lapply(names(pairs), function(mid){
  cl <- pairs[mid]
  pwm1 <- denovo.PWMs[[cl]]
  pwm2 <- (exp(as.matrix(cisbp.pwms[[mid]])) / 4) %>%
    apply(2, function(x){ x / sum(x) })
  cat(
    "similarity of ",
    paste("de novo motif", match(cl, C)), " and CisBP ", mid, ": ",
    PWMEnrich::motifSimilarity(m1 = PFM(pwm1), m2 = PFM(pwm2)),
    "\n", sep = ""
  )
})

