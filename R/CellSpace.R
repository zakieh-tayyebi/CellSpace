
#' The CellSpace Class
#'
#' The CellSpace class stores CellSpace embedding and
#' related information needed for performing downstream analyses.
#'
#' @slot project title of the project
#' @slot emb.file the .tsv output of CellSpace containing the embedding matrix for cells and k-mers
#' @slot cell.emb the embedding matrix for cells
#' @slot kmer.emb the embedding matrix for k-mers
#' @slot motif.emb the embedding matrix for TF motifs
#' @slot meta.data data frame containing meta-information about each cell
#' @slot dim the dimensions of the CellSpace embeddings
#' @slot k the length of DNA k-mers
#' @slot similarity the similarity function in hinge loss
#' @slot p the embedding of an entity equals the sum of its M feature embedding vectors divided by M^p
#' @slot label cell label prefix
#' @slot neighbors list containing nearest neighbor graphs
#' @slot reductions list containing dimensional reductions
#' @slot misc list containing miscellaneous objects
#'
#' @name CellSpace-class
#' @rdname CellSpace-class
#' @exportClass CellSpace
#'
setClass("CellSpace", slots = list(
  project = "character",
  emb.file = "character",
  cell.emb = "matrix",
  kmer.emb = "matrix",
  motif.emb = "list",
  meta.data = "data.frame",
  dim = "integer",
  k = "integer",
  similarity = "character",
  p = "numeric",
  label = "character",
  neighbors = "list",
  reductions = "list",
  misc = "list"
))

setMethod(
  "show", "CellSpace",
  function(object){
    cat("An object of class \"CellSpace\"\n")
    if(object@project != "") cat("Project:", object@project, "\n")
    cat("CellSpace model: ", object@dim, "-dimensional embedding for ",
        nrow(object@cell.emb), " cells and all DNA k-mers (k=", object@k, ")\n", sep = ""
    )
    if(length(object@motif.emb) > 0){
      cat("TF motif embedding matrices: ")
      cat(names(object@motif.emb), sep = ", ")
      cat("\n")
    }
    if(ncol(object@meta.data) > 0){
      cat("Cell meta-data: ")
      if(length(colnames(object@meta.data)) > 10){
        cat(head(colnames(object@meta.data), 3), sep = ", ")
        cat(", ..., ")
        cat(tail(colnames(object@meta.data), 2), sep = ", ")
      } else cat(colnames(object@meta.data), sep = ", ")
      cat("\n")
    }
    if(length(object@neighbors) > 0){
      cat("Nearest neighbor graphs created: ")
      cat(names(object@neighbors), sep = ", ")
      cat("\n")
    }
    if(length(object@reductions) > 0){
      cat("Dimensional reductions calculated: ")
      cat(names(object@reductions), sep = ", ")
      cat("\n")
    }
    if(length(object@misc) > 0){
      cat("miscellaneous: ")
      cat(names(object@misc), sep = ", ")
      cat("\n")
    }
  }
)

setMethod(
  "$", "CellSpace",
  function(x, name){
    if(name %in% colnames(x@meta.data)){
      return(x@meta.data[, name])
    } else return(NULL)
  }
)

setMethod(
  "$<-", "CellSpace",
  function(x, name, value){
    x@meta.data[, name] <- value
    return(x)
  }
)

#' CellSpace
#'
#' Generates an object from the CellSpace class.
#'
#' @param emb.file the .tsv output of CellSpace containing the embedding matrix for cells and k-mers
#' @param cell.names vector of unique cell names
#' @param meta.data a \code{data.frame} containing meta-information about each cell
#' @param project title of the project
#' @param similarity the similarity function in hinge loss
#' @param p the embedding of an entity equals the sum of its M feature embedding vectors divided by M^p
#' @param label cell label prefix
#'
#' @return a new \code{CellSpace} object
#'
#' @name CellSpace
#' @rdname CellSpace
#' @export
#'
CellSpace <- function(
  emb.file,
  cell.names = NULL,
  meta.data = NULL,
  project = NULL,
  similarity = "cosine",
  p = 0.5,
  label = "__label__"
){
  project <- ifelse(is.null(project), "", project)

  emb.file <- normalizePath(emb.file)
  emb <- read.table(emb.file, header = F, row.names = 1, sep = "\t")
  dim <- ncol(emb)
  colnames(emb) <- paste0("CS", 1:dim)

  cell.label <- paste0(label, "C")
  cell.labels <- rownames(emb)[grep(cell.label, rownames(emb))]
  cell.idx <- sort(as.integer(gsub(cell.label, "", cell.labels)))
  cell.emb <- as.matrix(emb[paste0(cell.label, cell.idx), ])
  if(!is.null(cell.names)){
    if(length(cell.names) != nrow(cell.emb) || any(duplicated(cell.names))){
      warning("Unable to change cell names: \'cell.names\' must be a character vector with ", nrow(cell.emb), " unique values!")
    } else rownames(cell.emb) <- cell.names
  }

  kmer.emb <- as.matrix(emb[!grepl(label, rownames(emb)), ])
  k <- unique(nchar(rownames(kmer.emb)))
  if(length(k) > 1) stop("Varying k-mer lengths!")

  if(is.null(meta.data)){
    meta.data <- data.frame(row.names = rownames(cell.emb), check.rows = F, check.names = F)
  } else {
    if(!is.data.frame(meta.data)) meta.data <- as.data.frame(meta.data, check.names = F)
    rownames(meta.data) <- rownames(cell.emb)
    if(nrow(meta.data) != nrow(cell.emb))
      stop("The number of rows in \'meta.data\' does not match the cell embedding matrix!")
  }

  new("CellSpace",
      project = project,
      emb.file = emb.file,
      cell.emb = cell.emb,
      kmer.emb = kmer.emb,
      motif.emb = list(),
      meta.data = meta.data,
      dim = dim,
      k = k,
      similarity = similarity,
      p = p,
      label = label,
      neighbors = list(),
      reductions = list(),
      misc = list()
  )
}

#' find_neighbors
#'
#' Built a nearest neighbor graph and shared nearest neighbor graph from the CellSpace embedding.
#'
#' @importFrom Seurat FindNeighbors
#'
#' @param object a \code{CellSpace} object
#' @param n.neighbors the number of nearest neighbors for the KNN algorithm
#' @param emb the embedding matrix used to create the nearest neighbor graphs
#' @param emb.name prefix for the graph names that will be added to the \code{neighbors} slot
#' @param ... arguments passed to \code{Seurat::FindNeighbors}
#'
#' @return a \code{CellSpace} object containing nearest neighbor and shared nearest neighbor graphs in the \code{neighbors} slot
#'
#' @name find_neighbors
#' @rdname find_neighbors
#' @export
#'
find_neighbors <- function(
  object,
  n.neighbors = 30,
  emb = object@cell.emb,
  emb.name = "cells",
  ...
){
  graphs <- FindNeighbors(
    emb, distance.matrix = F, k.param = n.neighbors, l2.norm = F,
    nn.method = "annoy", annoy.metric = object@similarity, ...
  )
  names(graphs) <- paste(emb.name, names(graphs), sep = "_")
  object@neighbors[names(graphs)] <- graphs
  return(object)
}

#' find_clusters
#'
#' Find clusters in a nearest neighbor graph built from the CellSpace embedding.
#'
#' @importFrom Seurat FindClusters
#'
#' @param object a \code{CellSpace} object
#' @param graph name of the nearest neighbor graph in the \code{neighbors} slot used to find clusters
#' @param ... arguments passed to \code{Seurat::FindClusters}
#'
#' @return a \code{CellSpace} object with the cell clusters added to the \code{meta.data} slot
#'
#' @name find_clusters
#' @rdname find_clusters
#' @export
#'
find_clusters <- function(object, graph = "cells_snn", ...){
  if(!(graph %in% names(object@neighbors)))
    stop("\'", graph, "\' not available! Run \'find_neighbors\' to create nearest neighbor graphs.")
  cl <- FindClusters(object@neighbors[[graph]], ...)
  colnames(cl) <- gsub("res\\.", "Clusters.res_", colnames(cl))
  object@meta.data[, colnames(cl)] <- cl
  return(object)
}

#' merge_small_clusters
#'
#' Merge cells from small clusters with the nearest clusters.
#'
#' @param object a \code{CellSpace} object
#' @param clusters a vector of cluster labels or the name of a column, in the \code{meta.data} slot, containing cluster labels
#' @param min.cells any cluster with fewer cells than \code{min.cells} will be merged with the nearest cluster
#' @param graph a nearest neighbor graph, or the name of a nearest neighbor graph in the \code{neighbors} slot, used to find clusters
#'
#' @return new cluster labels
#'
#' @name merge_small_clusters
#' @rdname merge_small_clusters
#' @export
#'
merge_small_clusters <- function(
    object,
    clusters,
    min.cells = 10,
    graph = "cells_snn",
    seed = 1
){
  if(length(graph) == 1 && graph %in% names(object@neighbors)){
    graph <- object@neighbors[[graph]]
  } else if(class(graph) != "Graph")
    stop("'graph' must be a nearest neighbor graph of class 'Graph'")

  if(length(clusters) == 1 && clusters %in% colnames(object@meta.data)){
    clusters <- object@meta.data[, clusters]
  } else if(length(clusters) == nrow(graph)){
    if(class(clusters) != "factor")
      clusters <- factor(clusters)
  } else stop("'clusters' must be the name of a column in 'object@meta.data', or a vector of the same length as the number of cells.")

  small.clusters <- names(which(table(clusters) < min.cells)) %>% intersect(levels(clusters))
  cluster_names <- setdiff(levels(clusters), small.clusters)
  connectivity <- vector(mode = "numeric", length = length(cluster_names))
  names(connectivity) <- cluster_names

  new.ids <- clusters
  for(i in small.clusters){
    i.cells <- which(clusters == i)
    for(j in cluster_names){
      j.cells <- which(clusters == j)
      subSNN <- graph[i.cells, j.cells]
      if(is.object(subSNN)){
        connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
      } else connectivity[j] <- mean(subSNN)
    }
    set.seed(seed = seed)
    m <- max(connectivity, na.rm = T)
    mi <- which(connectivity == m, arr.ind = T)
    closest_cluster <- sample(names(connectivity[mi]), 1)
    new.ids[i.cells] <- closest_cluster
  }

  return(factor(new.ids, levels = cluster_names))
}

#' run_UMAP
#'
#' Compute a UMAP embedding from the CellSpace embedding.
#'
#' @importFrom Seurat RunUMAP
#'
#' @param object a \code{CellSpace} object
#' @param emb the embedding matrix used to compute the UMAP embedding
#' @param graph name of the nearest neighbor graph in the \code{neighbors} slot used to compute the UMAP embedding, used only if \code{emb} is NULL
#' @param name name of the lower-dimensional embedding that will be added to the \code{reductions} slot
#' @param ... arguments passed to \code{Seurat::RunUMAP}
#'
#' @return a \code{CellSpace} object containing a UMAP embedding in the \code{reductions} slot
#'
#' @name run_UMAP
#' @rdname run_UMAP
#' @export
#'
run_UMAP <- function(
  object,
  emb = object@cell.emb,
  graph = NULL,
  name = "cells_UMAP",
  ...
){
  if(!is.null(graph)){
    if(graph %in% names(object@neighbors)){
      umap <- RunUMAP(object = object@neighbors[[graph]], assay = "CellSpace", ...)
    } else stop("\'", graph, "\' not available! Run \'find_neighbors\' to create nearest neighbor graphs.")
  } else if(!is.null(emb)){
    umap <- RunUMAP(
      object = emb,
      metric = object@similarity,
      assay = "CellSpace",
      ...
    )
  } else stop("'emb' or 'graph' must be provided!")
  object@reductions[[name]] <- umap@cell.embeddings
  return(object)
}

#' cosine_similarity
#'
#' Cosine similarity in the embedding space.
#'
#' @param x an embedding matrix
#' @param y NULL, in which case \code{y=x}, or an embedding matrix with compatible dimensions to \code{x}
#'
#' @name cosine_similarity
#' @rdname cosine_similarity
#' @export
#'
cosine_similarity <- function(x, y = NULL){
  if(!is.matrix(x)) x <- matrix(x, nrow = 1)
  if(!is.null(y) && !is.matrix(y)) y <- matrix(y, nrow = 1)

  normx <- sqrt(rowSums(x ^ 2))
  if(is.null(y)){
    y <- x
    normy <- normx
  } else normy <- sqrt(rowSums(y ^ 2))

  s <- tcrossprod(x, y) / (normx %o% normy)
  return(s)
}

#' embedding_distance
#'
#' Distance in the embedding space based on cosine similarity.
#'
#' @param x an embedding matrix
#' @param y NULL, in which case \code{y=x}, or an embedding matrix with compatible dimensions to \code{x}
#' @param distance the distance metric, either 'cosine' or 'angular', to compute from the cosine similarity
#'
#' @name embedding_distance
#' @rdname embedding_distance
#' @export
#'
embedding_distance <- function(x, y = NULL, distance = "cosine"){
  s <- cosine_similarity(x = x, y = y)
  idx1 <- which(s >  1); s[idx1] <-  1
  idx2 <- which(s < -1); s[idx2] <- -1

  if(distance == "cosine"){ ds <- 1 - s
  } else if(distance == "angular"){ ds <- acos(s) / pi
  } else stop("The distance metric must be \'cosine\' or \'angular\'!\n")
  return(ds)
}

#' DNA_sequence_embedding
#'
#' Map a DNA sequence to the embedding space.
#'
#' @importFrom Biostrings reverseComplement DNAStringSet
#'
#' @param object a \code{CellSpace} object
#' @param seq a DNA sequence
#'
#' @name DNA_sequence_embedding
#' @rdname DNA_sequence_embedding
#' @export
#'
DNA_sequence_embedding <- function(object, seq){
  if(!is.character(seq)) seq <- as.character(seq)
  sl <- nchar(seq)
  if(sl < object@k){
    warning("The sequence \'", seq, "\' is shorter than CellSpace k-mers (k=", object@k, ")!")
    return(rep(NA, object@dim))
  }

  b <- 1:(sl - object@k + 1)
  kmers <- substring(text = toupper(seq), first = b, last = b + object@k - 1)
  kmers.rc <- toupper(as.character(reverseComplement(DNAStringSet(kmers))))
  kmer.rownames <- c(kmers[kmers %in% rownames(object@kmer.emb)],
                     kmers.rc[kmers.rc %in% rownames(object@kmer.emb) & kmers.rc != kmers])
  kn <- length(kmer.rownames)

  if(kn != sl - object@k + 1) stop("Invalid DNA sequence \'", seq, "\'!")
  if(kn == 1) return(object@kmer.emb[kmer.rownames, ])
  else return(colSums(object@kmer.emb[kmer.rownames, ]) / (kn ^ object@p))
}

#' motif_embedding
#'
#' Map a motif to the embedding space.
#'
#' @importFrom Biostrings as.matrix
#'
#' @param object a \code{CellSpace} object
#' @param PWM \code{PFMatrix} or \code{PWMatrix}
#'
#' @name motif_embedding
#' @rdname motif_embedding
#' @export
#'
motif_embedding <- function(object, PWM){
  freq.mtx <- as.matrix(PWM)
  consensus <- paste(rownames(freq.mtx)[apply(freq.mtx, 2, which.max)], collapse = "")
  DNA_sequence_embedding(object = object, seq = consensus)
}

#' add_motif_db
#'
#' Compute the CellSpace embedding and activity scores of transcription factor motifs.
#'
#' @param object a \code{CellSpace} object
#' @param PWMs \code{PFMatrixList} or \code{PWMatrixList}
#' @param db.name the name of the transcription factor motif database
#'
#' @return a \code{CellSpace} object containing the motif embedding matrix, in the \code{motif.emb} slot, and the corresponding similarity Z-scores, in the \code{misc} slot
#'
#' @name add_motif_db
#' @rdname add_motif_db
#' @export
#'
add_motif_db <- function(object, PWMs, db.name){
  object@motif.emb[[db.name]] <- lapply(PWMs, function(motif.pwm){
    motif_embedding(object, PWM = motif.pwm)
  }) %>% do.call(what = rbind) %>% na.omit()

  object@misc[[paste(db.name, "scores", sep = "_")]] <- cosine_similarity(
    x = object@cell.emb,
    y = object@motif.emb[[db.name]]
  ) %>% scale()

  return(object)
}

