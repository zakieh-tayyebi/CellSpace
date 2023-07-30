library(ArchR)
library(Matrix)
library(matrixStats)
library(irlba)
library(dplyr)

IterativeLSI <- function(
  counts, # a peak-by-cell count matrix of class dgCMatrix
  outDir = "IterativeLSI/", # where to save the LSI embedding and clustering results of every iteration
  iterations = 2, # number of LSI iterations (at least 2)
  varFeatures = 25000, # number of top variable peaks
  groupNames = NULL, # cell labels provided if you want the random sampling of cells to preserve group proportions
  sampleCellsPre = 10000, # An integer specifying the number of cells to sample in iterations prior to the last in order to perform a sub-sampled LSI and sub-sampled clustering
  sampleCellsFinal = NULL, # An integer specifying the number of cells to sample in order to perform a sub-sampled LSI in final iteration
  filterQuantile = 0.995,
  clusterParams = list(resolution = 2, maxClusters = 10),
  LSIMethod = 2,
  scaleTo = 10000,
  dimsToUse = 1:30,
  scaleDims = T,
  outlierQuantiles = c(0.02, 0.98),
  selectionMethod = "var",
  corCutOff = 0.75,
  filterBias = T,
  seed = 1
){
  if(class(counts) != "dgCMatrix")
    stop("'counts' must be a sparse matrix of class 'dgCMatrix'!")

  if(is.null(colnames(counts)) || any(duplicated(colnames(counts))) ||
     is.null(rownames(counts)) || any(duplicated(rownames(counts))))
    stop("'counts' must be an event-by-cell matrix with unique row and column names!")

  if(is.null(groupNames))
    groupNames <- rep("G1", ncol(counts))

  set.seed(seed)
  dir.create(outDir, showWarnings = F, recursive = T)

  # Handle number of cells to sample #####
  cellNames <- colnames(counts)
  if(!is.null(sampleCellsPre)){
    if(length(cellNames) < sampleCellsPre){
      sampleCellsPre <- NULL
    }
  }
  if(!is.null(sampleCellsFinal)){
    if(length(cellNames) < sampleCellsFinal){
      sampleCellsFinal <- NULL
    }
  }

  # LSI Iteration 1 #####
  nFeature <- varFeatures
  rmTop <- floor((1 - filterQuantile) * nrow(counts))
  message("\nSelecting top ", varFeatures, " most accessible features (excluding top ", rmTop, ") for LSI iteration 1.")

  totalAcc <- data.frame(
    rowSums = Matrix::rowSums(counts),
    idx = 1:nrow(counts),
    row.names = rownames(counts)
  )
  topIdx <- head(order(totalAcc$rowSums, decreasing = T), nFeature + rmTop)[-(1:rmTop)]
  topFeatures <- totalAcc[sort(topIdx), ]

  cellDF <- data.frame(
    cellDepth = log10(Matrix::colSums(counts) + 1),
    groupNames = groupNames,
    idx = 1:ncol(counts),
    row.names = colnames(counts)
  )

  j <- 1
  message("Running LSI iteration ", j, ".")

  outLSI <- .LSIPartialMatrix(
    counts = counts,
    featureDF = topFeatures,
    cellDF = cellDF,
    LSIMethod = LSIMethod,
    dimsToUse = dimsToUse,
    scaleTo = scaleTo,
    outlierQuantiles = outlierQuantiles,
    sampleCells = sampleCellsPre,
    projectAll = F
  )
  outLSI$scaleDims <- scaleDims
  saveRDS(outLSI, file = paste0(outDir, "/iteration", j, "-LSI.rds"))

  message("Identify clusters after LSI iteration ", j, ".")
  clusterDF <- .LSICluster(
    outLSI = outLSI,
    cellDF = cellDF,
    filterBias = filterBias,
    dimsToUse = dimsToUse,
    scaleDims = scaleDims,
    corCutOff = corCutOff,
    clusterParams = clusterParams,
    j = j
  )
  saveRDS(clusterDF, file = paste0(outDir, "/iteration", j, "-clusters.rds"))

  # LSI Iteration 2+ #####
  variableFeatures <- topFeatures

  while(j < iterations){
    j <- j + 1

    message("\nSelecting top ", varFeatures, " most variable features for LSI iteration ", j, ".")
    variableFeatures <- .identifyVarFeatures(
      counts = counts,
      outLSI = outLSI,
      clusterDF = clusterDF,
      prevFeatures = variableFeatures,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures
    )

    message("Running LSI iteration ", j, ".")
    outLSI <- .LSIPartialMatrix(
      counts = counts,
      featureDF = variableFeatures,
      cellDF = cellDF,
      LSIMethod = LSIMethod,
      scaleTo = scaleTo,
      dimsToUse = dimsToUse,
      outlierQuantiles = outlierQuantiles,
      sampleCells = if(j != iterations) sampleCellsPre else sampleCellsFinal,
      projectAll = j == iterations
    )
    outLSI$scaleDims <- scaleDims
    saveRDS(outLSI, file = paste0(outDir, "/iteration", j, "-LSI.rds"))

    clusterDF <- .LSICluster(
      outLSI = outLSI,
      cellDF = cellDF,
      dimsToUse = dimsToUse,
      scaleDims = scaleDims,
      corCutOff = corCutOff,
      filterBias = filterBias,
      j = j,
      clusterParams = clusterParams
    )
    saveRDS(clusterDF, file = paste0(outDir, "/iteration", j, "-clusters.rds"))
  }

  return(outLSI)
}

.identifyVarFeatures <- function(
  counts,
  outLSI,
  clusterDF,
  prevFeatures,
  totalAcc,
  scaleTo,
  firstSelection,
  selectionMethod,
  varFeatures
){

  groupMat <- sapply(SimpleList(split(clusterDF$cellNames, clusterDF$clusters)), function(group.idx){
    Matrix::rowSums(counts[totalAcc$idx, group.idx])
  }, simplify = T, USE.NAMES = T)

  nFeature <- varFeatures
  if(tolower(selectionMethod) == "var"){

    # Log-Normalize
    groupMat <- log2(t(t(groupMat) / colSums(groupMat)) * scaleTo + 1)
    feature.var <- matrixStats::rowVars(groupMat)
    idx <- sort(head(order(feature.var, decreasing = T), nFeature))
    variableFeatures <- totalAcc[idx, ]


  } else if(tolower(selectionMethod) == "vmr"){

    # Variance-to-Mean Ratio
    feature.vmr <- matrixStats::rowVars(groupMat) / rowMeans(groupMat)
    idx <- sort(head(order(feature.vmr, decreasing = T), nFeature))
    variableFeatures <- totalAcc[idx, ]


  } else stop("Feature selection method is not valid!")

  return(variableFeatures)
}

.LSIPartialMatrix <- function(
  counts,
  featureDF,
  cellDF,
  LSIMethod,
  dimsToUse,
  scaleTo,
  outlierQuantiles,
  sampleCells,
  projectAll,
  projection.steps = 10
){

  if(is.null(sampleCells)){

    message("- Computing LSI for all ", nrow(cellDF), " cells.")

    # Construct Matrix
    mat <- counts[featureDF$idx, cellDF$idx]
    if(!all(rownames(cellDF) == colnames(mat)))
      stop("Names of cells don't match!")

    # Perform LSI
    outLSI <- .computeLSI(
      mat = mat,
      LSIMethod = LSIMethod,
      scaleTo = scaleTo,
      nDimensions = max(dimsToUse),
      outlierQuantiles = outlierQuantiles
    )

  }
  else{

    message("- Computing partial LSI for ", sampleCells, " sampled cells.")

    sampledCellNames <- .sampleBySample(
      cellNames = rownames(cellDF),
      sampleNames = cellDF$groupNames,
      cellDepth = cellDF$cellDepth,
      sampleCells = sampleCells,
      outlierQuantiles = outlierQuantiles,
      factor = 2
    )
    sampledCells <- cellDF[sampledCellNames, ]
    mat <- counts[featureDF$idx, sampledCells$idx]
    if(!all(rownames(sampledCells) == colnames(mat)))
      stop("Names of sampled cells don't match!")

    print(table(sampledCells$groupNames))

    # Perform LSI on Partial Sampled Matrix
    outLSI <- .computeLSI(
      mat = mat,
      LSIMethod = LSIMethod,
      scaleTo = scaleTo,
      nDimensions = max(dimsToUse),
      outlierQuantiles = outlierQuantiles
    )

    if(projectAll){

      message("- Projecting LSI in ", projection.steps, " steps for the rest of the cells.")
      rest.idx <- setdiff(cellDF$idx, sampledCells$idx)
      steps <- c(floor(seq(1, length(rest.idx), length.out = projection.steps + 1))[-(projection.steps + 1)], length(rest.idx) + 1)
      pLSI <- lapply(1:projection.steps, function(i){
        message("-- Step ", i, " with ", steps[i + 1] - steps[i], " cells:")
        first <- steps[i]
        last <- steps[i + 1] - 1
        .projectLSI(mat = counts[featureDF$idx, rest.idx[first:last]], LSI = outLSI)
      }) %>% Reduce("rbind", .)
      matSVD <- rbind(outLSI$matSVD, pLSI)

      outLSI$exlcude.projection <- !(rownames(cellDF) %in% rownames(pLSI))
      outLSI$matSVD <- matSVD[rownames(cellDF), ]
    }
  }

  outLSI$LSIFeatures <- featureDF
  outLSI$corToDepth <- list(
    scaled = abs(cor(.scaleDims(outLSI[[1]]), cellDF[rownames(outLSI[[1]]), "cellDepth"]))[, 1],
    none = abs(cor(outLSI[[1]], cellDF[rownames(outLSI[[1]]), "cellDepth"]))[, 1]
  )

  return(outLSI)
}

.computeLSI <- function(
  mat,
  LSIMethod,
  scaleTo,
  nDimensions,
  outlierQuantiles
){ # TF IDF LSI adapted from flyATAC

    # Binarize Matrix
    mat@x[mat@x > 0] <- 1

    message("-- Cleaning up the matrix.")

    # Clean up zero columns
    colSm <- Matrix::colSums(mat)
    exclude <- colSm == 0
    if(any(exclude)){
      mat <- mat[, !exclude, drop = F]
      colSm <- colSm[!exclude]
    }

    # Remove outlying columns
    cn <- colnames(mat)
    filterOutliers <- 0
    if(!is.null(outlierQuantiles)){
      qCS <- quantile(colSm, probs = c(min(outlierQuantiles), max(outlierQuantiles)))
      idxOutlier <- which(colSm <= qCS[1] | colSm >= qCS[2])
      if(length(idxOutlier) > 0){
        matO <- mat[, idxOutlier, drop = F]
        mat <- mat[, -idxOutlier, drop = F]
        mat2 <- mat[, 1:10, drop = F] # A 2nd Matrix to Check Projection is Working
        colSm <- colSm[-idxOutlier]
        filterOutliers <- 1
      }
    }

    # Clean up zero rows
    rowSm <- Matrix::rowSums(mat)
    exclude2 <- rowSm == 0
    if(any(exclude2)){
      mat <- mat[!exclude2, ]
      rowSm <- rowSm[!exclude2]
    }

    message("-- Computing TF-IDF.")

    # TF - Normalize
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    # TF-IDF
    if(LSIMethod == 1 | tolower(LSIMethod) == "tf-log(idf)"){ #Adapted from Casanovich et al.

      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    } else if(LSIMethod == 2 | tolower(LSIMethod) == "log(tf-idf)"){ #Adapted from Stuart et al.

      #IDF
      idf   <- as(ncol(mat) / rowSm, "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * scaleTo + 1)


    } else if(LSIMethod == 3 | tolower(LSIMethod) == "log(tf-log(idf))"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    } else stop("Invalid LSI method!")

    message("-- Calculating SVD and LSI.")
    # Calculate SVD then LSI
    svd <- irlba::irlba(mat, nDimensions, nDimensions)
    svdDiag <- matrix(0, nrow = nDimensions, ncol = nDimensions)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI", 1:ncol(matSVD))

    # Return Object
    out <- SimpleList(
        matSVD = matSVD,
        rowSm = rowSm,
        nCol = length(colSm),
        exclude = exclude,
        exclude2 = exclude2,
        svd = svd,
        scaleTo = scaleTo,
        nDimensions = nDimensions,
        LSIMethod = LSIMethod,
        outliers = NA
    )

    if(filterOutliers == 1){
      message("-- Checking if LSI projection works.")
      # Quick Check LSI-Projection Works
      pCheck <- .projectLSI(mat = mat2, LSI = out)
      pCheck2 <- out[[1]][rownames(pCheck), ]
      pCheck3 <- lapply(1:ncol(pCheck), function(x){
        cor(pCheck[, x], pCheck2[, x])
      }) %>% unlist
      if(min(pCheck3) < 0.95) stop("cor<0.95 of re-projection!")

      message("-- Project LSI for ", ncol(matO), " outlying cells.")
      # Project LSI Outliers
      out$outliers <- colnames(matO)
      outlierLSI <- .projectLSI(mat = matO, LSI = out)
      allLSI <- rbind(out[[1]], outlierLSI)
      allLSI <- allLSI[cn, , drop = F] # Re-Order Correctly to original
      out[[1]] <- allLSI
    }

   return(out)
}

.projectLSI <- function(
  mat,
  LSI,
  returnModel = F
){

    # Get Same Features
    if(any(LSI$exclude2)) mat <- mat[!LSI$exclude2, ]

    # Binarize Matrix
    mat@x[mat@x > 0] <- 1

    message("--- Cleaning up the matrix.")

    # Clean up zero columns
    colSm <- Matrix::colSums(mat)
    exclude <- colSm == 0
    if(any(exclude)){
      mat <- mat[, !exclude]
      colSm <- colSm[!exclude]
    }

    message("--- Computing TF-IDF.")

    # TF - Normalize
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))

    # TF-IDF
    if(LSI$LSIMethod == 1 | tolower(LSI$LSIMethod) == "tf-log(idf)"){ #Adapted from Casanovich et al.

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-LogIDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    } else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){ #Adapted from Stuart et al.

      #IDF
      idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat

      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)


    } else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "log(tf-log(idf))"){

      #LogTF
      mat@x <- log(mat@x + 1)

      #LogIDF
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")

      #TF-IDF
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat


    } else stop("Invalid LSI method!")

    # Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat), arr.ind = T)
    if(length(idxNA) > 0) mat[idxNA] <- 0

    message("--- Calculating V and LSI.")

    # Calculate V
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1 / LSI$svd$d)

    # LSI Diagonal
    svdDiag <- matrix(0, nrow = LSI$nDimensions, ncol = LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI", 1:ncol(matSVD))

    if(returnModel){
        X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
        out <- list(matSVD = matSVD, V = V, X = X)
    } else out <- matSVD

  rm(mat, V, svdDiag, matSVD)

  return(out)
}

.LSICluster <- function(
  outLSI,
  cellDF,
  corCutOff,
  dimsToUse,
  scaleDims,
  clusterParams,
  j,
  filterBias
){

  if(scaleDims){
    dimsPF <- dimsToUse[which(outLSI$corToDepth$scaled[dimsToUse] <= corCutOff)]
  } else{
    dimsPF <- dimsToUse[which(outLSI$corToDepth$none[dimsToUse] <= corCutOff)]
  }

  if(length(dimsPF) != length(dimsToUse)){
    message("- Filtering ", length(dimsToUse) - length(dimsPF), " dims with cor>", corCutOff, " to log10(depth+1)")
  }
  if(length(dimsPF) < 2){
    stop("Dimensions to use (after filtering for correlation to depth) fewer than 2!")
  }

  # Time to compute clusters
  parClust <- clusterParams
  if(scaleDims){ parClust$input <- as.matrix(.scaleDims(outLSI$matSVD)[, dimsPF, drop = F])
  } else parClust$input <- as.matrix(outLSI$matSVD[, dimsPF, drop = F])

  if(filterBias){
    parClust$testBias <- T
    parClust$filterBias <- T
  }
  parClust$biasVals <- cellDF[rownames(outLSI$matSVD), "cellDepth"]

  clusters <- do.call(ArchR::addClusters, parClust) %>% suppressMessages()
  parClust$input <- NULL
  parClust$biasVals <- "cellDepth"
  nClust <- length(unique(clusters))
  message("- Identified ", nClust, " clusters:")
  message("-- method=\'Seurat\'")
  message("-- resolution=", parClust$resolution)
  message("-- maxClusters=", parClust$maxClusters)
  if(parClust$filterBias) message("-- biasVals=\'", parClust$biasVals, "\'")

  df <- DataFrame(
    cellNames = rownames(outLSI$matSVD), clusters = clusters,
    cellDF[rownames(outLSI$matSVD), c("idx", "groupNames", "cellDepth")]
  )
  metadata(df)$parClust <- parClust
  print(table(df[, c("groupNames", "clusters")]))

  return(df)
}

.scaleDims <- function(m, scaleMax = 2, min = -scaleMax, max = scaleMax, limit = F){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

.filterSample <- function(
    x = NULL,
    n = NULL,
    vals = x,
    outlierQuantiles = NULL,
    factor = 2
){
  if(!is.null(outlierQuantiles)){
    quant <- quantile(vals, probs = c(min(outlierQuantiles) / factor, 1 - ((1-max(outlierQuantiles)) / factor)))
    idx <- which(vals >= quant[1] & vals <= quant[2])
  } else idx <- seq_along(x)

  if(length(idx) >= n){ return(sample(x = x[idx], size = n))
  } else return(sample(x = x, size = n))
}

.sampleBySample <- function(
    cellNames = NULL,
    cellDepth = NULL,
    sampleNames = NULL,
    sampleCells = NULL,
    outlierQuantiles = NULL,
    factor = 2
){
  if(sampleCells < length(cellNames)){

    sampleN <- ceiling(sampleCells * table(sampleNames) / length(sampleNames))
    splitCells <- split(cellNames, sampleNames)
    splitDepth <- split(cellDepth, sampleNames)

    sampledCellNames <- lapply(seq_along(splitCells), function(x){
      .filterSample(
        x = splitCells[[x]],
        n = sampleN[names(splitCells)[x]],
        vals = splitDepth[[x]],
        outlierQuantiles = outlierQuantiles,
        factor = factor
      )
    }) %>% unlist %>% sort

    return(sampledCellNames)

  } else return(cellNames)
}

