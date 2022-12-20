# ArchR functions

library(ArchR)
library(Matrix)
library(matrixStats)
library(irlba)
library(dplyr)


computeLSI <- function(
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

projectLSI <- function(
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

scaleDims <- function(
    m,
    scaleMax = 2,
    min = -scaleMax,
    max = scaleMax,
    limit = F
){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

