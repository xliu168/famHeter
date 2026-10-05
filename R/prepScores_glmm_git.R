prepScores <- function (Z, formula, family = gaussian(), SNPInfo = NULL, snpNames = "Name", 
                        aggregateBy = "gene", kins = NULL, mito.kins = NULL, sparse = TRUE, 
                        data = parent.frame(), verbose = FALSE) 
{
  if (is.null(SNPInfo)) {
    warning("No SNP Info file provided")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda", sep = "/"))
    aggregateBy = "SKATgene"
  }
  
  if (!is.null(kins) || !is.null(mito.kins)) {
    n <- if (!is.null(kins)) nrow(kins) else nrow(mito.kins)
    if (!is.null(kins)) stopifnot(nrow(kins) == nrow(Z))
    if (!is.null(mito.kins)) stopifnot(nrow(mito.kins) == nrow(Z))
    
    if (!is.null(kins)) {
      data$id <- if (is.null(colnames(kins))) 1:ncol(kins) else colnames(kins)
    } else {
      data$id <- if (is.null(colnames(mito.kins))) 1:ncol(mito.kins) else colnames(mito.kins)
    }
    kinList <- list()
    if (!is.null(kins))      kinList <- c(kinList, list(kins))
    if (!is.null(mito.kins)) kinList <- c(kinList, list(mito.kins))
    
    kinList_conv <- lapply(kinList, function(K) as(as.matrix(K), "dgCMatrix"))
    
    if (sparse && !is.null(kins)) {
      kins_thresh <- kins
      kins_thresh[kins_thresh < 2 * 2^(-6)] <- 0
      kins_thresh <- Matrix::forceSymmetric(kins_thresh)
      kinList_conv[[1]] <- as(as.matrix(kins_thresh), "dgCMatrix")
    }
    
       nullmodel <- GMMAT::glmmkin(fixed = formula, data = data, 
                                kins = kinList_conv, id = "id", family = family, method = "REML")
       res <- as.vector(nullmodel$scaled.residuals)
    X1 <- nullmodel$X

    P <- nullmodel$P
    Sigma_i <- nullmodel$Sigma_i
    Sigma_iX <- nullmodel$Sigma_iX
    cov_beta <- nullmodel$cov
    
    n <- nrow(X1)
  } 
  else {
    nullmodel <- stats::glm(formula = formula, family = family, data = data)
    res <- stats::residuals(nullmodel, type = "response")
    X1 <- stats::model.matrix(nullmodel)
    n <- nrow(X1)
    P <- NULL
    Sigma_i <- NULL
    Sigma_iX <- NULL
    cov_beta <- NULL
  }
  
  invisible(seqMeta:::check_format_skat(Z, SNPInfo, nullmodel, aggregateBy, snpNames))
  mysnps <- colnames(Z)
  SNPInfo[, aggregateBy] <- as.character(SNPInfo[, aggregateBy])
  SItoZ <- which(colnames(Z) %in% SNPInfo[, snpNames])
  which.snps.Z <- colnames(Z) %in% SNPInfo[, snpNames]
  ZtoSI <- match(SNPInfo[, snpNames], mysnps[which.snps.Z])
  nsnps <- sum(!is.na(ZtoSI))
  if (nsnps == 0) stop("No column names in Z match SNP names in the SNP Info file!")
  
  if (verbose) {
    cat("\n Scoring... Progress:\n")
    pb <- txtProgressBar(min = 0, max = nsnps, style = 3)
    pb.i <- 0
  }
  
  maf0 <- colMeans(Z, na.rm = TRUE)[which.snps.Z]
  maf0[is.nan(maf0)] <- -1
  maf <- maf0[ZtoSI]
  names(maf) <- SNPInfo[, snpNames]
  
  scores <- apply(Z[, which.snps.Z, drop = FALSE], 2, function(z) {
    if (any(is.na(z))) {
      if (all(is.na(z))) z <- rep(0, length(z))
      z[is.na(z)] <- mean(z, na.rm = TRUE)
    }
    if (verbose) {
      assign("pb.i", get("pb.i", env) + 1, env)
      if (get("pb.i", env) %% ceiling(nsnps / 100) == 0)
        setTxtProgressBar(get("pb", env), get("pb.i", env))
    }
    sum(res * z)
  })[ZtoSI]
  scores[is.na(scores)] <- 0
  names(scores) <- SNPInfo[, snpNames]
  if (verbose) close(pb)
  scores[maf == 0] <- 0
  maf[!(SNPInfo[, snpNames] %in% colnames(Z))] <- -1
  scores <- split(scores, SNPInfo[, aggregateBy])
  maf <- split(maf, SNPInfo[, aggregateBy])
  
  if (is.null(kins) && is.null(mito.kins)) {
    X1 <- sqrt(nullmodel$family$var(nullmodel$fitted)) * X1
  }
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())()
  }
    if (is.null(kins) && is.null(mito.kins)) {
    AX1 <- with(
      svd(X1),
      v[, d > 0, drop = FALSE] %*%
        ((1 / d[d > 0]) * t(u[, d > 0, drop = FALSE]))
    )
  }
  
  ngenes <- length(unique(SNPInfo[, aggregateBy]))
  if (verbose) {
    cat("\n Calculating covariance... Progress:\n")
    pb <- txtProgressBar(min = 0, max = ngenes, style = 3)
    pb.i <- 0
  }
  
  re <- tapply(SNPInfo[, snpNames], SNPInfo[, aggregateBy], function(snp.names) {
    inds <- match(snp.names, colnames(Z))
    mcov <- matrix(0, length(snp.names), length(snp.names))
    if (length(na.omit(inds)) > 0) {
      Z0 <- Z[, na.omit(inds), drop = FALSE]
      storage.mode(Z0) <- "numeric"
      if (!is.null(kins) || !is.null(mito.kins)) {
        Z0m <- as.matrix(Z0)

          if (!is.null(P)) {
          PZ <- P %*% Z0m
          mcov[!is.na(inds), !is.na(inds)] <- as.matrix(
            crossprod(Z0m, as.matrix(PZ))
          )
        } else {
         
          if (is.null(Sigma_i) || is.null(Sigma_iX) || is.null(cov_beta)) {
            stop(
              "GMMAT null model did not provide P or the ",
              "Sigma_i/Sigma_iX/cov components required for ",
              "score covariance calculation."
            )
          }

          Sigma_iZ <- Sigma_i %*% Z0m
          first_term <- crossprod(Z0m, as.matrix(Sigma_iZ))

          GSigma_iX <- crossprod(Z0m, as.matrix(Sigma_iX))
          second_term <- GSigma_iX %*% as.matrix(cov_beta) %*% t(GSigma_iX)

          mcov[!is.na(inds), !is.na(inds)] <- as.matrix(
            first_term - second_term
          )
        }
      } else {
        Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted)) * Z0
        mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0) %*% X1) %*% (AX1 %*% Z0)
      }
    }
    rownames(mcov) <- colnames(mcov) <- snp.names
    if (verbose) {
      assign("pb.i", get("pb.i", env) + 1, env)
      if (get("pb.i", env) %% ceiling(ngenes / 100) == 0)
        setTxtProgressBar(get("pb", env), get("pb.i", env))
    }
    Matrix::forceSymmetric(Matrix::Matrix(mcov, sparse = TRUE))
  }, simplify = FALSE)
  
  if (!is.null(kins) || !is.null(mito.kins)) {
  	sey <- 1
	} else {
  	sey <- sqrt(var(res) * (nrow(X1) - 1) / (nrow(X1) - ncol(X1)))
  	if (family$family == "binomial") sey <- 1
	}  
  for (k in 1:length(re)) {
    re[[k]] <- list(scores = scores[[k]], cov = re[[k]], n = n, maf = maf[[k]], sey = sey)
  }
  attr(re, "family") <- family$family
  class(re) <- "seqMeta"
  return(re)
}
