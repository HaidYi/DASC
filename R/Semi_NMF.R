require(NMF)
require(foreach)
require(doParallel)

#'  Initialization of the semi-NMF
#'
#' @import NMF
#' @param model Object of class: NMFfit
#' @param target gene expression matrix
#' @return \code{model} The initial objective of class: NMFfit
#'
#' @author Haidong Yi, Ayush T. Raman

Ini_SemiNMF <- function(target, rank) {
  N <- ncol(target)
  km <- kmeans(t(target), rank)
  G0 <- matrix(0, nrow = N, ncol = rank)
  for (i in 1:N) {
    G0[i, km$cluster[i]] <- 1
  }
  G0 <- G0 + 0.2
  coef <- t(G0)
  basis <- target %*% G0 %*% solve(t(G0) %*% G0)
  return(list(coef=coef, basis=basis))
}


#' Update G in Semi-NMF
#'
#' @param X Data expression matrix need to be factorized
#' @param mf The basis matrix
#' @param mg The co-efficient matrix
#' @return \code{G} The basis matrix
#'
#' @details
#' By definition, G is a graph adjacency matrix. The \code{update_G} updates G 
#' after every iteration.
#'
#' @author Haidong Yi, Ayush T. Raman
#' @export
#' @examples
#' X <- matrix(1:12,nrow=4)
#' mf <- matrix(1:8,nrow=4)
#' mg <- matrix(1:6,ncol=2)
#' mg <- update_G(X,mf,mg)

update_G <- function(X, mf, mg) {
  n <- ncol(X)
  k <- ncol(mf)
  G <- matrix(nrow = n, ncol = k)
  GFFP <- t(mf) %*% mf
  GFFN <- t(mf) %*% mf
  XF <- t(X) %*% mf
  for (i in 1:k) {
    for (j in 1:k) {
      GFFP[i, j] <- 0.5 * (abs(GFFP[i, j]) + GFFP[i, j])
      GFFN[i, j] <- 0.5 * (abs(GFFN[i, j]) - GFFN[i, j])
    }
  }
  GFFP <- mg %*% GFFP
  GFFN <- mg %*% GFFN
  for (i in 1:n) {
    for (j in 1:k) {
      num <- 0.5 * (abs(XF[i, j]) + XF[i, j]) + GFFN[i, j]
      den <- 0.5 * (abs(XF[i, j]) - XF[i, j]) + GFFP[i, j]
      G[i, j] <- mg[i, j] * sqrt(num / den)
    }
  }
  G
}


#' Get the error of Semi-NMF using frobenius norm
#'
#' @details
#' This is a customerized function defined in terms of \code{\link[NMF]{nmf}}.
#' For more information, please go through the NMF vignette
#' \url{https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf}
#'
#' @import NMF
#' @param model Object of class: NMFfit
#' @param target gene expression matrix
#' @return The result of semi-NMF for the current iteration
#'
#' @author Haidong Yi, Ayush T. Raman
#'

Loss_Fro <- function(model, target) {
  norm( (target - model$basis %*% model$coef), type = "f")
}


#' Predict function of semi-NMF to determine which category
#' are the samples assigned.
#' 
#' @import NMF
#' @param x transpose of the coefficient matrix 
#' @return a factor vector representing the classification
#' @author Haidong Yi

pred.nmf <- function(x) {
    if( !is.matrix(x) ) stop("semi-NMF: only works on matrices")
    # for each column return the (row) index of the maxium
    return ( as.factor(apply(x, 2L, function(v) which.max(abs(v)))) )
}

#' Main function of semi-NMF
#'
#' @import NMF
#' @param target gene expression matrix
#' @param model Object of class: NMFfit
#' @param iternum Number of iterations
#' @return \code{model} the result from the semi-NMF algorithm
#' @author Haidong Yi, Ayush T. Raman

Semi_NMF <- function(target, rank, iternum = 100, tol = 0.01) {
  n <- 0
  model <- Ini_SemiNMF(target, rank)
  while (TRUE) {
    model$coef <- t(update_G(target, model$basis, t(model$coef)))
    if (length(which(model$coef < 1e-20) > 0)) break
    model$basis <- target %*% t(model$coef) %*%
      solve(model$coef %*% t(model$coef))
    if (n > iternum || Loss_Fro(model, target) < tol) break
    n <- n + 1
  }
  model
}


#'
#'
#' @param consen
#' @param rank
#'
#'

pred.consen <- function(consen, rank) {
  # build the tree from consensus matrix
  h <- hclust(as.dist(1-consen), method='average')
  # extract the membership from the tree
  cl <- stats::cutree(h, k = rank)
  as.factor(cl)
}


#'
#'
#' @param target
#' @param rank
#' @param nrun
#' @return
#' @export

sNMF <- function(target, rank, nrun=30) {
  n_cores <- detectCores(logical = FALSE) - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  n_sample <- ncol(target)
  
  snmf <- foreach(rank = rank, .inorder = TRUE, .packages = "NMF",
          .export = c("Semi_NMF", "Ini_SemiNMF", "update_G",
                      "pred.nmf", "Loss_Fro", "pred.consen")) %dopar% {
            consen <- matrix(0, nrow = n_sample, ncol = n_sample)
            for (i in 1:nrun) {
              model <- Semi_NMF(target, rank)
              sample_label <- pred.nmf(model$coef)
              consen <- consen + connectivity(sample_label)
            }
            consen <- consen / nrun
            list(consensus = consen, class = pred.consen(consen, rank), dispersion = dispersion(consen))
          }
  stopImplicitCluster()
  names(snmf) <- as.character(rank)
  snmf
}


