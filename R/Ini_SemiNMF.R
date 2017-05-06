#'  Initialization of the semi-NMF
#'
#' @import NMF
#' @param model Object of class: NMFfit
#' @param target gene expression matrix
#' @return \code{model} The initial objective of class: NMFfit
#'
#' @author Haidong Yi, Ayush T. Raman

Ini_SemiNMF <- function(model, target) {
    N <- ncol(target)
    ans <- kmeans(t(target), nbasis(model))
    G0 <- matrix(0, nrow = N, ncol = nbasis(model))
    for (i in 1:N) {
        G0[i, ans$cluster[i]] <- 1
    }
    G0 <- G0 + 0.2
    coef(model) <- t(G0)
    basis(model) <- target %*% G0 %*% solve(t(G0) %*% G0)
    return(model)
}
