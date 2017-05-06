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
    H <- basis(model)
    G <- t(coef(model))
    norm( (target - H %*% t(G)), type = "f")
}
