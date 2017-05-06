#' Outputs Adjacency matrix from the factor vector
#'
#' @param col_data A factor vector
#' @return \code{adjacency} Adjacency matrix of \code{col_data}
#' @export
#' @author Haidong Yi, Ayush T. Raman
#'
#' @examples
#' batch.factor <- c(rep('human',13),rep('mouse',13))
#' batch.factor <- as.factor(batch.factor)
#' adj <- trans_Laplace(batch.factor)
#'

trans_ADJ <- function(col_data) {
    adjacency <- matrix(0, nrow = length(col_data),
                            ncol = length(col_data))
    iter <- length(levels(col_data))
    for (i in 1:iter) {
        FACTOR <- which(col_data == levels(col_data)[i])
        N <- length(FACTOR)
        for (j in 1:(N - 1)) {
            for (k in (j + 1):N) {
                adjacency[FACTOR[j], FACTOR[k]] <- 1
            }
        }
    }
    adjacency <- adjacency + t(adjacency)
    return(adjacency)
}
