require(Biobase)
require(NMF)
require(cvxclustr)

#' Batch factor detection via DASC (Data-adaptive Shrinkage and Clustering-DASC)
#'
#' @param edata the normalized target matrix, a data.frame The row is gene,
#' the column is sample
#' @param pdata Phenotypic data summarizes information about the samples
#' @param factor A factor vector which controls the convex clustering
#' @param method Algorithm to use: 'admm' or 'ama'
#' @param type An integer indicating the norm used: 1 = 1-norm 2 = 2-norm
#' 3 = 2-norm^2
#' @param lambda A double number A regularization parameter in the convex
#' optimization
#' @param rank integer sequence
#' @param nrun the iteration numbers of Semi-NMF
#' @param spanning parameter is assigned as false
#' @param annotation An annotation of the dataset
#' @return outputs the result of semi-NMF. It classifies each sample to its
#' batch factor.
#' @details
#' The \code{DASC} function is the main function of our algorithm DASC
#' (Data-adaptive Shrinkage and Clustering-DASC) package. The DASC includes 
#' two main steps
#' \itemize{ \item Data-adaptive shrinkage using convex clustering shrinkage
#' (Implemented by convex optimization.);
#' \item Extract batch factors using matrix factorization.
#' }
#'
#' @import Biobase
#' @import cvxclustr
#' @import NMF
#'
#' @export
#' @examples
#' library(NMF)
#' library(cvxclustr)
#' library(Biobase)
#' dat <- data.frame(matrix(rnbinom(n=200, mu=100, size=1/0.5), ncol=4))
#' pdat <- data.frame(sample = colnames(dat), type = c(rep('A',2), rep('B',2)))
#' rownames(pdat) <- colnames(dat)
#' res <- DASC(edata=dat, pdata=pdat, factor=pdat$type, method='ama', type=3,
#' lambda = 1, rank = 2, nrun = 50, spanning = FALSE,
#' annotation='simulated dataset')
#'
#' @author Haidong Yi, Ayush T. Raman
#'
#' @seealso \code{\link[cvxclustr]{cvxclust_path_ama}} and
#' \code{\link[cvxclustr]{cvxclust_path_admm}} for the detailed algorithm
#'

DASC <- function(edata, pdata, factor, method="ama", type=3, lambda, rank,
                    nrun, spanning=FALSE, annotation){
    if (!is.null(type) && !(type %in% c(1, 2, 3))) {
        stop("type must be '1', '2', '3', or NULL")
    }
    if (!is.null(method) && !(method %in% c("ama", "admm"))) {
        stop("method must be 'ama', 'admm', or NULL")
    }
    edata <- as.matrix(edata)
    Zero <- apply(edata, 1, sd)
    Zero.num <- which(Zero < 0.001)
    if (length(Zero.num) > 0) {
        names(Zero.num) <- NULL
        edata <- edata[-Zero.num, ]
    }
    edata <- log(1 + edata)
    if (type == 3) {
        Laplace <- trans_Laplace(as.factor(factor))
        Udata <- edata %*% solve(diag(ncol(edata)) + lambda * Laplace)
        Udata <- as.matrix(Udata)
        Bdata <- edata - Udata
    } else {
        ADJ <- trans_ADJ(as.factor(factor))
        if (spanning) {
            ADJ <- Sptree(ADJ)
        }
        w <- adj2vector(ADJ, nrow(ADJ))
        sol <- cvxclust(edata, w, lambda, method = method, type = type)
        Bdata <- edata - sol$U[[1]]
    }
    Zero <- apply(Bdata, 1, sd)
    Zero.num <- which(Zero < 0.001)
    if (length(Zero.num) > 0) {
        names(Zero.num) <- NULL
        Bdata <- Bdata[-Zero.num, ]
    }
    metadata <- data.frame(labelDescription = names(pdata),
                            row.names = names(pdata))
    pdata <- new("AnnotatedDataFrame", data = pdata, varMetadata = metadata)
    data_set <- ExpressionSet(assayData = edata, phenoData = pdata,
                                annotation = annotation)
    data.nmf <- nmf(data_set, rank, Semi_NMF, nrun = nrun, .opt = "v",
                        objective = Loss_Fro, seed = Ini_SemiNMF, mixed = TRUE)
    return(data.nmf)
}

#' Transform the adjacency matrix to a vector
#'
#' @param Adjacency the adjacency matrix of factor
#' @param n number of samples
#' @return \code{w} the vector of the adjacency matrix
#' @author Haidong Yi, Ayush T. Raman
#'
#' @export
#' @examples
#' W <- matrix(c(0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0),nrow=4)
#' w <- adj2vector(W,4)
#'

adj2vector <- function(Adjacency, n) {
    w <- double(n * (n - 1) / 2)
    iter <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            w[iter] <- Adjacency[i, j]
            iter <- iter + 1
        }
    }
    return(w)
}

#' Representing node in this subtype
#'
#' @param v the index of the node
#' @param X the saved vector with the information of the parent of every node
#' @return \code{r} the parent index of the node
#' @author Haidong Yi, Ayush T. Raman
#'
#' @export
#' @examples
#' nodes <- c(2,3,4,4)
#' get_father(2, nodes)
#'

get_father <- function(v, X) {
    r <- v
    while (X[r] != r) {
        r <- X[r]
    }
    i <- v
    j <- numeric()
    while (i != r) {
        j <- X[i]
        X[i] <- r
        i <- j
    }
    return(r)
}


#' Get Spanning tree from adjacency matrix
#'
#' @param ADJ the adjacency matrix of the factor
#' @return \code{ADJ} the spaning tree of the adjacency matrix
#' @author Haidong Yi, Ayush T. Raman
#'
#' @export
#' @examples
#' W <- matrix(c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0), nrow=4)
#' Sptree(W)
#'

Sptree <- function(ADJ) {
    Rownum <- nrow(ADJ)
    Colnum <- ncol(ADJ)
    father <- numeric()
    for (i in 1:Rownum) {
        father[i] <- i
    }
    for (i in 2:Rownum) {
        for (j in 1:(i - 1)) {
            if (ADJ[i, j] > 0) {
                father <- merge(i, j, father)
            }
        }
    }
    sum <- numeric()
    pre <- numeric()
    j <- 1
    for (i in 1:Rownum) {
        pre[i] <- get_father(i, father)
        if (pre[i] == i) {
            sum[j] <- i
            j <- j + 1
        }
    }
    ADJ <- matrix(0, Rownum, Colnum)
    for (i in 1:Colnum) {
        if (i != pre[i]) {
            ADJ[i, pre[i]] <- 1
            ADJ[pre[i], i] <- 1
        }
    }
    ADJ
}

#' Combine two trees into one
#'
#' @param x the index of the node
#' @param y the index of the node
#' @param X the saved vector with the information of the parent of every node
#' @return \code{X} A updated X vector with updates on father of every node
#' @author Haidong Yi, Ayush T. Raman
#' @details
#'
#' During the traversal of the graph matrix, merge function joins two 
#' disjoint sets into a single subset. It is a union step of Disjoint-set 
#' algorithm by Bernard A. Galler and Michael J. Fischer. For further details, 
#' please refer to: 
#' \url{https://en.wikipedia.org/wiki/Disjoint-set_data_structure}
#'

merge <- function(x, y, X) {
    fx <- get_father(x, X)
    fy <- get_father(y, X)
    if (fx < fy) {
        X[fx] <- fy
    } else {
        X[fy] <- fx
    }
    return(X)
}
