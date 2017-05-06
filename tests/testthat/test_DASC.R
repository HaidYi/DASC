library(DASC)
context("DASC functionality")

test_that("adj2vector", {
    W <- matrix(c(0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 4)
    w <- adj2vector(W, 4)
})

test_that("DASC", {
    dat <- data.frame(matrix(rnbinom(n = 200, mu = 100, size = 1 / 0.5),
                                ncol = 4))
    pdat <- data.frame(sample = colnames(dat), type = c(rep("A", 2),
                                                            rep("B", 2)))
    rownames(pdat) <- colnames(dat)
    if (.Platform$OS.type == "windows") {
        res <- DASC(edata = dat, pdata = pdat, factor = pdat$type,
                     method = "ama", type = 3, lambda = 1, rank = 2,
                     nrun = 1, annotation = "simulated dataset")
    }else{
        res <- DASC(edata = dat, pdata = pdat, factor = pdat$type,
                     method = "ama", type = 3, lambda = 1, rank = 2,
                     nrun = 50, annotation = "simulated dataset")
    }
})

test_that("trans_Laplace", {
    factors <- data.frame(Sample = c("Sample1", "Sample2", "Sample3",
                                      "Sample4"), type = c(rep("A", 2),
                                                            rep("B", 2)))
    trans_Laplace(as.factor(factors$type))
})

test_that("Sptree", {
    W <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0), nrow = 4)
    Sptree(W)
})
