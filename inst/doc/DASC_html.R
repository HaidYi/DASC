## ----style, echo = FALSE, results = 'asis'---------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE, message=FALSE--------------------------------------
library(knitr)
opts_chunk$set(comment=NA, fig.align="center", warning=FALSE)

## Libraries
require(Biobase)
require(NMF)
require(cvxclustr)

## ---- message=FALSE, eval=TRUE---------------------------------------------
library(DASC)
data("esGolub")
samples <- c(20,21,28,30)
dat <- exprs(esGolub)[1:100,samples]
pdat <- pData(esGolub)[samples,]

## use nrun = 50 or more for better convergence of results
res <- DASC(edata = dat, pdata = pdat, factor = pdat$Cell, 
                        method = 'ama', type = 3, lambda = 1, 
                        rank = 2:3, nrun = 5, annotation='esGolub Dataset')
res

## ---- message=FALSE, eval=TRUE---------------------------------------------
## libraries
set.seed(99999)
library(DESeq2)
library(ggplot2)
library(pcaExplorer)

## dataset
rawCounts <- stanfordData$rawCounts
metadata <- stanfordData$metadata

## ---- message=FALSE, eval=TRUE---------------------------------------------
## Using a smaller dataset
idx <- which(metadata$tissue %in% c("adipose", "adrenal", "sigmoid"))
rawCounts <- rawCounts[,idx]
metadata <- metadata[idx,]

## ---- message=FALSE, eval=TRUE---------------------------------------------
head(rawCounts)
head(metadata)

## ---- message=FALSE, eval=TRUE---------------------------------------------
## Normalizing the dataset using DESeq2
dds <- DESeqDataSetFromMatrix(rawCounts, metadata, design = ~ species+tissue)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)
lognormalizedCounts <- log2(dat + 1)

## ---- message=FALSE, eval=TRUE---------------------------------------------
## PCA plot using 
rld.dds <- rlog(dds)
pcaplot(rld.dds, intgroup=c("tissue","species"), ntop=1000, pcX=1, pcY=2)

## ---- message=FALSE, eval=TRUE---------------------------------------------
res <- DASC(edata = dat, pdata = metadata, factor = metadata$tissue,
                method = 'ama', type = 3, lambda = 1, rank = 2:3, nrun = 10,
                annotation = 'Stanford Dataset')

## ---- message=FALSE, eval=TRUE---------------------------------------------
## Consensus plot
annotation <- data.frame(Batch = res$`2`$class, Tissue = as.character(metadata$tissue)) 
consensusmap(res$`2`$consensus, annCol = annotation, main = "rank = 2")

## ---- message=FALSE, eval=TRUE---------------------------------------------
## Batches -- dataset has 6 batches
sample.clust <- data.frame(sample.name = colnames(lognormalizedCounts), 
                            clust = as.vector(res$`2`$class), 
                            batch = metadata$seqBatch)
ggplot(data = sample.clust, aes(x=c(1:6), y=clust, color=factor(clust))) + 
    geom_point(size = 4) + xlab("Sample Number") + ylab("Cluster Number")

## ---- message=FALSE, eval=FALSE--------------------------------------------
#  
#  ## not running this part of the code for building package
#  ## Using entire dataset
#  rawCounts <- stanfordData$rawCounts
#  metadata <- stanfordData$metadata
#  dds <- DESeqDataSetFromMatrix(rawCounts, metadata, design = ~ species+tissue)
#  dds <- estimateSizeFactors(dds)
#  dat <- counts(dds, normalized = TRUE)
#  lognormalizedCounts <- log2(dat + 1)
#  
#  ## PCA Plot
#  rld.dds <- rlog(dds)
#  pcaplot(rld.dds, intgroup=c("tissue","species"), ntop = nrow(rld.dds),
#          pcX = 1, pcY = 2)
#  
#  ## Running DASC
#  res <- DASC(edata = dat, pdata = metadata, factor = metadata$tissue,
#                  method = 'ama', type = 3, lambda = 1, rank = 2:10,
#                  nrun = 100, annotation = 'Stanford Dataset')
#  
#  ## Consensus plot
#  annotation <- data.frame(Batch = res$`10`$class, Tissue = metadata$tissue)
#  consensusmap(res$`10`$consensus, annCol = annotation, main = "rank = 10")
#  
#  ## Clustering samples based on batches
#  sample.clust <- data.frame(sample.name = colnames(lognormalizedCounts),
#                              clust = as.vector(res$`10`$class),
#                              batch = metadata$seqBatch)
#  ggplot(data = sample.clust, aes(x=c(1:26), y=clust, color=factor(clust))) +
#      geom_point(size = 4) + xlab("Sample Number") + ylab("Cluster Number")

## ---- out.width = "800px", echo=FALSE--------------------------------------
knitr::include_graphics("PCA_Plot_Fig1.png")
knitr::include_graphics("Consensus_Plot_Fig2.png")
knitr::include_graphics("ClusterNumber_Fig4.png")

## ---- message=FALSE--------------------------------------------------------
sessionInfo()

