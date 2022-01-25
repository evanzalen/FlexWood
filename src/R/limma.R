#' ---
#' title: "Batch correction ZE and SE"
#' author: "Elena van Zalen"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(gplots)
  library(ggplot2)
  library(here)
  library(hyperSpec)
  library(limma)
  library(matrixStats)
  library(parallel)
  library(plotly)
  library(tibble)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Data
load(here("data/analysis/salmon/dds.rda"))
levels(dds$BioID) <- c("T89_STA","T89_ROT")

#' # Normalise
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

#' # Batch effect
#' ## Estimation
mat <- assay(vsd)
mat <- mat-min(mat)
mat <- mat[rowSums(mat) > 0, ]

test <- removeBatchEffect(mat, dds$TW, design = model.matrix(~dds$BioID))
#pca
#heatmap

#' We select FMG and ZE from mature seeds (29Seed), that correspond to
#' Stage B8 and B9 of the ZE
sel <- dds$BioID == "T89_ROT" | dds$TW %in% c("yes","no")

mat <- mat[,sel]

BioID = dds$BioID[sel]

contrasts(BioID) <- contr.sum(levels(BioID))

design <- model.matrix(~BioID)

fit <- lmFit(mat, design)

#' beta is the actual batch correction for each gene
beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0

#' cross-validation with the removeBatchEffect function
stopifnot(all((mat - beta %*% t(design[,-1,drop=FALSE]))[1:6,1:6] == 
                removeBatchEffect(mat,dds$BioID[sel])[1:6,1:6]))

#' ## Correction 
load(here("analysis/salmon/ZE-SE-dds.rda"))
levels(dds.sz$Experiment) <- c("ZE","SE")
vsd <-varianceStabilizingTransformation(dds.sz,blind=FALSE)

fullbatch <- vsd$Experiment
contrasts(fullbatch) <- contr.sum(levels(fullbatch))
fullbatch <- model.matrix(~fullbatch)[, -1, drop = FALSE]

assay(vsd) <- assay(vsd) - beta
