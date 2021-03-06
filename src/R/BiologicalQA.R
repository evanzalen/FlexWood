#' ---
#' title: "FlexWood Biological QA"
#' author: "Elena van Zalen"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(pander)
  library(plotly)
  library(RColorBrewer)
  library(scatterplot3d)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information
#' ```{r CHANGEME1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#' ```
samples <- readr::read_csv(here("doc/Samples.csv"),
                      col_types = cols(
                        col_character(),
                        col_factor(),
                        col_factor(),
                        col_factor(),
                        col_character()
                      )) 

#' tx2gene translation table
#' ```{r CHANGEME2,eval=FALSE,echo=FALSE}
#' # This file is necessary if your species has more than one transcript per gene.
#' #
#' # It should then contain two columns, tab delimited, the first one with the transcript
#' # IDs and the second one the corresponding
#' ```
tx2gene <- suppressMessages(read_delim(here("reference/annotation/tx2gene.tsv.gz"), delim = "\t",
                                 col_names = c("TXID", "GENE")))

#' # Raw data
filelist <- list.files(here("data/salmon"), 
                          recursive = TRUE, 
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info
#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are in the same order as 
#' # described in the samples object. If not, then they need to be sorted by making sure 
#' they match to the order of the filelist. This was the case for this filelist so the first
#' line orders the rows from the meta data accordingly.
#' ```
samples <- samples[match(sub("_sortmerna_trimmomatic", "", basename(dirname(filelist))),
                          samples$sampleID),]
stopifnot(all(match(sub("_sortmerna_trimmomatic", "", basename(dirname(filelist))),
                    samples$sampleID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$Name

#' Read the expression at the gene level
#' ```{r CHANGEME4, eval = FALSE, echo = FALSE}
#' If the species has only one transcript per gene, replace with the following
#' counts <- suppressMessages(round(tximport(files = filelist, type = "salmon", txOut=TRUE)$counts))
#' ```
counts <- suppressMessages(round(tximport(files = filelist, 
                                  type = "salmon",
                                  tx2gene = tx2gene)$counts))

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts), digits = 1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by CHANGEME
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x = colnames(counts), y = colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat, aes(x, y, fill = samples$Flexure)) + geom_col() + 
  scale_y_continuous(name = "reads") +
  theme(axis.text.x = element_text(angle = 90, size = 4), axis.title.x = element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value = log10(rowMeans(counts))), aes(x = value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name = "mean raw counts (log10)")

#' The same is done for the individual samples. 
#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Flexure = samples$Flexure[match(ind, samples$Name)]) %>% 
  mutate(Tissue = samples$Tissue[match(ind, samples$Name)])

ggplot(dat, aes(x = values, group = ind, col = Tissue)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name ="per gene raw counts (log10)")

#' ## Export
#dir.create(here("data/analysis/salmon"), showWarnings = FALSE, recursive = TRUE)
#write.csv(counts, file = here("data/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Flexure)
colnames(dds) <- as.character(colData(dds)$Name)
saveRDS(dds, file = here("data/analysis/salmon/dds.rds"))
save(dds, file = here("data/analysis/salmon/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main = "Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
saveRDS(vst, file = here("data/analysis/salmon/vst.rds"))
write.csv(vst, file = here("data/analysis/salmon/normalised-gene-expression_data_with S339.csv"))

#' * Validation
#' 
#' The variance stabilisation worked adequately
meanSdPlot(vst[rowSums(vst) > 0, ])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2, ]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar = 2

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables as needed.
#' ```
nlevel = nlevels(dds$Tissue) * nlevels(dds$Flexure) # 4 levels

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x = 1:length(percent), y = cumsum(percent)), aes(x = x, y = y)) +
  geom_line() + scale_y_continuous("variance explained (%)", limits = c(0, 100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept = nvar, colour = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = cumsum(percent)[nvar], colour = "red", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = nlevel, colour = "orange", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = cumsum(percent)[nlevel], colour = "orange", linetype = "dashed", size = 0.5)
  
#' ### 3 first dimensions
mar = c(5.1, 4.1, 4.1, 2.1)

#' The PCA shows that a large fraction of the variance is 
#' explained by both variables.
scatterplot3d(pc$x[, 1],
              pc$x[, 2],
              pc$x[, 3],
              xlab = paste("Comp. 1 (", percent[1], "%)", sep = ""),
              ylab = paste("Comp. 2 (", percent[2], "%)", sep = ""),
              zlab = paste("Comp. 3 (", percent[3], "%)", sep = ""),
              color = pal[as.integer(dds$Tissue)],
              pch = c(17:19)[as.integer(dds$Flexure)])
legend("topleft",
       fill = pal[1:nlevels(dds$Tissue)],
       legend = levels(dds$Tissue))

legend("topright",
       pch = 17:19,
       legend = levels(dds$Flexure))

par(mar = mar)

#' ### 2D
pc.dat <- bind_cols(PC1 = pc$x[, 1],
                    PC2 = pc$x[, 2],
                    samples)

p <- ggplot(pc.dat, aes(x = PC1, y = PC2, col = Tissue, shape = Flexure, text = Name)) + 
  geom_point(size = 2) + 
  ggtitle("Principal Component Analysis", subtitle = "variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis = list(title = paste("PC1 (", percent[1], "%)", sep = "")),
         yaxis = list(title = paste("PC2 (", percent[2], "%)", sep = "")))
# add replicate number, plot number instead of shapes

#' ### Heatmap
#' 
#' Filter for noise
#' 
conds <- factor(paste(samples$Flexure))
sels <- rangeFeatureSelect(counts = vst,
                           conditions = conds,
                           nrep = 5)
vst.cutoff <- 2

#' * Heatmap of "all" genes
#' Not all genes would fit in the plot so this plot only contains a small part of the entire data.
#' (Added different margins and turned the Col labels 45 degrees)
par(mar = c(7, 4, 4, 2) + 0.1)
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff + 1]], ]))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X, method = "ward.D2")},
          labRow = NA, trace = "none",
          labCol = conds,
          margins = c(12, 8),
          srtCol = 45,
          col = hpal)

plot(as.hclust(hm$colDendrogram), xlab = "", sub = "")


#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
