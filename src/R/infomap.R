#' ---
#' title: ""
#' author: "Elena van Zalen"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Setup
#' ## Environment
#' Set the working dir
setwd("~/Git/FlexWood/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/Flexwood")
#' ```

#' ## Libraries
suppressPackageStartupMessages({
  library(igraph)
  library(RColorBrewer)
  library(readr)
  library(here)
  library(gplots)
  library(data.table)
  library(devtools)
  library(tidyverse)
  library(parallel)
  library(scatterplot3d)
  library(plot3D)
  library(vsn)
  library(hyperSpec)
  library(grid)
  library(VennDiagram)
  library(treemap)
  library(pander)
  library(writexl)
})

#' Functions
plotEigengene <- function(data, gene, condition, time,timeUnits="Time",
                          inverse = F, title=T){
  require(ggplot2)
  require(RColorBrewer)
  require(reshape2)
  
  d <- data[,which(colnames(data) %in% gene)]
  pca <- prcomp(scale(d))
  pc1 <- pca$x[, 1]
  if(sum(sign(cor(pca$x[,1,drop = FALSE], d))) < 0) {
    pc1 <- pc1 * -1
  }
  if(inverse)
    pc1 <- pc1 * -1
  pc1i<- pc1 * -1
  
  if(title)
    header=paste("Eigengene plot for gene. Genes:", gene)
  else
    header=gene
  myplot <- ggplot( data.frame(x = time, y = scale(d), g = condition),
                    aes(x = x, y = y, group = g)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", fill = "grey", alpha = 0.75) +
    stat_summary(fun.data = mean_se, geom = "line", aes(col = g), lwd = 1) + #      plot_output_list <- lapply(shiftedColors, function(color) {
    xlab(timeUnits) +
    ylab("Relative expression") + 
    labs(color = "Condition") + 
    ggtitle(header)
  
  return (myplot)
  
}
#' Palette
pal <- brewer.pal(8,"Dark2")

#' # Infomap analysis
#' Read in the original network (not reheadered)
bgr <- read.delim2("/mnt/picea/projects/aspseq/nstreet/v2/network/aspwood/EdgeList.txt",
                   stringsAsFactors = FALSE,
                   header=TRUE,sep="\t")
bg <- unique(unlist(bgr[,2:1])) # background later used for gopher and the directionality in the network

#' Read in the infomap of the network
Network.imap <- read.delim(here("/mnt/picea/projects/aspseq/nstreet/v2/network/aspwood/infomap/resolve.txt"), header = FALSE, sep ="",
                              col.names = c("Path", "lvl1", "lvl2", "lvl3", "dc", "GeneIndex", "Genes"))

## =================== The first level ====================== ##
#' Check the total number of clusters and the total number of genes within these clusters
nr_clus_lvl1 <- unique(Network.imap$lvl1) #2 clusters at lvl1

#' Then we obtain the frequencies (number of genes per cluster)
freq_lvl1 <- as.data.frame(table(Network.imap$lvl1)) #22727 in the first cluster, this is not good
#' and the cluster names as well (at the moment they are numbers)
freq_lvl1$Var1 <- names(table(Network.imap$lvl1))

#' Check the number of cluster in the lowest granularity (the largest clusters/modules)
barplot(table(Network.imap[, 2]), main = "Gene network", xlab = "Main cluster size")

## =================== The second level ====================== ##
Network.imap$lvl1_2 <- paste0(Network.imap$lvl1, ":", Network.imap$lvl2)

nr_clus_lvl2 <- unique(Network.imap$lvl1_2) # 57
freq_lvl2 <- as.data.frame(table(Network.imap$lvl1_2))
freq_lvl2$Var1 <- names(table(Network.imap$lvl1_2))

#' Sum the second column of the first 20 clusters
counts_lvl2 <-freq_lvl2[order(freq_lvl2$Freq, decreasing = T), ]
sum(counts_lvl2[1:20, 2]) #most genes are in the first 20 clusters (22907)

#' The same, but zoomed in the first ~20
par(mar = c(5.1, 3.1, 3.1, 7.1))
pander(table(Network.imap[, 3], xlim = c(0, 20)))
barplot(table(Network.imap[, 3]), xlim = c(0, 20))


## =================== The third level ====================== ##
Network.imap$lvl1_2_3 <- paste0(Network.imap$lvl1, "", Network.imap$lvl2, "", Network.imap$lvl3)

nr_clus_lvl3 <- unique(Network.imap$lvl1_2_3) # 57
freq_lvl3 <- as.data.frame(table(Network.imap$lvl1_2_3))
freq_lvl3$Var1 <- names(table(Network.imap$lvl1_2_3))

#' FDN of GOI
XylemFDN <- list.files(path = here("doc/FirstDegreeNeighbours_All/"), 
                      full.names = TRUE) %>% 
  lapply(read_csv) %>% bind_rows
FDNxylem <- unique(XylemFDN[,1])
FDN_xylem_infomap <- Network.imap[Network.imap[,7] %in% FDNxylem$FDN, ]
barplot(table(FDN_xylem_infomap[,2]), xlim = c(0, 20))


goiPhloem <- read_csv(here("data/analysis/DE/unique_phloem_BioID-S339.csv"), show_col_types=FALSE)
FDN_phloem_infomap <- Network.imap[Network.imap[,7] %in% goiPhloem$x, ]
barplot(table(FDN_phloem_infomap[,3]), xlim = c(0, 20))

FDN_phloem_infomap$Genes[FDN_phloem_infomap$lvl1_2 == 1:2]

### =============================== Enrichments ====================================== ###
#' Obtaining the enrichment for the clusters of interest
source("~/Git/UPSCb-common/src/R/gopher.R")
suppressPackageStartupMessages(library(jsonlite))

#' For the enrichment:
#' Subsetting the infomap dataframe so only the clusters of interest are present
clustvec1 <- coi1$Var1[1:30]
Network_All_coi1 <- NetworkAll.imap[NetworkAll.imap$lvl1 %in% clustvec1, ]
nr_clus_coi1 <- unique(Network_All_coi1$lvl1)

#' Gene Ontology (GO) and Mapmap enrichment
All_Enr_lvl1 <- lapply(unique(Network_All_coi1$lvl1),
                       function(cl){
                         res_lvl1 <- list()
                         enr<-gopher(Network_All_coi1$Genes[Network_All_coi1$lvl1 == cl],
                                     task = list("go","kegg", "pfam"),
                                     background = bg,
                                     url = "potra2")
                         if (length(rownames(enr$go)) > 0) {
                           go <- enr$go[,c("name","padj")]
                           res_lvl1$go <- go
                         }
                         if (length(rownames(enr$kegg)) > 0) {
                           kegg <- enr$kegg[,c("name","padj")]
                           kegg$name <- gsub("\\.", " ", kegg$name)
                           res_lvl1$kegg <- kegg
                         }
                         if (length(rownames(enr$pfam)) > 0) {
                           pfam <- enr$pfam[,c("name","padj")]
                           pfam$name <- gsub("\\.", " ", pfam$name)
                           res_lvl1$pfam <- pfam
                         }
                         res_lvl1
                       }) 

FDN_Enr_cl1 <- gopher(FDN_All_infomap$Genes[FDN_All_infomap$lvl1 == 1],
                      task = list("go", "kegg", "pfam"),
                      background = bg,
                      url = "potra2")
treemap(FDN_Enr_cl1$pfam, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = pal, fontsize.labels=c(15,12), 
        title = "FDN - Level 1 PFAM - Cluster 1", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)
#' This function will write an R object with the Gene Ontology (GO), pfam and kegg
#' enrichments if there are any. If there are no enrichments, it will not write 
#' it into the object.

#' Name the lists with their proper cluster number
names(All_Enr_lvl1) <- clustvec1

### ==================================== Treemaps ==================================== ###
#' Level 1
paltr <- brewer.pal(9,"YlGnBu")
treemap(All_Enr_lvl1$`3`$go, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = paltr, fontsize.labels=c(15,12), 
        title = "AllNetwork - Level 1 GO - Cluster 3", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)

treemap(All_Enr_lvl1$`2`$go, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = paltr, fontsize.labels=c(15,12), 
        title = "AllNetwork - Level 1 GO - Cluster 2", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)
treemap(All_Enr_lvl1$`1`$go, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = paltr, fontsize.labels=c(15,12), 
        title = "AllNetwork - Level 1 GO - Cluster 1", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)
treemap(All_Enr_lvl1$`4`$go, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = paltr, fontsize.labels=c(15,12), 
        title = "AllNetwork - Level 1 GO - Cluster 4", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)
treemap(All_Enr_lvl1$`5`$go, index = "name", vSize = "padj", type = "dens", 
        vColor = "padj", palette = paltr, fontsize.labels=c(15,12), 
        title = "AllNetwork - Level 1 GO - Cluster 5", inflate.labels = FALSE, lowerbound.cex.labels = 0, 
        bg.labels = "NA", position.legend = "none", file = FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
