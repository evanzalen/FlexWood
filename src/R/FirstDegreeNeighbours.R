#' ---
#' title: "Gene Of Interest First Degree Neighbour"
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
  library(here)
  library(igraph)
  library(readr)
  library(RColorBrewer)
  library(VennDiagram)
})

#' * Palette
pal=brewer.pal(8,"Dark2")

#' * Data
aspwoodedge <- read.table("/mnt/picea/projects/aspseq/nstreet/v2/network/aspwood/EdgeList.txt",
                       header=FALSE,as.is=TRUE)[,-3]

edgelist <- read.table("~/Git/FlexWood/data/analysis/DE/FDNxylem-s339.txt",
                       header=TRUE,as.is=TRUE)
edgelist <- read.table("~/Git/FlexWood/data/analysis/DE/FDNphloem-s339.txt",
                       header=TRUE,as.is=TRUE)

graf <- graph.edgelist(as.matrix(edgelist[,1:2]))
clust <- clusters(graf)

comp <- components(graf, "strong")

selectionX <-names(comp$membership)[comp$membership == which(comp$csize == 385)]

strongComp <- list(selectionX)
write.csv(strongComp, "data/analysis/DE/igraph_strong_components.csv")


source(here("UPSCb-common/src/R/gopher.R"))
suppressPackageStartupMessages(library(jsonlite))
bg <- unlist(unique(aspwoodedge[,1:2]))
enrichselec <- gopher(selectionX, background = bg, task = "go", url = "potra2")


#' We can check the structure of the graf
#' 
#' There is one graph with 34,312 genes
clusters(graf)$csize

#' We can look at the neighborhood of the gene of interests
goi <- read_csv(here("data/analysis/DE/unique_phloem_BioID-S339.csv"), show_col_types=FALSE)[,2]

stopifnot(all(goi$x %in% get.vertex.attribute(graf,"name")))
stopifnot(all(goi$x %in% edgelist$V1))

#' Extract the first degree neighbors of these genes from the network
subgraf1 <- make_ego_graph(graf,1,match(goi$ID,get.vertex.attribute(graf,"name")))

#' relatively few genes
barplot(table(sapply(lapply(subgraf1,clusters),"[[","csize")),
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")

sizes <- sapply(lapply(subgraf1,clusters),"[[","csize")
names(sizes) <- goi$Symbol
barplot(sizes,
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")

#' gene overlap
lst <- lapply(subgraf1,get.vertex.attribute,"name")
names(lst) <- goi$Symbol
grid.newpage()
grid.draw(venn.diagram(lst,filename=NULL,fill=pal[1:length(lst)]))
grid.newpage()
grid.draw(venn.diagram(lst,filename=NULL,fill=pal[1:length(lst)],category.names=goi$ID))

#' There are only two genes in common: 
#' 
#' One in between FT1 and GA20ox (Potra2n15c28070)
#' * Potra2n18c32697 which is a myb transcription factor
names(lst) <- goi$ID
intersect(lst[["Potra2n15c28070"]],lst[["Potra2n8c17315"]])

#' One in between AGL8 and GA20ox (Potra2n15c28070, as Potra2n12c24804 has no overlap whatsoever)
#' * Potra2n11c23377 a COP1-interacting protein
intersect(lst[[which(goi$ID=="Potra2n1c2860")]],lst[["FT1"]])

#' These 2 genes connect to 108!
clusters(make_ego_graph(graf,1,
                        get.vertex.attribute(graf,"name") %in% c("Potra2n18c32697","Potra2n11c23377"))[[1]])$csize

#' Check the second degree neighbors
subgraf2 <- make_ego_graph(graf,2,match(goi$ID,get.vertex.attribute(graf,"name")))

#' which explodes the number of genes
sizes <- sapply(lapply(subgraf2,clusters),"[[","csize")
names(sizes) <- goi$Symbol
barplot(sizes,
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")

#' A much larger overlap
lst <- lapply(subgraf2,get.vertex.attribute,"name")
names(lst) <- goi$Symbol
grid.newpage()
grid.draw(venn.diagram(lst,filename=NULL,fill=pal[1:length(lst)]))
grid.newpage()
grid.draw(venn.diagram(lst,filename=NULL,fill=pal[1:length(lst)],category.names=goi$ID))

#' combine all these networks together
fdn1 <- Reduce("%u%",subgraf1)
fdn2 <- Reduce("%u%",subgraf2)

#' Look at how many clusters we get and how many nodes are involved 
#' 1 cluster and 7153 genes, a third of all genes. That's certainly too many.
clusters(fdn2)$csize

#' So, let's look at a static plot of the first degree neighbours
sapply(subgraf1,plot)
plot(fdn1)

#' Let's export the data for visualization
dir.create(here("data/Combined_Data/analysis/fdn"),showWarnings=FALSE)
write_graph(fdn1,format = "graphml",
            file=here("data/Combined_Data/analysis/fdn/firstDegreeNeighbors.graphml"))

#' And as a goi file for looking them up in the cold-warm data
write(get.vertex.attribute(fdn1,"name"),
      file=here("data/Combined_Data/analysis/fdn/firstDegreeNeighbors.txt"))

#' # Session Info
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```

