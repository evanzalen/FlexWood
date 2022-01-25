#' ---
#' title: "Differential Expression of Flexwood data"
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
    library(data.table)
    library(DESeq2)
    library(dplyr)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(igraph)
    library(magrittr)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})


dfnew <- df[,1:604]

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal = brewer.pal(8, "Dark2")
hpal <- colorRampPalette(c("blue", "white", "red"))(100)
mar <- par("mar")

#' TODO remember to add a function to check for correlation between logFC and 
#' transcript length - check for the effective length

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values across your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds = dds, vst = vst, gene_id = gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=MDay,y=value,col=MGenotype,group=MGenotype)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              annot = tibble(),
                              labels=dds$Name,
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter)
    } else {
        res <- results(dds,contrast=contrast,filter = filter)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds, normalized = TRUE), probs = max.theta), digits = 5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
        message(sprintf("There are %s genes that are DE with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s",
                        sum(sel),padj,
                        lfc,expression_cutoff))
    }
    
    res %<>% as.data.frame() %>% rownames_to_column("ID")
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel])==0
        if (sum(val) >0){
            warning(sprintf("There are %s DE genes that have no vst expression in the selected samples",sum(val)))
            sel[sel][val] <- FALSE
        }    
        
        if(export){
          if(nrow(annot) > 0) {
          res %<>% left_join(annot, by = "ID")}
          if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write_csv(res,file = file.path(default_dir,paste0(default_prefix,"results.csv")))
            write_csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=res[sel, "ID"],
                up=res[sel & res$log2FoldChange > 0, "ID"],
                dn=res[sel & res$log2FoldChange < 0, "ID"]))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))
    
    # sanity
    if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
        message(paste("No enrichment for",task))
    } else {
    
        # write out
        if(export){
            write_tsv(enrichment[[task]],
                      path=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    path=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
        
        if(plot){
            sapply(go.namespace,function(ns){
                titles <- c(BP="Biological Process",
                            CC="Cellular Component",
                            MF="Molecular Function")
                suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                               namespace=ns,
                                                               de=de,title=paste(default_prefix,titles[ns]))},
                                          error = function(e) {
                                              message(paste("Treemap plot failed for",ns, 
                                                            "because of:",e))
                                              return(NULL)
                                          }))
            })
        }
    }
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/salmon/dds-S339.rda"))

#' ## Normalisation for visualisation

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"), showWarnings = FALSE)
save(vst, file = here("data/analysis/DE/vst-aware-S339.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware-S339.tsv"))

#' Getting the annotation for the DE genes itself.
PotraV2Annotation_file <- "ftp://plantgenie.org/Data/PopGenIE/Populus_tremula/v2.2/annotation/blast2go/Potra22_blast2go_GO_export.txt"
Potra_blast2go <- read_tsv(PotraV2Annotation_file, 
                           show_col_types = FALSE) %>%
  filter(grepl("\\.1$", `Sequence Name`)) %>% 
  mutate(ID = gsub("\\.\\d+$", "", `Sequence Name`)) %>% 
  select("ID", matches("Annotation|Enzyme|InterPro"))

#' ## Gene of interests
#' ```{r goi, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you can plot the expression pattern of your gene of
#' interest. You need to have a list of genes in a text file, one geneID per line
#' The ID should exist in your vst data.
#' Note that for the plot to work, you also need to edit the first function (`line_plot`)
#' at the top of this file
#' ```
#goi <- read_lines(here("doc/goi.txt"))
#stopifnot(all(goi %in% rownames(vst)))
#dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' ## Results
#' #' ```{r res, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to define the contrast you want to 
#' study - see the example in the next block. 
#' 
#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' 
#' The `label` argument is typically one (or a combination) of the metadata stored in colData
#' 
#' The function allows for looking at the independent filtering results using `debug=TRUE`
#' 
#' If you are not satisfied with the default from DESeq2, you can set your own cutoff using `expression_cutoff`
#' 
#' You can change the default output file prefix using `default_prefix`
#' 
#' You can select the set of samples to be added to the `heatmap`, using the `sample_sel` argument. It takes a logical vector.
#' 
#' ```

#' ```{r contrast, echo=FALSE,eval=FALSE}
contrast1L <- extract_results(dds = dds, vst = vst, contrast = "BioID_T89_ROT_vs_T89_STA")
contrast2L <- extract_results(dds = dds, vst = vst, contrast = "TW_no_vs_yes")
#' ```

#' ```
#' The tissue type still has a too big effect on the DESeq to distinguish what exactly are the DE genes when only focusing on ROT and STA.
#' Therefore we rerun the dds on the tissue specific data separately to take away the effect the tissue type has.
#' ```
ddsp <- dds[, dds$Tissue == "Phloem"]
ddsx <- dds[, dds$Tissue == "Xylem"]

vsdp <- varianceStabilizingTransformation(ddsp, blind = FALSE)
vstp <- assay(vsdp)
vstp <- vstp - min(vstp)
vsdx <- varianceStabilizingTransformation(ddsx, blind = FALSE)
vstx <- assay(vsdx)
vstx <- vstx - min(vstx)


#' ## Differential Expression
ddsp <- DESeq(ddsp)
ddsx <- DESeq(ddsx)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(ddsp)
plotDispEsts(ddsx)

#' Check the different contrasts
resultsNames(ddsp)
resultsNames(ddsx)

contrast1pL <- extract_results(dds = ddsp, vst = vstp, contrast = "BioID_T89_ROT_vs_T89_STA", annot = Potra_blast2go)
contrast1xL <- extract_results(dds = ddsx, vst = vstx, contrast = "BioID_T89_ROT_vs_T89_STA", annot = Potra_blast2go)

contrast2pL <- extract_results(dds = ddsp, vst = vstp, contrast = "TW_no_vs_yes", annot = Potra_blast2go)
contrast2xL <- extract_results(dds = ddsx, vst = vstx, contrast = "TW_no_vs_yes", annot = Potra_blast2go)

#' ### Venn Diagram
#' ```{r venn, echo = FALSE, eval = FALSE}
#' CHANGEME - Here, you typically would have run several contrasts and you want
#' to assess their overlap plotting VennDiagrams.
#' 
#' In the examples below, we assume that these results have been saved in a list
#' called `res.list`
#' ```
res.list <- list(c(contrast1p, contrast1x))

#' ```{r venn2, echo=FALSE,eval=FALSE}
grid.newpage()
grid.draw(draw.pairwise.venn(area1 = length(contrast1pL$all), area2 = length(contrast2pL$all),
                             cross.area = length(intersect(contrast1pL$all, contrast2pL$all)), category = c("BioID", "TW"),
                             fill = c("#6baed6", "#fd8d3c"), cat.pos = c(225, 180), euler.d = TRUE, sep.dist = 0.09,
                             rotation.degree = 45, filename = NULL
))

grid.newpage()
grid.draw(draw.pairwise.venn(area1 = length(contrast1xL$all), area2 = length(contrast2xL$all),
                             cross.area = length(intersect(contrast2pL$all, contrast2xL$all)), category = c("BioID", "TW"),
                             fill = c("#6baed6", "#fd8d3c"), cat.pos = c(225, 180), euler.d = TRUE, sep.dist = 0.09,
                             rotation.degree = 45, filename = NULL
))
#' ```

#' ## Analysis in aspect of a Network
#' ```
#' Now we load in the Aspwood network to see where our DE genes are located and what genes are connected to them.
#' 
#' 
#' ```
aspwood_file <- "ftp://anonymous@plantgenie.org/Data/PopGenIE/Populus_tremula/v2.2/Expression/Network/network_seidr.txt"
aspwood <- read_tsv(aspwood_file,
              col_names = c("dataset", "Source", "Target", "direction", "x1", "x2", "irp"),
               col_types = NULL, show_col_types = FALSE)
as.factor(aspwood$dataset)

#' Expression values
aspwood_expr_file <- "ftp://anonymous@plantgenie.org/Data/PopGenIE/Populus_tremula/v2.2/Expression/AspWood_tpm.txt"
aspwood_expr <- read_tsv(aspwood_expr_file,
                      col_names = c("gene", "sample", "tpm"),
                      col_types = "ccd") %>% 
  separate(aspwood_expr, col = c("gene", "sample", "tpm"), sep = "\t", remove = FALSE, convert = TRUE)

#' Lets look at the genes of interest (GOI) a little closer.
vst_PBID <- vstp[rownames(vstp) %in% contrast1pL$all,]
vst_PLig <- vstp[rownames(vstp) %in% contrast2pL$all,]
vst_XBID <- vstx[rownames(vstx) %in% contrast1xL$all,]
vst_XLig <- vstx[rownames(vstx) %in% contrast2xL$all,]

condsP <- factor(paste(ddsp$BioID, ddsp$TW))
condsX <- factor(paste(ddsx$Name, ddsx$TW))

selsp <- rangeFeatureSelect(counts = vstp[rownames(vstp) %in% PuniqueBID,],
                           conditions = condsP,
                           nrep = 5)
selspL <- rangeFeatureSelect(counts = vstp[rownames(vstp) %in% PuniqueLig,],
                            conditions = condsP,
                            nrep = 5)

selsx <- rangeFeatureSelect(counts = vstx[rownames(vstx) %in% XuniqueBID,],
                            conditions = condsX,
                            nrep = 5)
selsxL <- rangeFeatureSelect(counts = vstx[rownames(vstx) %in% XuniqueLig,],
                            conditions = condsX,
                            nrep = 5)

vst.cutoff <- 2

par(mar=c(7, 4, 4, 2)+0.1)
hm_PBID <- heatmap.2(t(scale(t(vst_PBID[selsp[[vst.cutoff + 1]], ]))),
                 distfun = pearson.dist,
                 hclustfun = function(X){hclust(X, method = "ward.D2")},
                 labRow = NA, trace = "none",
                 labCol = condsP,
                 margins = c(12, 8),
                 srtCol = 45,
                 col = hpal)
hmPLig <- heatmap.2(t(scale(t(vst_PLig[selspL[[vst.cutoff + 1]], ]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X, method = "ward.D2")},
                labRow = NA, trace = "none",
                labCol = condsP,
                margins = c(12, 8),
                srtCol = 45,
                col = hpal)

hmXBID <- heatmap.2(t(scale(t(vst_XBID[selsx[[vst.cutoff + 1]], ]))),
                 distfun = pearson.dist,
                 hclustfun = function(X){hclust(X, method = "ward.D2")},
                 labRow = NA, trace = "none",
                 labCol = condsX,
                 margins = c(12, 8),
                 srtCol = 45,
                 col = hpal)
hmXLig <- heatmap.2(t(scale(t(vst_XLig[selsxL[[vst.cutoff + 1]], ]))),
                    distfun = pearson.dist,
                    hclustfun = function(X){hclust(X, method = "ward.D2")},
                    labRow = NA, trace = "none",
                    labCol = condsX,
                    margins = c(12, 8),
                    srtCol = 45,
                    col = hpal)

#' Are the genes within the Aspwood network?
GOI <- commonDE %in% aspwood$gene

table(PuniqueBID %in% bg)
table(PuniqueLig %in% bg)
table(XuniqueBID %in% bg)
table(XuniqueLig %in% bg)

#' TODO
#' - get FDN of 32 genes - need edgelist for this
#' - create heatmap of common and 2 unique genesets
#' - enrich the genes  

#' ## Enriching the genes
#' Setting up Gofer
source("~/Git/UPSCb-common/src/R/gopher.R")
suppressPackageStartupMessages(library(jsonlite))
bg <- rownames(vst)

#' Running the enrichments
enr_commonDE_P <- gopher(commonDE_P,
                       task = list("go", "kegg", "pfam"),
                       background = bg, url = "potra2")
enr_commonDE_FDN <- gopher(commonDE_FDN,
                       task = list("go", "kegg", "pfam"),
                       background = bg, url = "potra2")

enr_PuniqueBID <- gopher(PuniqueBID,
                    task = list("go"),
                    background = bg, url = "potra2")
enr_PuniqueLig <- gopher(PuniqueLig,
                     task = list("go"),
                     background = bg, url = "potra2")
write.csv(enr_PuniqueBID[["go"]], here("data/analysis/DE/GO_DE_phloem_BioID.csv"))
write.csv(enr_PuniqueLig[["go"]], here("data/analysis/DE/GO_DE_phloem_Lignin.csv"))

enr_XuniqueBID <- gopher(XuniqueBID,
                     task = list("go"),
                     background = bg, url = "potra2")
enr_XuniqueLig <- gopher(XuniqueLig,
                     task = list("go"),
                     background = bg, url = "potra2")
write.csv(enr_XuniqueBID[["go"]], here("data/analysis/DE/GO_DE_xylem_BioID.csv"))
write.csv(enr_XuniqueLig[["go"]], here("data/analysis/DE/GO_DE_xylem_Lignin.csv"))


  
#' # First Degree Neighbours of Genes of Interest
#' Returns a dataframe of 1st degree neighbours
#' from an edge list and a collection of genes.
#' It assumes an undirected network.
#'
#' @param edgeList data frame with at least 2 columns
#' 1 with source genes and 1 with target genes. the order
#' is not important
#' @param genes a vector of the genes of interest to look for
#'
#' @return data frame with 2 columns, FDN (first degree neighbours)
#' and GOI (gene of interest)
#'
#' @examples t <- getFDN(edgeList, "MA_104203g0010")
#' write.table(t, file="myFile.tsv", sep='\t', row.names = F, quote = F)
getFDN <- function(edgelist, genes)
{
  
  res <- lapply(genes, function(gene){
    s2t <- edgelist[edgelist[1] == gene,][2]
    t2s <- edgelist[edgelist[2] == gene,][1]
    union(s2t,t2s)
  })
  names(res) <- genes
  res <- setNames(unlist(res, use.names=F),rep(names(res), lengths(res)))
  myRes <- data.frame(FDN = res, GOI = names(res))
}
#FDNtest <- getFDN(edgelist[,1:2], contrast1pL$all)

#' Phloem
edgelist <- read_tsv("/mnt/picea/projects/aspseq/nstreet/v2/network/aspwood/EdgeList.txt", 
                     col_names = c("Source", "Target", "Direction", "irp_score"))[,-3]

sourcephloem <- as.data.frame(edgelist[edgelist$Source %in% contrast1pL$all,])
targetphloem <- as.data.frame(edgelist[edgelist$Target %in% contrast1pL$all,][,c(2:1,3)])
names(targetphloem) <- c("Source", "Target", "irp_score")

FDNphloem_s339 <- bind_rows(sourcephloem, targetphloem)
write.csv(FDNphloem_s339, here("data/analysis/DE/FDNphloem-s339.csv"))

#' Xylem
sourcexylem <- as.data.frame(edgelist[edgelist$Source %in% contrast1xL$all,])
targetxylem <- as.data.frame(edgelist[edgelist$Target %in% contrast1xL$all,][,c(2:1,3)])
names(targetxylem) <- c("Source", "Target", "irp_score")

FDNxylem_s339 <- (bind_rows(sourcexylem, targetxylem))
FDNxylem_s339 <- distinct(FDNxylem_s339[,1:2], .keep_all = FALSE)
write.csv(FDNxylem_s339, here("data/analysis/DE/FDNxylem-s339.csv"))


grafxylem <- graph.edgelist(as.matrix(FDNxylem_s339[,1:2]))
clustxylem <- clusters(grafxylem)

names(FDNxylem_s339) <- c("GOI", "FDN", "IRP")

boolmatrixylem <- sapply(split(FDNxylem_s339$FDN,FDNxylem_s339$GOI), function(v){contrast1xL$all %in%v})
boolXylem <- as.matrix(sapply(split(FDNxylem_s339$FDN, FDNxylem_s339$GOI), function(v, goi){goi %in% v}, unique(FDNxylem_s339$GOI)))

HMBOOL <- as.numeric(as.matrix(str(boolXylem)))
hmBool <- heatmap.2(t(HMBOOL),
                    distfun = pearson.dist,
                    hclustfun = function(X){hclust(X, method = "ward.D2")},
                    labRow = NA, trace = "none",
                    margins = c(12, 8),
                    srtCol = 45,
                    col = hpal)
#' ### Gene Ontology enrichment
#' ```{r go, echo=FALSE, eval=FALSE}
#' Once you have obtained a list of candidate genes, you most probably want
#' to annotate them.
#' 
#' In the following example, we first identify the background; _i.e._ the
#' population of expressed genes. We select the genes expressed in a least
#' 2 replicate of one condition at a cutoff of `exp`.
#' 
#' Next we run the enrichment, in the example against `athaliana` using 
#' the gofer3 REST API (interfaced through the gopher.R script loaded at the
#' beginning of this fil).
#' 
#' Finally we export the go enrichment as a complete table.
#' We used to export another table consisting
#' of only the `id` and `padj` columns for using as input for _e.g._
#' REVIGO; but since flash is EOL and REVIGO not updated, we instead rely on 
#' the RtoolBox treemap.
#' 
#' In addition we now also export the list of genes that most likely resulted in
#' the corresponding go enrichment.
#' 
#' Make sure to change the `url` to match your species
#' 
#' ```
background <- rownames(vst)[featureSelect(vst, dds$BioID, exp = CHANGEME)]

enr.list <- lapply(res.list, function(r){
    lapply(r,gopher,background=background,task="go",url="athaliana")
})

dev.null <- lapply(names(enr.list),function(n){
    lapply(names(enr.list[[n]]),function(de){
        extractEnrichmentResults(enr.list[[n]][[de]],
                                 diff.exp=de,
                                 genes=res.list[[n]][[de]],
                                 default_prefix=paste(n,de,sep="-"),
                                 url="athaliana")
    })
})

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


