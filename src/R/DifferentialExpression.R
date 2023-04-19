#' ---
#' title: "Differential Expression"
#' author: "CHANGEME"
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
    library(gplots)
    library(here)
    library(hyperSpec)
    library(MatrixGenerics)
    library(pheatmap)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
#' ```{r edit1, echo=FALSE,eval=FALSE}
#' CHANGEME - here you need to change the variables in the 
#' plot to display the expression values accross your samples
#' The example below has 2 variables MGenotype and MDay. These 
#' need replacing by the variable(s) of interest in your project
#' ```
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
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
                              export=TRUE,default_dir=here(),
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
                    round(max.theta*100, digits=5),
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
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    file=here(default_dir,
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

#' 4. extract enrichment results from DEGs clusters determined by pheatmap and plotting each cluster in a separate heatmap.
enrich_heatmap_clusters <- function(nclust, datasetname, gene_mat, ann_ribbon, ann_colour, col_gaps, gopher_url, background){
  # plotting a heatmap per cluster with kmeans
  message("Plotting kmeans heatmap")
  kphmap <- pheatmap(gene_mat,
                     main = datasetname,
                     cluster_rows = TRUE, cluster_cols = FALSE,
                     clustering_distance_rows = "correlation", clustering_method = "ward.D2",
                     col = hpal, fontsize_row = 7, fontsize_col = 7, angle_col = 90,
                     show_rownames = TRUE, show_colnames = TRUE,
                     annotation_col = ann_ribbon, annotation_colors = ann_colour,
                     gaps_col = col_gaps, kmeans_k = nclust,
                     filename = here("results-Aspwood/Phloem", paste0("Kmeans_Heatmap_", datasetname, ".png"))
  )
  
  # defining the clusters by the predetermined nr of clusters from a pheatmap from the entire set of DEGs
  ct <- kphmap$kmeans$cluster
  # For Aspwood, write this in a dataframe:
  write_csv(enframe(kphmap$kmeans$cluster, name = "PotraID", value = "Cluster"), file = here("results-Aspwood", paste0("Dataframe_", datasetname, ".csv")))
  # apply the following to each one of these clusters
  return(sapply(1:nclust, 
                function(i, ct, gene_mat, ann_ribbon, ann_colour, col_gaps, gopher_url, background){
                  
                  cl = names(ct)[ct == i]
                  
                  # plotting a heatmap per cluster
                  message(paste0("Plotting heatmap for cluster ", i))
                  phmap <- pheatmap(gene_mat[cl,],
                                    main = paste0("Cluster ", i, " ", datasetname),
                                    cluster_rows = TRUE, cluster_cols = FALSE,
                                    clustering_distance_rows = "correlation", clustering_method = "ward.D2",
                                    col = hpal, fontsize_row = 7, fontsize_col = 7, angle_col = 90,
                                    show_rownames = TRUE, show_colnames = TRUE,
                                    annotation_col = ann_ribbon, annotation_colors = ann_colour,
                                    gaps_col = col_gaps,
                                    filename = here("results-Aspwood/Phloem", paste0("Heatmap_Cluster_", i, "_", datasetname, ".png"))
                  )
                  # enrichment of clusters
                  message(paste0("Enriching cluster ", i))
                  enr <- gopher(cl,
                                task = list("go","kegg", "pfam"),
                                background = background,
                                url = gopher_url)
                  if (length(rownames(enr$go)) > 0) {
                    enr$go <- enr$go[,c("name","padj")] %>% mutate(type = "GO")
                    
                  }
                  if (length(rownames(enr$kegg)) > 0) {
                    enr$kegg <- enr$kegg[,c("name","padj")] %>% mutate(type = "KEGG")
                    
                  }
                  if (length(rownames(enr$pfam)) > 0) {
                    enr$pfam <- enr$pfam[,c("name","padj")] %>% mutate(type = "PFAM")
                    enr$pfam$name <- gsub("\\.", " ", enr$pfam$name)
                    
                  }
                  message("Writing enrichments into file")
                  write_csv(reduce(enr, bind_rows), file = here("results-Aspwood/Phloem", paste0("Enrichment_Cluster_", i, "_", datasetname, ".csv")))
                  
                }, ct, gene_mat, ann_ribbon, ann_colour, col_gaps, gopher_url, background) 
  )}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/salmon/dds.rda"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst, file = here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' Getting the annotation for the DE genes.
PotraV2Annotation_file <- "ftp://plantgenie.org:980/Data/PopGenIE/Populus_tremula/v2.2/annotation/blast2go/Potra22_blast2go_GO_export.txt"
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

#' The `contrast` can be given
#' by name, as a list (numerator/denominator) or as a vector of weight (e.g. c(0,1));
#' read the DESeq2 vignette for more info
#' Getting the contrast results per tissue type and annotating each gene with our blast2go file.
#' The results are automatically stored in a "DE-genes.csv" file in the default directory here("data/analysis/DE")
#' specified in the second function coded at the beginning of this script.
contrast1pL <- extract_results(dds = ddsp, vst = vstp, contrast = "Flexure_T89_ROT_vs_T89_STA", 
                               annot = Potra_blast2go, default_dir = here(), default_prefix = "data/analysis/DE/export/phloem_DEGs_")
contrast1xL <- extract_results(dds = ddsx, vst = vstx, contrast = "Flexure_T89_ROT_vs_T89_STA", 
                               annot = Potra_blast2go, default_dir = here(), default_prefix = "data/analysis/DE/export/xylem_DEGs_")

#' ### Venn Diagram
#' Assessing overlap of DEGs between tissues.
#' 32 genes overlap between the two different tissue types.
grid.newpage()
grid.draw(draw.pairwise.venn(area1 = length(contrast1pL$all), area2 = length(contrast1xL$all),
                             cross.area = length(intersect(contrast1pL$all, contrast1xL$all)), category = c("Phloem", "Xylem"),
                             fill = c("#6baed6", "#fd8d3c"), cat.pos = c(110, 310), euler.d = TRUE, sep.dist = 0.09,
                             rotation.degree = 45, filename = NULL
))
#grid.draw(venn.diagram(list(phloem = contrast1pL$all, xylem = contrast1xL$all), filename = NULL, fill = pal[1:2]))


#' ## Analysis of tissue type DEGs in Aspwood network.
#' First we enrich each set of DEGs. Then we extract the expression values for these two lists of genes from the Aspwood network
#' and plot the heatmaps. Thereafter we cluster the genes based on the heatmap expression pattern and enrich each cluster.
#' Then we do the same for the Norwood orthologues of the genes.
#' ```
#' Setting up Gopher and using the Aspwood network as background (bg). 
#' Loading the bg might take a while because it's a big file.
#' Running Gopher takes a minute or 2 as well, depending on how big the DEG list is.
#' ```
#' Enriching the tissue specific DEGs.
aspwood_expr <- read_delim("ftp://anonymous@plantgenie.org:980/Data/PopGenIE/Populus_tremula/v2.2/Expression/AspWood_tpm.txt", 
                           delim = "\t")

#' Visualising expression values of tissue specific DEGs using the expression values from the Aspwood network in a heatmap.
#' We transform the Aspwood expression file into a matrix with genes as rownames and the samples as column names, 
#' using the tpm as values. The original K9 samples were removed from the matrix to keep consistency with publications.
#' K5 was renamed to K4
aspwood_matrix <- as.matrix(pivot_wider(aspwood_expr,names_from = sample_name, values_from = expression) %>% 
                              column_to_rownames(var = "gene_id"))
colnames(aspwood_matrix) <- gsub("-.*-", "-", colnames(aspwood_matrix))

sds <- rowSds(log1p(aspwood_matrix))
plot(density(sds))
abline(v=0.15,lty=2,col="grey")

bg <- rownames(aspwood_matrix)[sds>0.15]

#' Extracting the genes from the DEGs list that have expression in Aspwood.
DEGphloem_Aspwood <- aspwood_matrix[bg[bg %in% contrast1pL$all],]
# 201 out of 211 DE phloem genes have expression in Aspwood
DEGxylem_Aspwood <- aspwood_matrix[bg[bg %in% contrast1xL$all],] 
# 348 out of 355 DE xylem genes have expression in Aspwood

#See which DE are in flexwood but not in Aspwood
difDEp <- setdiff(contrast1pL$all, rownames(DEGphloem_Aspwood))
difDEx <- setdiff(contrast1xL$all, rownames(DEGxylem_Aspwood))

#' Running enrichment on the entire lists per tissue.
#' First on the DE from the FlexWood data with Aspwood as a background.
enr_ph <- gopher(rownames(contrast1pL),
                 task = list("go", "kegg", "pfam"),
                 background = bg, url = "potra2")
enr_ph[["go"]]

enr_xy <- gopher(rownames(contrast1xL),
                 task = list("go", "kegg", "pfam"),
                 background = bg, url = "potra2")
enr_xy[["go"]]

#' Second on the DE genes from FlexWood found in the Aspwood network.
enr_phloem <- gopher(rownames(DEGphloem_Aspwood),
                     task = list("go", "kegg", "pfam"),
                     background = bg, url = "potra2")
enr_phloem[["go"]]

enr_xylem <- gopher(rownames(DEGxylem_Aspwood),
                    task = list("go", "kegg", "pfam"),
                    background = bg, url = "potra2")
enr_xylem[["go"]]


#' Getting the labeling of Aspwood cell type to create a ribbon for the heatmap. 
#' Obtained from the supplementary figures from the Aspwood paper.
ribbon <- read.table(here("doc/Aspwood_ribbon.txt"),
                     header = FALSE,
                     sep = "\t",
                     col.names = c("","CellType"),
                     row.names = 1)
ann_colours <- list(
  CellType = c(Phloem = "#99CC33", Cambium = "#FF9900", "Expanding xylem" = "#3399CC", "Lignified xylem" = "#993300"))

#'Plotting the obtained expression values in a heatmap.
genelistscalep <- t(scale(t(DEGphloem_Aspwood)))
genelistscalep[genelistscalep > 2] <- 2
genelistscalep[genelistscalep < -2] <- -2

phm_p <- pheatmap(genelistscalep,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  clustering_distance_rows = "correlation",
                  clustering_method = "ward.D2",
                  col = hpal,
                  fontsize_col = 7,
                  angle_col = 90,
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  annotation_col = ribbon,
                  annotation_colors = ann_colours,
                  cutree_rows = 10,
                  gaps_col = c(25,51,79))

genelistscalex <- t(scale(t(DEGxylem_Aspwood)))
genelistscalex[genelistscalex > 2] <- 2
genelistscalex[genelistscalex < -2] <- -2
phm_x <- pheatmap(genelistscalex,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  clustering_distance_rows = "correlation",
                  clustering_method = "ward.D2",
                  col = hpal,
                  fontsize_col = 7,
                  angle_col = 90,
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  annotation_col = ribbon,
                  annotation_colors = ann_colours,
                  cutree_rows = 10,
                  gaps_col = c(25,51,79))



#' Taking clusters from the heatmap with the cutree function to enrich each cluster of genes individually.

#' Running the enrichments per heatmap cluster.
#' The function can be found at the top of this script listed as function nr. 4. We use the results and inputs from the preliminary 
#' heatmap as arguments for this function. Any time this function is run the cluster nr is assigned randomly, 
#' so it might be that you will obtain different names for the clusters. This is just how the cutree command works. 
#' The only way we found around this is to plot the kmeans heatmap and use those cluster names to directly enrich and plot 
#' the individual clusters. All results are stored in one folder under "results-Aspwood" for phloem and for xylem respectively.
Enr_phloem_clusters <- enrich_heatmap_clusters(nclust = 9, datasetname = "Aspwood_Phloem", gene_mat = genelistscalep, 
                                                    ann_ribbon = ribbon, ann_colour = ann_colours,
                                                    col_gaps = c(25,51,79), gopher_url = "potra2", background = bg)

#' Luckily for xylem we did not need to do an extra split and we could run everything at once without having to change the number of clusters.
Enr_xylem_clusters <- enrich_heatmap_clusters(nclust = 10, datasetname = "Aspwood_Xylem", gene_mat = genelistscalex, 
                                              ann_ribbon = ribbon, ann_colour = ann_colours,
                                              col_gaps = c(25,51,79), gopher_url = "potra2", background = bg)



#' ## Analysis of DEGs in Picea abies 
#' Obtaining Picea abies orthologues for the DEGs with the best Diamond hits.
potra_piabi_diamond_file <- "ftp://anonymous@plantgenie.org:980/Data/Cross-Species/Orthologs/BEST_DIAMOND/potra_piabi_best_diamond.tsv"
potra_piabi_diamond <- read_tsv(potra_piabi_diamond_file,
                                col_names = c("Potra", "Piabi"))

table(potra_piabi_diamond$Potra %in% row.names(DEGphloem_Aspwood)) 
#182 out of 201 genes have an orthologue for the phloem DEGs in the AspWood dataset
table(potra_piabi_diamond$Potra %in% row.names(DEGxylem_Aspwood)) 
#315 out of 348 genes have an orthologue for the xylem DEGs in the AspWood dataset

DEGphloem_spruceOrtho <- potra_piabi_diamond$Piabi[potra_piabi_diamond$Potra %in% row.names(DEGphloem_Aspwood)]
DEGxylem_spruceOrtho <- potra_piabi_diamond$Piabi[potra_piabi_diamond$Potra %in% row.names(DEGxylem_Aspwood)]

#' Getting Norwood expression values for the selected gene lists. We removed the Late Wood samples because there is no 
#' equivalent to it in Populus tremula.
norwood_expr <- read_delim("ftp://anonymous@plantgenie.org:980/Data/PlantGenIE/Picea_abies/v1.0/annotation/expression/expression_norwood_vst.txt",
                           delim = ",", col_names = c("id","sample", "log2"))[,1:3]

norwood_matrix <- pivot_wider(norwood_expr, names_from = sample, values_from = log2) %>% 
                              column_to_rownames(var = "id") %>%
                              dplyr::select(-c("T1-20", "T2-18", "T3-20"))
norwood_matrix <- as.matrix(norwood_matrix[, order(colnames(norwood_matrix))])

#' Filtering the background 
sdsNW <- rowSds(log1p(norwood_matrix))
plot(density(sdsNW))
abline(v=0.21,lty=2,col="grey")

bgNW <- rownames(norwood_matrix)[sdsNW > 0.21]

#' Calculating the heatmaps with the Norwood expression values.
#' For phloem only 44 out of the 182 spruce orthologs are in the Norwood network, 
#' for xylem this was 90 out of 315 genes. 
DEGphloem_Norwood <- norwood_matrix[bgNW[bgNW %in% DEGphloem_spruceOrtho],]
# we have 44
DEGxylem_Norwood <- norwood_matrix[bgNW[bgNW %in% DEGxylem_spruceOrtho],]
# we have 90

#' The code above that enriches each cluster individually, also outputs a dataframe with the PotraIDs and their respective cluster number. 
#' We will use this dataframe to create an overview with the Potra DEGs, their respective heatmap clusters, 
#' their Picea abies orthologue and whether or not this orthologue is expressed in Norwood. 
Aspwood_dataframe_phloem <- read_csv(here("results-Aspwood/Phloem/Dataframe_Aspwood_Phloem.csv")) %>% 
  arrange(Cluster) %>% 
  left_join(potra_piabi_diamond, by = c("PotraID" = "Potra")) %>% 
  mutate("ExpressionNorwood" = Piabi %in% bgNW) %>% 
  write_csv(file = here("results-Aspwood/Phloem/Dataframe_Aspwood_Phloem.csv"))

Aspwood_dataframe_xylem <- read_csv(here("results-Aspwood/Xylem/Dataframe_Aspwood_Xylem.csv")) %>% 
  arrange(Cluster) %>% 
  left_join(potra_piabi_diamond, by = c("PotraID" = "Potra")) %>% 
  mutate("ExpressionNorwood" = Piabi %in% bgNW) %>% 
  write_csv(file = here("results-Aspwood/Xylem/Dataframe_Aspwood_Xylem.csv"))

#' ## Analysis of DE Norwood orthologues
#' Setting up the ribbon for the Norwood heatmaps.
#' SCW = Secondary Cell Wall
#' PCD = Programmed Cell Death
ribbonNW <- read.table(here("doc/Norwood_ribbon.txt"),
                       header = FALSE,
                       sep = "\t",
                       col.names = c("","CellType"),
                       row.names = 1)
ann_coloursNW <- list(
  CellType = c("Cambium" = "#3399CC", "Expanding Xylem"= "#993300", "SCW" = "#993300", "PCD" = "#996600"))

#' The samples in Norwood don't look too well structured, so we decided on a cut-off of 8 clusters with cutree_rows = 8.
#' This number can be changed to adjust the clusters to liking. 
#' The code that follows after can be run without any changes to adapt for the different clustering.
genelistscale.nwp <- t(scale(t(DEGphloem_Norwood)))
genelistscale.nwp[genelistscale.nwp > 2] <- 2
genelistscale.nwp[genelistscale.nwp < -2] <- -2
phm_nwp <- pheatmap(genelistscale.nwp,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    clustering_distance_rows = "correlation",
                    clustering_method = "ward.D2",
                    col = hpal,
                    fontsize_col = 7,
                    angle_col = 90,
                    show_rownames = FALSE,
                    show_colnames = TRUE,
                    annotation_col = ribbonNW,
                    annotation_colors = ann_coloursNW,
                    cutree_rows = 8,
                    gaps_col = c(17,32))

genelistscale.nwx <- t(scale(t(DEGxylem_Norwood)))
genelistscale.nwx[genelistscale.nwx > 2] <- 2
genelistscale.nwx[genelistscale.nwx < -2] <- -2
phm_nwx <- pheatmap(genelistscale.nwx,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    clustering_distance_rows = "correlation",
                    clustering_method = "ward.D2",
                    col = hpal,
                    fontsize_col = 7,
                    angle_col = 90,
                    show_rownames = FALSE,
                    show_colnames = TRUE,
                    annotation_col = ribbonNW,
                    annotation_colors = ann_coloursNW,
                    cutree_rows = 8,
                    gaps_col = c(17,32),filename = "results/Heatmap_Norwood_Xylem.png")


#' Enrichment of the heatmap clusters
Enr_phloem_clustersNW <-enrich_heatmap_clusters(nclust = 8, datasetname = "Norwood_Phloem", gene_mat = genelistscale.nwp, 
                                                ann_ribbon = ribbonNW, ann_colour = ann_coloursNW, 
                                                col_gaps = c(17,32), gopher_url = "pabies", background = bgNW)

Enr_xylem_clustersNW <-enrich_heatmap_clusters(nclust = 8, datasetname = "Norwood_Xylem", gene_mat = genelistscale.nwx, 
                                               ann_ribbon = ribbonNW, ann_colour = ann_coloursNW, 
                                               col_gaps = c(17,32), gopher_url = "pabies", background = bgNW)



#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```