# Flex wood

## Setup
```{bash}
ln -s /mnt/picea/projects/aspseq/emellerowicz/t89-rot-sta data
```

## Folders
data/             - Contains the raw data and the output of the FastQC, SortMeRNA, Trimmomatic and Salmon alignment.
doc/              - Contains the metadata of the samples and the ribbons used in pheatmap plots.
MultiQC/          - Contains results of the MultiQC analysis.
pipeline/         - Contains bash scripts used for the MultiQC analysis and salmon alignment.
results/          - Contain the lists of Differential Expressed Genes (DEGs) within the Flexwood data per tissue type. 
                    This folder also contain the dds object created for the DE analysis and the R script used to generate the two lists.
results-Aspwood/  - This folder contains the results of the analysis of DEGs within the AspWood dataset. This includes the heatmaps 
                    plotted from AspWood TPM values, each heatmap cluster plotted separately and the cluster specific GO 
                    enrichment results for both tissue types.
results-Norwood/  - This folder contains the same results as in the results-Aspwood folder, but for the Norwood expression values                            instead.
src/R/            - Contains the QA and DE scripts used to generate the final results. The DE script is an expanded version of the DE                        script in the results folder and includes the analysis of the DEGs within the AspWood and NorWood datasets to                            generate the results in their respective results folders.
UPSCb-common/     - Is linked to a repository containing source code used within UPSC. 
    
