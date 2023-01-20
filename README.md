# Flex wood

## Setup
```{bash}
ln -s /mnt/picea/projects/aspseq/emellerowicz/t89-rot-sta data
```

## Folders
doc/              - Contains the metadata of the samples and the ribbons used in pheatmap plots.
MultiQC/          - Contains results of the MultiQC analysis.
pipeline/         - Contains bash scripts used for the MultiQC analysis and salmon alignment.
src/R/            - Contains the QA and DE scripts used to generate the results. The DE script is an expanded version of the DE                              script in the results folder and includes the analysis of the DEGs within the AspWood and NorWood datasets to                            generate the results in their respective results folders.
UPSCb-common/     - Is linked to a repository containing source code used within UPSC. 
    
