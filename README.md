# Co-transcriptome analysis of plants and commensal bacteria
This repository contains some key data and scripts used in Nobori et al., 2020, bioRxiv.

## **Scripts**
**drawing_ko_expression_plots.R**
This script is for generating boxplots of commensal genes with various functional (KEGG) annotations shown in Fig. 2, 4, S5, S6, S9, S10, S11, and S12.

**plant-associated_genes.R**
This script is for Fig. 3B and Fig. 4.
For each commensal strain, expression fold changes (vs in vitro) of plant-associated (PA) genes or non-PA genes are plotted (Fig. 3B). Statistical tests are performed for the enrichment of genes that are both PA and up-regulated or nonPA and down-regulated (Fig. 3B). KEGG enrichment analysis is performed for genes that are PA/nonPA and up/down-regulated (used for Fig. 4).

**plant_bacteria_integrated_correlation.R**
This script is to generate the correlation heatmap in Fig. 6.

## **Data**
