# Co-transcriptome analysis of plants and commensal bacteria
This repository contains some key data and scripts used in Nobori et al., 2021, bioRxiv.

## **Scripts**
**[drawing_ko_expression_plots.R](scripts/drawing_ko_expression_plots.R)**\
This script is for generating boxplots of commensal genes for various functions (KEGG annotations) (Fig. 2, 4, S5, S6, S9, S10, S11, and S12).

**[plant-associated_genes.R](scripts/plant-associated_genes.R)**\
This script is for Fig. 3B and Fig. 4.
For each commensal strain, expression fold changes (vs in vitro) of plant-associated (PA) genes or non-PA genes are plotted (Fig. 3B). Statistical tests are performed for the enrichment of genes that are both PA and up-regulated or nonPA and down-regulated (Fig. 3B). KEGG enrichment analysis is performed for genes that are PA/nonPA and up/down-regulated (used for Fig. 4).

**[plant_bacteria_integrated_correlation.R](scripts/plant_bacteria_integrated_correlation.R)**\
This script is to generate the correlation heatmap in Fig. 6.

## **Data**
**[bacterial_genome_files](data/bacterial_genome_files/), [bacterial_gff_files](data/bacterial_gff_files/)**\
The data were used for mapping/counting RNA-seq reads\ 
Obtained from:\
https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA297956 \
https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA297942 \
https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA298127

**[RNA-seq count files](data/count_files/)**\
Raw count files. To integrate bacterial transcriptome data, the genes of each strain should be annotated with orthologous groups and/or KEGG orthology available [here](data/commensal_annotation_files/).

**[commensal_annotation_files](data/commensal_annotation_files/)**\
The files summarize gene IDs, Orthologous Groups (OGs), KEGG Orthology (KOs).

**[normalized_bacterial_RNA-seq_data](data/normalized_bacterial_RNA-seq_data/)**\
Normalized RNA-seq data for each strain with differentially expressed gene (DEG) information. Data of different bacterial strains can be integrated using OGs (an integrated data is in [processed_data](data/processed_data/)).

**[ko_database_for_each_strain](data/ko_database_for_each_strain/)**\
KEGG annotations for the genome of each strain. These data are used for KO enrichment analyses.

**[plant-associated_genes](data/plant-associated_genes/)**\
Statistical summaries of a comparative genomics study for plant-associated genes (obtained from Levy et al., 2017, Nature Genetics).

**[processed_data](data/processed_data/)**\
bacteria_RNA-seq_combined.txt: bacterial RNA-seq data combined using OGs.\
plant_RNA-seq_fitted_mean_DEG_atleast_one_strain.txt: expression of plant genes that are DEG in at least one condition\
These data were used for the analysis in Fig. 6C.




