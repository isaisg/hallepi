# hallepi
Repository associated with the publication

#### The effects of soil phosphorous content on microbiota are driven by the plant phosphate starvation response
Omri M. Finkel, Isai Salas-González, Gabriel Castrillo, Stijn Spaepen, Theresa F. Law, Corbin D. Jones, Jeffery L. Dangl
bioRxiv 608133; doi: https://doi.org/10.1101/608133

If you use have iny inquire about the scripts in this repository contact isai@email.unc.edu or open an issue ticket.

This repository is arranged in 3 folders:
scripts: Contains all the R scripts used for all the analyses presented in the article.<br />
cleandata: Contains RDS objects ready to be read in R. These objects are organized structures of data required to reproduce all the results and figures in the manuscript <br />
rawdata: Contains tab delimited files with adittional data. <br />

### Scripts
We utilized prefixes to separate the scripts in different categories. Additionally we numbered the scripts in each category. To reproduce all analyses for a category please be sure to run the scripts according to the numeric order. <br />
amplicon_bacteria : These scripts analyze the bacterial census data <br />
amplicon_fungi: These scripts analyze the fungal census data <br />
amplicon_phenotypes: These scripts analyze the phenotypic data obtained from the plants grown in the Halle soils <br />
amplicon_syncom: These scripts analyze the data from the Synthetic Community experiments <br />
rnaseq_hallepi: These scripts analyze the transcriptomic data obtained from plants grown in the Halle soils <br />
rnaseq_syncom: These scripts analyze the transcriptomic data from the SynCom agar experiment <br />
amplicon_create_suptable: Scripts to create the supplementary tables that populate the manuscript <br />

### Cleandata
dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS: R object that contains the ASV Table, Metadata and Taxonomic profiles of the bacterial ASVs obtained from the Halle soil experiment. <br />
dat_hallepi_amplicon_fungi_otus_soil.2019.RDS: R object that contains the OTU Table, Metadata and Taxonomic profile  of the fungal OTUs obtained from the Halle soil experiment. <br />
dat_hallepi_syncom_bacteria_useq97.RDS: R object that contains the USeq table for the SynCom experiment <br />
dat_rnaseq_hallepi.RDS: R object that contains the Transcriptomic matrix (gene X samples) and Metadata for the RNA-Seq data derived from the Halle soil experiment. <br />
dat_rnaseq_syncom.RDS: R object that contains the Transcriptomic matrix (gene X sapmles) and Metadata for the RNA-Seq data derived from the SynCom experiment. <br />
dds_deseq2_hallepi_amplicon_bacteria_asv_maineffects.RDS: R object that contains the results of the GLM ran over the bacterial ASVS testing the main effects (Fraction,Genotype and Phosphate) across the Halle soil experiment. <br />
dds_deseq2_hallepi_amplicon_fungi_otus_maineffects.RDS: R object that contains the results of the GLM ran over the fungal OTUs testing the main effects (Fraction,Genotype and Phosphate) across the Halle soil experiment <br />
dds_deseq2_hallepi_amplicon_bacteria_asv_groupsgenotypephosphate.RDS: R object that contains the result of the GLM ran over the bacterial ASVS testing the Genotype:Phosphate interaction across the Halle soil experiment. <br />
dds_deseq2_hallepi_amplicon_fungi_otus_groupsgenotypephosphate.RDS: R object that contains the result of the GLM ran over the fungal OTUs testing the Genotype:Phosphate interaction across the Halle soil experiment <br />
