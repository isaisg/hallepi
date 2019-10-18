library(ohchibi)
library(DESeq2)
library(paletteer)
library(palettesPM)
library(scales)
library(UpSetR)

set.seed(seed = 130816)

#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_asvs_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$ASV

pval_thres <- 0.1

########## Determine ASVs that have genotype inside Pi treatments #######
res_inside_low_col_pp <- results(object = dds,
                                 contrast = c("group","phr1_phl1_low","Col_0_low"))
res_inside_low_col_phf <- results(object = dds,
                                  contrast = c("group","phf1_low","Col_0_low"))

res_inside_low_phf_pp <- results(object = dds,
                                 contrast = c("group","phf1_low","phr1_phl1_low"))


wt <- union(res_inside_low_col_pp %>%
              subset(padj < pval_thres & log2FoldChange < 0) %>%
              rownames,
            res_inside_low_col_phf %>%
              subset(padj < pval_thres & log2FoldChange < 0) %>%
              rownames
)

single <- union(res_inside_low_phf_pp %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames,
                res_inside_low_col_phf %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames
)

double <- union(res_inside_low_col_pp %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames,
                res_inside_low_phf_pp %>%
                  subset(padj < pval_thres & log2FoldChange < 0) %>%
                  rownames
)

lista <- list(wt = wt,single = single,double = double)
pdf(file = "../figures/upsetr_fungi_low_genotype_effect.pdf")
upset(data = fromList(lista))
dev.off()


wt_specific <- which(!(wt %in% c(single,double))) %>% wt[.]
single_specific <-  which(!(single %in% c(wt,double))) %>% single[.]
double_specific <-  which(!(double %in% c(wt,single))) %>% double[.]
int_all <- intersect(wt,single) %>% intersect(double)
int_wt_single <- intersect(wt,single)
int_wt_single <- which(!(int_wt_single %in% double)) %>% int_wt_single[.]
int_wt_double <- intersect(wt,double)
int_wt_double <- which(!(int_wt_double %in% single)) %>% int_wt_double[.]
int_single_double <- intersect(single,double)
int_single_double <- which(!(int_single_double %in% wt)) %>% int_single_double[.]

lista <- list(wt_specific = wt_specific, single_specific = single_specific, double_specific = double_specific,
              int_all = int_all, int_wt_single = int_wt_single, int_wt_double = int_wt_double,
              int_single_double = int_single_double)

saveRDS(object = lista,file = "../cleandata/list_fungi_asv_low_genotypeeffect_root.RDS")




### Add the abundance patterns 
#Retrieve abundance of each asv
Dat_otu <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
melted_ab_all <- Dat_otu$RelativeAbundance$Tab %>% melt
colnames(melted_ab_all) <- c("Id","DADA2_Header","RA")
map_ab <- Dat_otu$RelativeAbundance$Map 
melted_ab_all <- merge(melted_ab_all,map_ab, by = "DADA2_Header")
#Compute average RE in fraction 
melted_ab <- dcast(data = melted_ab_all,formula = Id ~ Fraction,fun.aggregate = mean,value.var = "RA")



#Check the total abundance in root of all those ASVs
melted_ab_sub <- lista %>% unlist %>% unique %>% match(melted_ab$Id) %>%
  melted_ab[.,] %>% droplevels
melted_ab_sub$Root %>% sum


rm(list=ls())
set.seed(seed = 130816)

#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_asvs_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$ASV

pval_thres <- 0.1


########## Determine ASVs that have genotype inside Pi treatments #######
res_inside_low_col_pp <- results(object = dds,
                                 contrast = c("group","phr1_phl1_low_Pi","Col_0_low_Pi"))
res_inside_low_col_phf <- results(object = dds,
                                  contrast = c("group","phf1_low_Pi","Col_0_low_Pi"))

res_inside_low_phf_pp <- results(object = dds,
                                 contrast = c("group","phf1_low_Pi","phr1_phl1_low_Pi"))


wt <- union(res_inside_low_col_pp %>%
              subset(padj < pval_thres & log2FoldChange < 0) %>%
              rownames,
            res_inside_low_col_phf %>%
              subset(padj < pval_thres & log2FoldChange < 0) %>%
              rownames
)

single <- union(res_inside_low_phf_pp %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames,
                res_inside_low_col_phf %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames
)

double <- union(res_inside_low_col_pp %>%
                  subset(padj < pval_thres & log2FoldChange > 0) %>%
                  rownames,
                res_inside_low_phf_pp %>%
                  subset(padj < pval_thres & log2FoldChange < 0) %>%
                  rownames
)

lista <- list(wt = wt,single = single,double = double)
pdf(file = "../figures/upsetr_fungi_low+Pi_genotype_effect.pdf")
upset(data = fromList(lista))
dev.off()


wt_specific <- which(!(wt %in% c(single,double))) %>% wt[.]
single_specific <-  which(!(single %in% c(wt,double))) %>% single[.]
double_specific <-  which(!(double %in% c(wt,single))) %>% double[.]
int_all <- intersect(wt,single) %>% intersect(double)
int_wt_single <- intersect(wt,single)
int_wt_single <- which(!(int_wt_single %in% double)) %>% int_wt_single[.]
int_wt_double <- intersect(wt,double)
int_wt_double <- which(!(int_wt_double %in% single)) %>% int_wt_double[.]
int_single_double <- intersect(single,double)
int_single_double <- which(!(int_single_double %in% wt)) %>% int_single_double[.]

lista <- list(wt_specific = wt_specific, single_specific = single_specific, double_specific = double_specific,
              int_all = int_all, int_wt_single = int_wt_single, int_wt_double = int_wt_double,
              int_single_double = int_single_double)

saveRDS(object = lista,file = "../cleandata/list_fungi_asv_low+Pi_genotypeeffect_root.RDS")




### Add the abundance patterns 
#Retrieve abundance of each asv
Dat_otu <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
melted_ab_all <- Dat_otu$RelativeAbundance$Tab %>% melt
colnames(melted_ab_all) <- c("Id","DADA2_Header","RA")
map_ab <- Dat_otu$RelativeAbundance$Map 
melted_ab_all <- merge(melted_ab_all,map_ab, by = "DADA2_Header")
#Compute average RE in fraction 
melted_ab <- dcast(data = melted_ab_all,formula = Id ~ Fraction,fun.aggregate = mean,value.var = "RA")



#Check the total abundance in root of all those ASVs
melted_ab_sub <- lista %>% unlist %>% unique %>% match(melted_ab$Id) %>%
  melted_ab[.,] %>% droplevels
melted_ab_sub$Root %>% sum
