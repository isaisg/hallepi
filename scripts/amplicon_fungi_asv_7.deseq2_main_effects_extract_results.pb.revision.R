library(ohchibi)
library(DESeq2)

set.seed(seed = 130816)
setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')

##Read
df_tax <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS") %$%
  df_tax
df_tax <- df_tax[,-2]


dds_res <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_asvs_maineffects.RDS")

########################## ASvs ###########################################
### Root vs Soil ###
dds <- dds_res$ASV

res_rs <- dds %>% results(contrast = c("Fraction","Root","Soil")) %>%
  as.data.frame
res_rs$Contrast <- rep("Root-Soil",nrow(res_rs))
res_rs$Id <- rownames(res_rs)

### Shoot vs Soil ###
res_ss <- dds %>% results(contrast = c("Fraction","Shoot","Soil")) %>%
  as.data.frame
res_ss$Contrast <- rep("Shoot-Soil",nrow(res_ss))
res_ss$Id <- rownames(res_ss)

### Shoot vs Root ###
res_sr <- dds %>% results(contrast = c("Fraction","Root","Shoot")) %>%
  as.data.frame

res_sr$Contrast <- rep("Root-Shoot",nrow(res_sr))
res_sr$Id <- rownames(res_sr)


res <- rbind(res_rs,res_ss,res_sr)
rownames(res) <- NULL
merged <- merge(res,df_tax, by = "Id")

#Correc the pvalues globally
merged$padj <- merged$pvalue %>% p.adjust(p = .,method = "fdr")


#Write the fraction effect data.frame
write.table(x = merged,file  = "../cleandata/df_dds_res_amplicon_fungi_asvs_fraction_asvlevel.tsv",
            quote = F,sep = "\t",row.names = F,col.names = T)

