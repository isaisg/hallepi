library(ohchibi)
library(DESeq2)



#Read taxonomy
df_tax <- read.table(file = "../rawdata/df_taxonomy_fungi_otus.tsv",
                     header = T,sep = "\t",quote = "",comment.char = "")
colnames(df_tax)[7] <- "OTU_Id"

dds_res <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_otus_maineffects.RDS")

########################## OTUs ###########################################
### Root vs Soil ###
dds <- dds_res$OTUs

res_rs <- dds %>% results(contrast = c("Fraction","Root","Soil")) %>%
  as.data.frame
res_rs$Contrast <- rep("Root-Soil",nrow(res_rs))
res_rs$OTU_Id <- rownames(res_rs)

### Shoot vs Soil ###
res_ss <- dds %>% results(contrast = c("Fraction","Shoot","Soil")) %>%
  as.data.frame
res_ss$Contrast <- rep("Shoot-Soil",nrow(res_ss))
res_ss$OTU_Id <- rownames(res_ss)

### Shoot vs Root ###
res_sr <- dds %>% results(contrast = c("Fraction","Root","Shoot")) %>%
  as.data.frame

res_sr$Contrast <- rep("Root-Shoot",nrow(res_sr))
res_sr$OTU_Id <- rownames(res_sr)


res <- rbind(res_rs,res_ss,res_sr)
rownames(res) <- NULL
merged <- merge(res,df_tax, by = "OTU_Id")

#Correc the pvalues globally
merged$padj <- merged$pvalue %>% p.adjust(p = .,method = "fdr")


#Write the fraction effect data.frame
write.table(x = merged,file  = "../cleandata/df_dds_res_amplicon_fungi_otus_fraction_asvlevel.tsv",
            quote = F,sep = "\t",row.names = F,col.names = T)

