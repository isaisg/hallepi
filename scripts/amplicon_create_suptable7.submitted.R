library(ohchibi)
library(DESeq2)


set.seed(seed = 130816)

#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_otus_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$OTU

pval_thres <- 0.1

########## Determine ASVs that requiere PSR #############
###### Pi effect inside each genotype ######
res_inside_col0 <- results(object = dds,contrast = c("group","Col_0_low_Pi","Col_0_low")) %>% as.data.frame
res_inside_phf1 <- results(object = dds,contrast = c("group","phf1_low_Pi","phf1_low")) %>% as.data.frame
res_inside_phr1phl1 <- results(object = dds,contrast = c("group","phr1_phl1_low_Pi","phr1_phl1_low")) %>% as.data.frame

wt <- res_inside_col0 %>% subset(padj < pval_thres) %>% rownames
single <- res_inside_phf1 %>% subset(padj < pval_thres) %>% rownames
double <- res_inside_phr1phl1 %>% subset(padj < pval_thres) %>% rownames



res_inside_col0$OTU_Id <- rownames(res_inside_col0)
res_inside_col0$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_col0))
res_inside_col0$Genotype <- rep("Col-0",nrow(res_inside_col0))
res_inside_col0$Significance <- rep("NoSignificant",nrow(res_inside_col0))
res_inside_col0$Significance[which(res_inside_col0$OTU_Id %in% wt)] <- "Significant"
res_inside_col0$Direction <- rep("NoSignificant",nrow(res_inside_col0))
low <- res_inside_col0 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% OTU_Id %>% 
  as.character
up <- res_inside_col0 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% OTU_Id %>% 
  as.character
res_inside_col0$Direction[which(res_inside_col0$OTU_Id %in% low)] <- "low"
res_inside_col0$Direction[which(res_inside_col0$OTU_Id %in% up)] <- "low+Pi"
res_inside_col0$Significance <- res_inside_col0$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_col0$Direction <- res_inside_col0$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_col0 <- with(res_inside_col0,order(Significance,Direction)) %>%
  res_inside_col0[.,]
res_inside_col0 <- res_inside_col0[,c(7,1:6,8:11)]
rownames(res_inside_col0) <- NULL


res_inside_phf1$OTU_Id <- rownames(res_inside_phf1)
res_inside_phf1$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_phf1))
res_inside_phf1$Genotype <- rep("phf1",nrow(res_inside_phf1))
res_inside_phf1$Significance <- rep("NoSignificant",nrow(res_inside_phf1))
res_inside_phf1$Significance[which(res_inside_phf1$OTU_Id %in% wt)] <- "Significant"
res_inside_phf1$Direction <- rep("NoSignificant",nrow(res_inside_phf1))
low <- res_inside_phf1 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% OTU_Id %>% 
  as.character
up <- res_inside_phf1 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% OTU_Id %>% 
  as.character
res_inside_phf1$Direction[which(res_inside_phf1$OTU_Id %in% low)] <- "low"
res_inside_phf1$Direction[which(res_inside_phf1$OTU_Id %in% up)] <- "low+Pi"
res_inside_phf1$Significance <- res_inside_phf1$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_phf1$Direction <- res_inside_phf1$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_phf1 <- with(res_inside_phf1,order(Significance,Direction)) %>%
  res_inside_phf1[.,]
res_inside_phf1 <- res_inside_phf1[,c(7,1:6,8:11)]
rownames(res_inside_phf1) <- NULL


res_inside_phr1phl1$OTU_Id <- rownames(res_inside_phr1phl1)
res_inside_phr1phl1$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Genotype <- rep("phr1/phl1",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Significance <- rep("NoSignificant",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Significance[which(res_inside_phr1phl1$OTU_Id %in% wt)] <- "Significant"
res_inside_phr1phl1$Direction <- rep("NoSignificant",nrow(res_inside_phr1phl1))
low <- res_inside_phr1phl1 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% OTU_Id %>% 
  as.character
up <- res_inside_phr1phl1 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% OTU_Id %>% 
  as.character
res_inside_phr1phl1$Direction[which(res_inside_phr1phl1$OTU_Id %in% low)] <- "low"
res_inside_phr1phl1$Direction[which(res_inside_phr1phl1$OTU_Id %in% up)] <- "low+Pi"
res_inside_phr1phl1$Significance <- res_inside_phr1phl1$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_phr1phl1$Direction <- res_inside_phr1phl1$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_phr1phl1 <- with(res_inside_phr1phl1,order(Significance,Direction)) %>%
  res_inside_phr1phl1[.,]
res_inside_phr1phl1 <- res_inside_phr1phl1[,c(7,1:6,8:11)]
rownames(res_inside_phr1phl1) <- NULL

all_res <- rbind(res_inside_col0,res_inside_phf1,res_inside_phr1phl1)

#Read taxonomy
df_tax <- read.table(file = "../rawdata/df_taxonomy_fungi_otus.tsv",header = T,sep = "\t")
colnames(df_tax)[7] <- "OTU_Id"
suptab <- merge(all_res,df_tax, by = "OTU_Id")


#Clean the taxonomy columns
suptab$Kingdom <- suptab$Kingdom %>% gsub(pattern = "k__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")
suptab$Phylum <- suptab$Phylum %>% gsub(pattern = "p__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Class <- suptab$Class %>% gsub(pattern = "c__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")
suptab$Order <- suptab$Order %>% gsub(pattern = "o__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Family <- suptab$Family %>% gsub(pattern = "f__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Genus <- suptab$Genus %>% gsub(pattern = "g__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")

suptab$Num <- suptab$OTU_Id %>% gsub(pattern = "OTU_",replacement = "") %>% as.numeric
suptab$Genotype <- suptab$Genotype %>% factor(levels = c("Col-0","phf1","phr1/phl1"))
suptab <- with(suptab,order(Num,Genotype)) %>% suptab[.,]
suptab <- suptab[,-18]

#Write table
write.table(x = all_res,file = "../cleandata/sup_table_7.csv",append = F,quote = F,sep = ",",
            row.names = F,col.names = T)
