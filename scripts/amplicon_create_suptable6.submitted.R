library(ohchibi)
library(DESeq2)


set.seed(seed = 130816)

#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_bacteria_asv_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$ASV

pval_thres <- 0.1

########## Determine ASVs that requiere PSR #############
###### Pi effect inside each genotype ######
res_inside_col0 <- results(object = dds,contrast = c("group","Col_0_low_Pi","Col_0_low")) %>% as.data.frame
res_inside_phf1 <- results(object = dds,contrast = c("group","phf1_low_Pi","phf1_low")) %>% as.data.frame 
res_inside_phr1phl1 <- results(object = dds,contrast = c("group","phr1_phl1_low_Pi","phr1_phl1_low")) %>% as.data.frame

wt <- res_inside_col0 %>% subset(padj < pval_thres) %>% rownames
single <- res_inside_phf1 %>% subset(padj < pval_thres) %>% rownames
double <- res_inside_phr1phl1 %>% subset(padj < pval_thres) %>% rownames

res_inside_col0$ASV_Id <- rownames(res_inside_col0)
res_inside_col0$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_col0))
res_inside_col0$Genotype <- rep("Col-0",nrow(res_inside_col0))
res_inside_col0$Significance <- rep("NoSignificant",nrow(res_inside_col0))
res_inside_col0$Significance[which(res_inside_col0$ASV_Id %in% wt)] <- "Significant"
res_inside_col0$Direction <- rep("NoSignificant",nrow(res_inside_col0))
low <- res_inside_col0 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% ASV_Id %>% 
  as.character
up <- res_inside_col0 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% ASV_Id %>% 
  as.character
res_inside_col0$Direction[which(res_inside_col0$ASV_Id %in% low)] <- "low"
res_inside_col0$Direction[which(res_inside_col0$ASV_Id %in% up)] <- "low+Pi"
res_inside_col0$Significance <- res_inside_col0$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_col0$Direction <- res_inside_col0$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_col0 <- with(res_inside_col0,order(Significance,Direction)) %>%
  res_inside_col0[.,]
res_inside_col0 <- res_inside_col0[,c(7,1:6,8:11)]
rownames(res_inside_col0) <- NULL


res_inside_phf1$ASV_Id <- rownames(res_inside_phf1)
res_inside_phf1$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_phf1))
res_inside_phf1$Genotype <- rep("phf1",nrow(res_inside_phf1))
res_inside_phf1$Significance <- rep("NoSignificant",nrow(res_inside_phf1))
res_inside_phf1$Significance[which(res_inside_phf1$ASV_Id %in% wt)] <- "Significant"
res_inside_phf1$Direction <- rep("NoSignificant",nrow(res_inside_phf1))
low <- res_inside_phf1 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% ASV_Id %>% 
  as.character
up <- res_inside_phf1 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% ASV_Id %>% 
  as.character
res_inside_phf1$Direction[which(res_inside_phf1$ASV_Id %in% low)] <- "low"
res_inside_phf1$Direction[which(res_inside_phf1$ASV_Id %in% up)] <- "low+Pi"
res_inside_phf1$Significance <- res_inside_phf1$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_phf1$Direction <- res_inside_phf1$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_phf1 <- with(res_inside_phf1,order(Significance,Direction)) %>%
  res_inside_phf1[.,]
res_inside_phf1 <- res_inside_phf1[,c(7,1:6,8:11)]
rownames(res_inside_phf1) <- NULL


res_inside_phr1phl1$ASV_Id <- rownames(res_inside_phr1phl1)
res_inside_phr1phl1$Contrast <- rep("low_Pi_vs_low",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Genotype <- rep("phr1/phl1",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Significance <- rep("NoSignificant",nrow(res_inside_phr1phl1))
res_inside_phr1phl1$Significance[which(res_inside_phr1phl1$ASV_Id %in% wt)] <- "Significant"
res_inside_phr1phl1$Direction <- rep("NoSignificant",nrow(res_inside_phr1phl1))
low <- res_inside_phr1phl1 %>% subset(Significance == "Significant" & log2FoldChange <0 ) %$% ASV_Id %>% 
  as.character
up <- res_inside_phr1phl1 %>% subset(Significance == "Significant" & log2FoldChange >0 ) %$% ASV_Id %>% 
  as.character
res_inside_phr1phl1$Direction[which(res_inside_phr1phl1$ASV_Id %in% low)] <- "low"
res_inside_phr1phl1$Direction[which(res_inside_phr1phl1$ASV_Id %in% up)] <- "low+Pi"
res_inside_phr1phl1$Significance <- res_inside_phr1phl1$Significance %>% factor(levels = c("Significant","NoSignificant"))
res_inside_phr1phl1$Direction <- res_inside_phr1phl1$Direction %>% factor(levels = c("low","low+Pi","NoSignificant"))
res_inside_phr1phl1 <- with(res_inside_phr1phl1,order(Significance,Direction)) %>%
  res_inside_phr1phl1[.,]
res_inside_phr1phl1 <- res_inside_phr1phl1[,c(7,1:6,8:11)]
rownames(res_inside_phr1phl1) <- NULL

all_res <- rbind(res_inside_col0,res_inside_phf1,res_inside_phr1phl1)

#Read taxonomy
df_tax <- read.table(file = "../rawdata/dada2_structure_maxee0_nopool.taxaclassification.tsv",header = T,sep = "\t")
df_tax <- df_tax[,-2]
all_res <- merge(all_res,df_tax, by = "ASV_Id")
all_res$Num <- all_res$ASV_Id %>% gsub(pattern = "ASV",replacement = "") %>% as.numeric
all_res$Genotype <- all_res$Genotype %>% factor(levels = c("Col-0","phf1","phr1/phl1"))
all_res <- with(all_res,order(Num,Genotype)) %>% all_res[.,]
all_res <- all_res[,-18]

#Write table
write.table(x = all_res,file = "../cleandata/sup_table_6.csv",append = F,quote = F,sep = ",",
            row.names = F,col.names = T)
