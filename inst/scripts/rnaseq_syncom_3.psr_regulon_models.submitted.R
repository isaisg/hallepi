library(ohchibi)
library(DESeq2)
#Set random seed
set.seed(130816)



Dat_raw <- readRDS(file = "../cleandata/dat_rnaseq_syncom.RDS")
Dat <- Dat_raw$Dat_rnaseq_syncom

### Create DEseq object ###
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,colData = Dat$Map,
                              design = ~ Experiment + group)

vsd <- vst(object = dds,blind = FALSE)
dds <- DESeq(dds)

#Extract contrasts of interests
res_axenic <- results(dds,contrast = c("group","NB_1000Pi","NB_50Pi")) %>%
  as.data.frame
res_syncom <- results(dds,contrast = c("group","Full_1000Pi","Full_50Pi")) %>%
  as.data.frame


#Low the psr genes
#Low the psr genes
regulon_psr <- read.table("../rawdata/regulon_psr.csv",header = F) %$%
  V1 %>% as.character



res_axenic <- match(regulon_psr,rownames(res_axenic)) %>%
  res_axenic[.,] %>% na.omit
res_syncom <- match(regulon_psr,rownames(res_syncom)) %>%
  res_syncom[.,] %>% na.omit

res_axenic %>% subset(padj < 0.1 & log2FoldChange < 0) %>% dim
res_syncom %>% subset(padj < 0.1 & log2FoldChange < 0 ) %>% dim

axenic <- res_axenic %>% subset(padj < 0.1 & log2FoldChange < 0) %>% rownames
full <- res_syncom %>% subset(padj < 0.1 & log2FoldChange < 0 ) %>% rownames

lista <- list(axenic = axenic,full = full)
saveRDS(object = lista,file = "../cleandata/list_psr_regulon_significance_insyncomrnaseq.RDS")

#Check the average fold change
2^(res_axenic$log2FoldChange %>% mean %>% abs)
2^(res_syncom$log2FoldChange %>% mean %>% abs)

