library(ohchibi)
library(RColorBrewer)
library(ggrepel)

set.seed(130816)
source(file = "plotting_parameters_hallepi.R")


Dat <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat <- Dat$RelativeAbundance

#Load enrichment profiles
res <- read.table(file = "../cleandata/df_dds_res_syncom_phosphateint.tsv",
           header = T,row.names = NULL,sep = "\t",quote = "",comment.char = "")
res$padj[which(is.na(res$padj))] <- 0

#Create a summary of both root and shoots
df_log <- aggregate(log2FoldChange~ Id,data = res,FUN = "mean")
df_padj <- aggregate(padj~ Id,data = res,FUN = "mean")
df_padj[which(df_padj$padj == 0),2] <- 1

merged <- merge(df_log,df_padj)
merged$padjtrans <- -log10(merged$padj)
merged$padjtrans[which(is.infinite(merged$padjtrans))] <- 0
merged$Significance <- rep("NoSignificant",nrow(merged))
merged$Significance[which(merged$padj < 0.1)] <- "Significant"
merged$Significance <- merged$Significance %>% factor
merged$Label <- rep(NA,nrow(merged))
chosen <- which(merged$Significance == "Significant")
merged$Label[chosen] <- merged$Id[chosen] %>% as.character
merged$Label <- merged$Label %>% gsub(pattern = "Sequence_",replacement = "USeq")


#Append the genus for each useq
df_useq2taxon <- read.table("../cleandata/df_useq2taxon.tsv",
                            header = T,sep = "\t",quote = "",comment.char = "")
colnames(df_useq2taxon)[2] <- "Label"
df_useq2taxon$Label <- df_useq2taxon$Label %>% as.character %>%
  gsub(pattern = "Sequence_",replacement = "USeq")

merged <- merge(merged,df_useq2taxon, by = "Label",all.x = TRUE)

p <- ggplot(data = merged,aes(log2FoldChange,padjtrans)) + 
  geom_vline(xintercept = 0 , size = 3,linetype = "dashed",color = "#D9D9D9") +
  geom_hline(yintercept = 0 , size = 3,linetype = "dashed",color = "#D9D9D9") +
  geom_point(aes(color = Significance,size = 8,alpha = 0.3)) + 
  geom_text_repel(aes(label = NLabels),size = 10,force = T) + 
  coord_flip()  +  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
        )  +
  ylab(label = "-log10 q.value")  +
  scale_color_manual(values = c("black","red"))
oh.save.pdf(p =p ,outname = "scatter_plot_interaction_syncom.pdf",outdir = "../figures/",
            height = 20,width = 15)


write.table(x = merged,file = "../data_figures/data_S7A.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


