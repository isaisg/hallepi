library(ohchibi)
library(DESeq2)
library(paletteer)
library(clusterProfiler)
library(org.At.tair.db)
library(egg)
library(pheatmap)

#Set random seed
set.seed(130816)


source('plotting_parameters_hallepi.R')


Dat_rnaseq <- readRDS(file = "../cleandata/dat_rnaseq_hallepi.RDS")
Dat <- Dat_rnaseq$Dat_rnaseq_halle

Dat <- Dat %>% subset.Dataset(subset = Phosphate == "low" |
                                Phosphate == "low_Pi",drop = T,clean = T)


# Create deseq object
Dat_sub <- Dat %>% subset.Dataset(Fraction == "Root",drop = T,clean = T)

dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design = ~group)

vsd <- vst(object = dds,blind = FALSE)
dds <- DESeq(dds,test = "LRT",reduced = ~1)

res <- results(dds) %>% as.data.frame
sig_genes <- res %>% subset(padj < 0.1) %>% rownames

Tab_vsd <- match(sig_genes,rownames(assay(vsd))) %>%
  assay(vsd)[.,]

### Scale the data
Tab_vsd_Z <- Tab_vsd %>% t %>% scale %>% t %>% as.matrix


#Deal with outliers
maxv <- quantile(Tab_vsd_Z,0.99) %>% as.numeric
minv <- quantile(Tab_vsd_Z,0.01) %>% as.numeric

Tab_vsd_Z[Tab_vsd_Z>maxv] <- maxv
Tab_vsd_Z[Tab_vsd_Z<minv] <- minv

#Summarize the columns
df_melted <- Tab_vsd_Z %>% melt
colnames(df_melted) <- c("genes","sample","value")
df_melted <- merge(df_melted,Dat_sub$Map,by = "sample")
Tab_vsd_Z <- acast(data = df_melted,value.var = "value",fun.aggregate = mean,
                   formula = genes~group)
# 
# 
pheatmap(Tab_vsd_Z,clustering_distance_rows = "euclidean",
          clustering_method = "ward.D2",cluster_cols = F,cutree_rows = 5)
#Cluster the following object#
genes_clust <- dist(x = Tab_vsd_Z,method = "euclidean") %>%
  hclust(method = "ward.D2") 
genes_order <- genes_clust$order %>% genes_clust$labels[.]
df_clusts <- genes_clust %>% cutree(tree = .,k = 5) %>%
  data.frame(Gene = names(.),Cluster = paste0("Cluster",.),row.names = NULL)
df_clusts <- df_clusts[,-1]

#Melt matrix
melted_chosen <- Tab_vsd_Z %>% melt
colnames(melted_chosen)[1] <- "Gene"
colnames(melted_chosen)[2] <- "group"


merged <- merge(melted_chosen,df_clusts,by = "Gene")
merged$Gene <- merged$Gene %>% factor(levels = genes_order)

merged$Genotype <- rep("phr1_phl1",nrow(merged))
merged$Genotype[merged$group %>% grep(pattern = "Col_0")] <- "Col_0"
merged$Genotype[merged$group %>% grep(pattern = "phf1")] <- "phf1"
merged$Genotype <- merged$Genotype %>% factor(levels = c("Col_0","phf1","phr1_phl1"))

merged$Phosphate <- rep("low",nrow(merged))
merged$Phosphate[merged$group %>% grep(pattern = "low_Pi")] <- "low_Pi"
merged$Phosphate <- merged$Phosphate %>% factor

order_clusters <- merged$Gene %>% levels %>% match(.,merged$Gene) %>%
  merged$Cluster[.] %>% unique %>% as.character
merged$Cluster <- merged$Cluster %>% factor(levels = order_clusters)


c2 <- merged %>% subset(Cluster == "Cluster2") %>% droplevels
c2$Cluster <- c2$Cluster %>% gsub(pattern = "2",replacement = "1")


c5 <- merged %>% subset(Cluster == "Cluster5") %>% droplevels
c5$Cluster <- c5$Cluster %>% gsub(pattern = "5",replacement = "2")

c4 <- merged %>% subset(Cluster == "Cluster4") %>% droplevels
c4$Cluster <- c4$Cluster %>% gsub(pattern = "4",replacement = "3")


c1 <- merged %>% subset(Cluster == "Cluster1") %>% droplevels
c1$Cluster <- c1$Cluster %>% gsub(pattern = "1",replacement = "4")

c3<- merged %>% subset(Cluster == "Cluster3") %>% droplevels
c3$Cluster <- c3$Cluster %>% gsub(pattern = "3",replacement = "5")

merged_obs <- rbind(c2,c5,c1,c3,c4)
merged_obs$Cluster <- merged_obs$Cluster %>% factor



p <- ggplot(data = merged_obs,aes(Phosphate,Gene, fill = value)) +
  geom_tile() + facet_grid(Cluster ~ Genotype,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",
                         palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2),na.value = "white")

oh.save.pdf(p = p,outname = "figure1_heatmap_phosphate_legend.pdf",
            height = 20,width = 10,outdir = "../figures/")  

p <- ggplot(data = merged_obs,aes(Phosphate,Gene, fill = value)) +
  geom_tile() + facet_grid(Cluster ~ Genotype,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",
                         palette = "kovesi.diverging_bwr_55_98_c37",
                         limits = c(-2,2),na.value = "white") +
  theme_ohchibi(size_axis_title.y = 0,size_axis_text.y = 0,
                size_strip_text.y = 30,size_axis_title.x = 0,font_family = "Arial") +
  theme(legend.position = "none",axis.ticks.y = element_blank()) +
  coord_cartesian(expand = FALSE)

#Read the psr regulon
psr_regulon <- read.table(file = "../rawdata/regulon_psr.csv") %$% V1 %>%
  as.character

df_reg <- data.frame(Gene = levels(merged_obs$Gene), 
                     Regulon = rep("No",length(levels(merged_obs$Gene))),stringsAsFactors = F)
df_reg$Regulon[which(df_reg$Gene %in% psr_regulon)] <- "Yes"
df_reg$Gene <- df_reg$Gene %>% factor(levels = levels(merged_obs$Gene))
df_reg$Type <- rep("PSR",nrow(df_reg)) %>% factor

df_reg <- merged_obs[,c(1,4)] %>% merge(df_reg,.,by = "Gene")
df_reg <- df_reg %>% unique
p_psr <- ggplot(data = df_reg,aes(Type,Gene, fill = Regulon)) +
  geom_raster() +  facet_grid(Cluster ~ .,space = "free",scales = "free") +
  scale_fill_manual(values = c("white","black")) + 
  theme_ohchibi(font_family = "Arial") +
  coord_cartesian(expand = FALSE)  +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),
        axis.title.x  = element_blank(),axis.text.y = element_blank(),
        axis.title.y = element_blank(),strip.text.y = element_blank())


oh.save.pdf(p = p_psr,outname = "figure1_bar_psr.pdf",
            height = 20,width = 10,outdir = "../figures/")  

p_psr <- p_psr + theme(legend.position = "none")

composition <- egg::ggarrange(p,p_psr,nrow = 1,widths = c(1,0.05))
oh.save.pdf(composition,outname = "figure1_heatmap_phosphate.pdf",
            height = 20,width = 10,outdir = "../figures/") 
#Get the identity of the genes
df_reg %>% subset(Regulon == "Yes") %>% droplevels %>%
  write.table(x = .,file = "../cleandata/df_7genes_rnaseq.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)


#Compute gene enrichment
lista <- list(
  Cluster1 = merged_obs %>% subset(Cluster == "Cluster1") %$% Gene %>% unique,
  Cluster2 = merged_obs %>% subset(Cluster == "Cluster2") %$% Gene %>% unique,
  Cluster3 = merged_obs %>% subset(Cluster == "Cluster3") %$% Gene %>% unique,
  Cluster4 = merged_obs %>% subset(Cluster == "Cluster4") %$% Gene %>% unique,
  Cluster5 = merged_obs %>% subset(Cluster == "Cluster5") %$% Gene %>% unique
)


cg <- compareCluster(geneCluster=lista,
                     fun="enrichGO",
                     keyType       = "TAIR",
                     OrgDb         = org.At.tair.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)


p <- dotplot(cg,showCategory = 20,font.size = 15) +
  geom_vline(xintercept = c(1.5,2.5),size = size_vline,color = "#D9D9D9") +
  theme(panel.grid.major.x = element_blank()) +
  scale_color_paletteer_c(package = "viridis",palette = "viridis")
oh.save.pdf(p = p,outname = "figure1_heatmap_phosphate_enrichment.pdf",
            height = 10,width = 10,outdir = "../figures/")  

saveRDS(object = lista,file = "../cleandata/list_clusters_genes_main.RDS")



### Create supplementary Table 2
res$Gene <- rownames(res)
res$Significance <- rep("NoSignificant",nrow(res))
res$Significance[which(res$Gene %in% sig_genes)] <- "Significant"

res$Cluster <- rep("NoSignificant",nrow(res))
res$Cluster[which(res$Gene %in% lista$Cluster1)] <- "Cluster1"
res$Cluster[which(res$Gene %in% lista$Cluster2)] <- "Cluster2"
res$Cluster[which(res$Gene %in% lista$Cluster3)] <- "Cluster3"
res$Cluster[which(res$Gene %in% lista$Cluster4)] <- "Cluster4"
res$Cluster[which(res$Gene %in% lista$Cluster5)] <- "Cluster5"

#Add the PSR belonging 
res$PSR_Regulon <- rep("No",nrow(res))
res$PSR_Regulon[which(res$Gene %in% psr_regulon)] <- "Yes"

res <- res[,c(7,1:6,8,9,10)]
rownames(res) <- NULL

#Order the table
res <- with(res,order(Cluster,padj)) %>%
  res[.,]
res <- res[,c(1,6:10)]
write.table(x = res,file = "../cleandata/sup_table_2.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)
