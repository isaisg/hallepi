library(ohchibi)
library(DESeq2)
library(paletteer)
library(palettesPM)
library(scales)
library(UpSetR)

source('plotting_parameters_hallepi.R')
set.seed(seed = 130816)

#Done
#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_bacteria_asv_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$ASV

pval_thres <- 0.1

########## Determine ASVs that requiere PSR #############
###### Pi effect inside each genotype ######
res_inside_col0 <- results(object = dds,contrast = c("group","Col_0_low_Pi","Col_0_low"))
res_inside_phf1 <- results(object = dds,contrast = c("group","phf1_low_Pi","phf1_low"))
res_inside_phr1phl1 <- results(object = dds,contrast = c("group","phr1_phl1_low_Pi","phr1_phl1_low"))

#Rbidn results for dataset object
res_inside_col0$Genotype <- "Col-0"
res_inside_phf1$Genotype <- "phf1"
res_inside_phr1phl1$Genotype <- "phr1/phl1"
towrite <- rbind(res_inside_col0,res_inside_phf1,res_inside_phr1phl1) %>% as.data.frame

write.table(x = towrite,file = "../data_figures/data_Fig3C_enrichedinlow.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)





wt <- res_inside_col0 %>% subset(padj < pval_thres & log2FoldChange < 0) %>% rownames
single <- res_inside_phf1 %>% subset(padj < pval_thres & log2FoldChange < 0) %>% rownames
double <- res_inside_phr1phl1 %>% subset(padj < pval_thres & log2FoldChange < 0) %>% rownames


lista <- list(wt = wt,single = single,double = double)
pdf(file = "../figures/upsetr_bacteria_required.pdf")
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

lista <- list(wt_specific = wt_specific, single_specific = single_specific, double_specific,
      int_all = int_all, int_wt_single = int_wt_single, int_wt_double = int_wt_double,
      int_single_double = int_single_double)

saveRDS(object = lista,file = "../cleandata/list_bacteria_asv_required_genotypephosphate_root.RDS")

################# Create dds object to plot ####################
#Create a df with zscore values
melted_z_all <- varianceStabilizingTransformation(dds) %>%assay %>%
  as.matrix %>% t %>% scale %>% t %>% melt
colnames(melted_z_all) <- c("Id","DADA2_Header","value")
Map <- dds %>% colData() %>% as.data.frame
melted_z_all <- merge(melted_z_all,Map, by = "DADA2_Header") 
melted_sum <- melted_z_all %>% 
  dcast(data = .,Id ~ Genotype+Phosphate,value.var = "value",fun.aggregate = mean) %>%
  melt
melted_sum$Genotype <- rep("phr1/phl1",nrow(melted_sum))
melted_sum$Genotype[melted_sum$variable %>% grep(pattern = "phf1")] <- "phf1"
melted_sum$Genotype[melted_sum$variable %>% grep(pattern = "Col")] <- "Col-0"
melted_sum$Genotype <- melted_sum$Genotype %>% factor(levels = c("Col-0","phf1","phr1/phl1"))
melted_sum$Phosphate <- rep("low",nrow(melted_sum))
melted_sum$Phosphate[melted_sum$variable %>% grep(pattern = "low_Pi")] <- "low+Pi"
melted_sum$Phosphate <- melted_sum$Phosphate %>% factor


#Cluster them 
asvs_required <- wt_specific
clust_required <- which(melted_sum$Id %in% asvs_required) %>% melted_sum[.,] %>%
  droplevels %>% acast(formula = Id~variable,value.var = "value",fun.aggregate = sum) %>% 
  dist %>% hclust(method = "ward.D")
ord_required <- clust_required$order %>% clust_required$labels[.]
melted_sum_required <- which(melted_sum$Id %in% asvs_required) %>% melted_sum[.,] %>%
  droplevels

### Add the abundance patterns 
#Retrieve abundance of each asv
Dat_asv <- readRDS(file = "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
melted_ab_all <- Dat_asv$RelativeAbundance$Tab %>% melt
colnames(melted_ab_all) <- c("Id","DADA2_Header","RA")
map_ab <- Dat_asv$RelativeAbundance$Map 
melted_ab_all <- merge(melted_ab_all,map_ab, by = "DADA2_Header")
#Compute average RE in fraction 
melted_ab <- dcast(data = melted_ab_all,formula = Id ~ Fraction,fun.aggregate = mean,value.var = "RA")

melted_sum_required <- merge(melted_sum_required,melted_ab, by = "Id")

#Taxonomy
df_tax <- read.table("../rawdata/dada2_structure_maxee0_nopool.taxaclassification.tsv",
                     header = T,sep = "\t",quote = "")
colnames(df_tax)[1] <- "Id"
melted_sum_required <- merge(melted_sum_required,df_tax , by = "Id")


### Create taxonomic classification 
mphyla <- palettesPM::pm.names.phyla()
melted_sum_required$mPhylum <- melted_sum_required$Phylum %>% as.character
melted_sum_required$mPhylum[which(!(melted_sum_required$Phylum %in% mphyla))] <- "Other"

#Is Burkholderia
#Create a burkholderiace column
melted_sum_required$Burk <- rep("NO",nrow(melted_sum_required))
melted_sum_required$Burk[which(melted_sum_required$Family == "Burkholderiaceae")] <- "YES"
melted_sum_required$Burk <- melted_sum_required$Burk %>% factor
melted_sum_required$IsBurk <- rep("IsBurk",nrow(melted_sum_required)) %>% factor


#Order according to the individual clusterings
melted_sum_required$Id <- melted_sum_required$Id %>% factor(levels = ord_required)

###### Figures for required 
p_heatmap <- ggplot(data = melted_sum_required,aes(Phosphate,Id, fill = value)) +
  geom_tile() +
  facet_grid(~Genotype,space = "free",scales = "free") +
  scale_fill_gradient2(midpoint = 0,mid = "white",high = "#E65037",
                       low = "#2480ff",guide = "colourbar",limits = c(-1,1),oob = squish) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 15,
                size_strip_text.x = 15,size_panel_border = 0.5)  +
  theme(axis.ticks = element_blank()) +
  coord_cartesian(expand = FALSE) + theme(legend.position = "none")

#Create heatmap to extract legend
p_legend <- ggplot(data = melted_sum_required,aes(Phosphate,Id, fill = value)) +
  geom_tile() +
  facet_grid(~Genotype,space = "free",scales = "free") +
  scale_fill_gradient2(midpoint = 0,mid = "white",high = "#E65037",
                       low = "#2480ff",guide = "colourbar",limits = c(-1,1),oob = squish) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 15,
                size_strip_text.x = 15,size_panel_border = 0.5)  +
  theme(axis.ticks = element_blank()) +
  coord_cartesian(expand = FALSE)
oh.save.pdf(p = p_legend,outname = "figure3_heatmap_legend.pdf",outdir = "../figures/")


melted_sum_required_forplot <- melted_sum_required[,c(1,6:9,11:19)] %>% unique
p_abundance <- ggplot(data = melted_sum_required_forplot,aes(Id,(Root))) + geom_bar(stat = "identity") + 
  coord_flip(expand = FALSE,ylim = c(0,0.001)) + 
  theme_ohchibi(size_axis_text.x = 12,size_axis_text.y =0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle = 270),
        axis.ticks.y = element_blank(),panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())   + theme(legend.position = "none")  

p_phyla <-  ggplot(data = melted_sum_required_forplot,mapping = aes(Id,Kingdom, fill = mPhylum)) + 
  geom_tile() +
  coord_flip(expand = FALSE) + 
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(legend.position = "none",axis.ticks = element_blank())  + 
  scale_fill_phyla() 


p_burk <- ggplot(data = melted_sum_required_forplot,mapping = aes(Id,IsBurk, fill = Burk)) +  geom_tile() +
  coord_flip(expand = FALSE) + 
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(legend.position = "none",axis.ticks = element_blank())  + 
  scale_fill_manual(values = c("white","black")) 

composition <- egg::ggarrange(p_phyla,p_burk,p_heatmap,p_abundance,nrow = 1,widths = c(0.1,0.1,1,0.4))

oh.save.pdf(p = composition,outname = "figure3_bacteria_asv_genotypebyphosphate_asvs_requiredpsr.pdf",
            outdir = "../figures/",width = 10,height = 20,pointsize = 12)

#Check the total abundance in root of all those ASVs
melted_ab_sub <- lista %>% unlist %>% unique %>% match(melted_ab$Id) %>%
  melted_ab[.,] %>% droplevels
melted_ab_sub$Root %>% sum

########## Determine ASVs that arre inhibyted by  PSR #############
rm(list=ls())
source('plotting_parameters_hallepi.R')
set.seed(seed = 130816)

#Load dataset
Dat_dds <- readRDS(file = "../cleandata/dds_deseq2_hallepi_amplicon_bacteria_asv_groupsgenotypephosphate.RDS")
dds <- Dat_dds$Root$ASV

pval_thres <- 0.1


###### Pi effect inside each genotype ######
res_inside_col0 <- results(object = dds,contrast = c("group","Col_0_low_Pi","Col_0_low"))
res_inside_phf1 <- results(object = dds,contrast = c("group","phf1_low_Pi","phf1_low"))
res_inside_phr1phl1 <- results(object = dds,contrast = c("group","phr1_phl1_low_Pi","phr1_phl1_low"))

#Rbidn results for dataset object
res_inside_col0$Genotype <- "Col-0"
res_inside_phf1$Genotype <- "phf1"
res_inside_phr1phl1$Genotype <- "phr1/phl1"
towrite <- rbind(res_inside_col0,res_inside_phf1,res_inside_phr1phl1) %>% as.data.frame

write.table(x = towrite,file = "../data_figures/data_Fig3C_enrichedinlowplusP.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)




wt <- res_inside_col0 %>% subset(padj < pval_thres & log2FoldChange > 0) %>% rownames
single <- res_inside_phf1 %>% subset(padj < pval_thres & log2FoldChange > 0) %>% rownames
double <- res_inside_phr1phl1 %>% subset(padj < pval_thres & log2FoldChange > 0) %>% rownames


lista <- list(wt = wt,single = single,double = double)
pdf(file = "../figures/upsetr_bacteria_inhibited.pdf")
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

lista <- list(wt_specific = wt_specific, single_specific = single_specific, double_specific,
              int_all = int_all, int_wt_single = int_wt_single, int_wt_double = int_wt_double,
              int_single_double = int_single_double)
saveRDS(object = lista,file = "../cleandata/list_bacteria_asv_inhibited_genotypephosphate_root.RDS")

################# Create dds object to plot ####################
#Create a df with zscore values
melted_z_all <- varianceStabilizingTransformation(dds) %>%assay %>%
  as.matrix %>% t %>% scale %>% t %>% melt
colnames(melted_z_all) <- c("Id","DADA2_Header","value")
Map <- dds %>% colData() %>% as.data.frame
melted_z_all <- merge(melted_z_all,Map, by = "DADA2_Header") 
melted_sum <- melted_z_all %>% 
  dcast(data = .,Id ~ Genotype+Phosphate,value.var = "value",fun.aggregate = mean) %>%
  melt
melted_sum$Genotype <- rep("phr1/phl1",nrow(melted_sum))
melted_sum$Genotype[melted_sum$variable %>% grep(pattern = "phf1")] <- "phf1"
melted_sum$Genotype[melted_sum$variable %>% grep(pattern = "Col")] <- "Col-0"
melted_sum$Genotype <- melted_sum$Genotype %>% factor(levels = c("Col-0","phf1","phr1/phl1"))
melted_sum$Phosphate <- rep("low",nrow(melted_sum))
melted_sum$Phosphate[melted_sum$variable %>% grep(pattern = "low_Pi")] <- "low+Pi"
melted_sum$Phosphate <- melted_sum$Phosphate %>% factor


#Cluster them 
asvs_required <- wt_specific
clust_required <- which(melted_sum$Id %in% asvs_required) %>% melted_sum[.,] %>%
  droplevels %>% acast(formula = Id~variable,value.var = "value",fun.aggregate = sum) %>% 
  dist %>% hclust(method = "ward.D")
ord_required <- clust_required$order %>% clust_required$labels[.]
melted_sum_required <- which(melted_sum$Id %in% asvs_required) %>% melted_sum[.,] %>%
  droplevels

### Add the abundance patterns 
#Retrieve abundance of each asv
Dat_asv <- readRDS(file = "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
melted_ab_all <- Dat_asv$RelativeAbundance$Tab %>% melt
colnames(melted_ab_all) <- c("Id","DADA2_Header","RA")
map_ab <- Dat_asv$RelativeAbundance$Map 
melted_ab_all <- merge(melted_ab_all,map_ab, by = "DADA2_Header")
#Compute average RE in fraction 
melted_ab <- dcast(data = melted_ab_all,formula = Id ~ Fraction,fun.aggregate = mean,value.var = "RA")

melted_sum_required <- merge(melted_sum_required,melted_ab, by = "Id")

#Taxonomy
df_tax <- read.table("../rawdata/dada2_structure_maxee0_nopool.taxaclassification.tsv",
                     header = T,sep = "\t",quote = "")
colnames(df_tax)[1] <- "Id"
melted_sum_required <- merge(melted_sum_required,df_tax , by = "Id")


### Create taxonomic classification 
mphyla <- palettesPM::pm.names.phyla()
melted_sum_required$mPhylum <- melted_sum_required$Phylum %>% as.character
melted_sum_required$mPhylum[which(!(melted_sum_required$Phylum %in% mphyla))] <- "Other"

#Is Burkholderia
#Create a burkholderiace column
melted_sum_required$Burk <- rep("NO",nrow(melted_sum_required))
melted_sum_required$Burk[which(melted_sum_required$Family == "Burkholderiaceae")] <- "YES"
melted_sum_required$Burk <- melted_sum_required$Burk %>% factor
melted_sum_required$IsBurk <- rep("IsBurk",nrow(melted_sum_required)) %>% factor


#Order according to the individual clusterings
melted_sum_required$Id <- melted_sum_required$Id %>% factor(levels = ord_required)

###### Figures for required 
p_heatmap <- ggplot(data = melted_sum_required,aes(Phosphate,Id, fill = value)) +
  geom_tile() +
  facet_grid(~Genotype,space = "free",scales = "free") +
  scale_fill_gradient2(midpoint = 0,mid = "white",high = "#E65037",
                       low = "#2480ff",guide = "colourbar",limits = c(-1,1),oob = squish) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 15,
                size_strip_text.x = 15,size_panel_border = 0.5)  +
  theme(axis.ticks = element_blank()) +
  coord_cartesian(expand = FALSE) + theme(legend.position = "none")


melted_sum_required_forplot <- melted_sum_required[,c(1,6:9,11:19)] %>% unique
p_abundance <- ggplot(data = melted_sum_required_forplot,aes(Id,(Root))) + geom_bar(stat = "identity") + 
  coord_flip(expand = FALSE,ylim = c(0,0.001)) + 
  theme_ohchibi(size_axis_text.x = 12,size_axis_text.y =0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.minor.x = element_blank(),axis.text.x = element_text(angle = 270),
        axis.ticks.y = element_blank(),panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())   + theme(legend.position = "none")  

p_phyla <-  ggplot(data = melted_sum_required_forplot,mapping = aes(Id,Kingdom, fill = mPhylum)) + 
  geom_tile() +
  coord_flip(expand = FALSE) + 
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(legend.position = "none",axis.ticks = element_blank())  + 
  scale_fill_phyla() 


p_burk <- ggplot(data = melted_sum_required_forplot,mapping = aes(Id,IsBurk, fill = Burk)) +  geom_tile() +
  coord_flip(expand = FALSE) + 
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.y = 0,size_strip_text.x = 0,
                size_lines_panel = 0,size_panel_border = 0.5)  +
  theme(legend.position = "none",axis.ticks = element_blank())  + 
  scale_fill_manual(values = c("white","black")) 

composition <- egg::ggarrange(p_phyla,p_burk,p_heatmap,p_abundance,nrow = 1,widths = c(0.1,0.1,1,0.4))
oh.save.pdf(p = composition,outname = "figure3_bacteria_asv_genotypebyphosphate_asvs_inhibitedpsr.pdf",
           outdir = "../figures/",width = 10,height = 20,pointsize = 12)

#Check the total abundance in root of all those ASVs
melted_ab_sub <- lista %>% unlist %>% unique %>% match(melted_ab$Id) %>%
  melted_ab[.,] %>% droplevels
melted_ab_sub$Root %>% sum

