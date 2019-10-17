library(ohchibi)
library(palettesPM)
library(paletteer)
library(DESeq2)
library(emmeans)

set.seed(130816)

size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4

setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts')

Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_piDO_useq97.RDS")



#Plot abundance of burkholderia
Dat_rab <- Dat_ori$RelativeAbundance
Dat_rab <- Dat_rab %>% subset.Dataset(Fraction == "Root" & SynCom == "Full" ,drop = T,clean = T)
melted <- Dat_rab$Tab %>% melt
colnames(melted) <- c("Id","ID_Matrix","RA")
melted <- merge(melted,Dat_rab$Map, by = "ID_Matrix")

melted_sub <- which(melted$Id %in% c("Sequence_16","Sequence_30")) %>%
  melted[.,] %>% droplevels

df_sub <- aggregate(RA~ID_Matrix,data = melted_sub,sum)
df_sub <- merge(df_sub,Dat_rab$Map, by = "ID_Matrix", all.x = TRUE)



p <- df_sub %>% subset(RA < 0.1) %>% 
  chibi.boxplot(Map = .,x_val = "Pi",y_val = "RA",
              style = "open",
              facet_formula = "Genotype",size_median = 2,
              size_point = 0,stroke_point = 0,size_axis_text.y = 35,
              size_axis_text.x = 30,median_color = "red",
              size_axis_title.y = 40,strip_text_size = 30,
              size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "Phosphate") + 
  xlab(label = "Syncom") +
  geom_sina(alpha = 0.15, color = "#414141",
            stroke = 0,size = 4) +
  #geom_text(data = Res_em,
  #          aes(x = Syncom,y = 0.03,label = Letters),
  #          inherit.aes = F,size = 7.5,family ="Arial") +
  coord_cartesian(expand = TRUE) +
  theme(
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,family = "Arial",face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 25,family = "Arial",face = "bold")
        )   + 
  ylab(label = "Relative abundance in root")
oh.save.pdf(p = p,outname = "pi_burkholderia_mutants.pdf",outdir = "../figures/")


## Generate model to confirm trends observed
Dat_raw <- Dat_ori$RawCounts

#Check inside root
Dat_sub <- Dat_raw %>% subset.Dataset(Fraction == "Root" & SynCom == "Full" ,drop = T,clean = T)
Dat_sub <- Dat_sub %>% collapse_by_taxonomy.Dataset(level = 7)

#Create Deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design =~ Genotype + Pi + Genotype:Pi)
dds <- DESeq(object = dds,test = "LRT",reduced = ~Genotype+Pi)

df_temp <- dds %>% results %>% as.data.frame 
df_int <- df_temp[rownames(df_temp) %in% 
                    c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design =~ Genotype + Pi )
dds <- DESeq(object = dds,test = "LRT",reduced = ~Genotype)

df_temp <- dds %>% results %>% as.data.frame 
df_pi <- df_temp[rownames(df_temp) %in% 
                    c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design =~ Genotype + Pi )
dds <- DESeq(object = dds,test = "LRT",reduced = ~Pi)

df_temp <- dds %>% results %>% as.data.frame 
df_geno <- df_temp[rownames(df_temp) %in% 
                   c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df <- rbind(df_geno,df_pi,df_int)
df$Type <- c("Genotype","Pi","Interaction")
rownames(df) <- NULL
#No interaction so check main effect
dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design =~ Genotype + Pi )
dds <- DESeq(object = dds)

df_temp <- dds %>% results(contrast = c("Genotype","phf1","Col-0")) %>%
  as.data.frame

df_temp[rownames(df_temp) %in% 
          c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("Genotype","phr1/phl1","Col-0")) %>%
  as.data.frame

df_temp[rownames(df_temp) %in% 
          c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]



df_temp <- dds %>% results(contrast = c("Pi","50uM","0uM")) %>%
  as.data.frame 
df_temp[rownames(df_temp) %in% 
          c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("Pi","1000uM","0uM")) %>%
  as.data.frame 
df_temp[rownames(df_temp) %in% 
          c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

df_temp <- dds %>% results(contrast = c("Pi","50uM","0uM")) %>%
  as.data.frame 
df_temp[rownames(df_temp) %in% 
          c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]







#Do a group model to compare pairs
Dat_sub$Map$group <- paste0(Dat_sub$Map$Genotype,"_",Dat_sub$Map$Pi)
dds <- DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                              design =~ group )
dds <- DESeq(object = dds)

#Do pairwise comparisons
#Col -0
df_temp <- dds %>% results(contrast = c("group","Col-0_50uM","phf1_50uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("group","Col-0_50uM","Col-0_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

df_temp <- dds %>% results(contrast = c("group","Col-0_1000uM","Col-0_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("group","Col-0_1000uM","Col-0_50uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

#Phf1
df_temp <- dds %>% results(contrast = c("group","phf1_50uM","phf1_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

df_temp <- dds %>% results(contrast = c("group","phf1_1000uM","phf1_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("group","phf1_1000uM","phf1_50uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

#PP
df_temp <- dds %>% results(contrast = c("group","phr1/phl1_50uM","phr1/phl1_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]

df_temp <- dds %>% results(contrast = c("group","phr1/phl1_1000uM","phr1/phl1_0uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]


df_temp <- dds %>% results(contrast = c("group","phr1/phl1_1000uM","phr1/phl1_50uM")) %>%
  as.data.frame
df_temp[rownames(df_temp) %in% c("Root; k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Burkholderia"),]
