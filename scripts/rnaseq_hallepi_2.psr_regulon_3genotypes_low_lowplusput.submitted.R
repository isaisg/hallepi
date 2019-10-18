library(ohchibi)
library(DESeq2)
library(paletteer)

#Set random seed
set.seed(130816)

#Done

source('plotting_parameters_hallepi.R')


Dat_rnaseq <- readRDS(file = "../cleandata/dat_rnaseq_hallepi.RDS")
Dat <- Dat_rnaseq$Dat_rnaseq_halle


#Subset dataset
Dat <- Dat %>% subset.Dataset(Fraction == "Root" &
                  (Phosphate == "low" | Phosphate == "low_Pi"),
                              drop = T,clean = T)

# Create deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,colData = Dat$Map,
                              design = ~ Plot + group)
vsd <- vst(object = dds,blind = FALSE)

#Melt object and z-scale it
melted_vsd <- vsd %>% assay %>% t  %>% scale %>% t %>% melt

colnames(melted_vsd) <- c("Gene","sample","value")

melted_vsd <- merge(melted_vsd,Dat$Map, by = "sample")

#Low the psr genes
regulon_psr <- read.table("../rawdata/regulon_psr.csv",header = F) %$%
  V1 %>% as.character

melted_vsd <- which(melted_vsd$Gene %in% regulon_psr) %>%
  melted_vsd[.,] %>% droplevels

melted_vsd$group <- paste(melted_vsd$Fraction,melted_vsd$Phosphate,melted_vsd$Genotype,sep = "_")

#Calculate the mean per gene so we dont overinflate the points
melted_psr <- dcast(data = melted_vsd,formula = Gene ~ group,fun.aggregate = mean,value.var = "value") %>%
  melt
colnames(melted_psr)[2] <- "group"

melted_psr <- melted_vsd[,c(6:8,10)] %>% unique %>% merge(melted_psr, . , by = "group",all.x = TRUE)

paleta <- palette_pi_soil[c(1,4)]

#Write dataset
write.table(x = melted_psr,file = "../data_figures/data_S1E.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



p <- chibi.boxplot(Map = melted_psr,style = "mix",
                   facet_formula = "Genotype",x_val = "Phosphate",y_val = "value",
                   col_val = "Phosphate",size_point = 0,stroke_point = 0,
                   median_colored_as_points = TRUE,
                   size_median = size_median,mpalette = paleta) + theme_hallepi_boxplot +
  geom_vline(xintercept = c(1.5),size = size_vline, color = color_vline) +
  #geom_text(data = df_test,
  #             aes(x = Soil,y = 2.7,label = Symbol,size = Size),
  #              inherit.aes = F,family ="Arial",color = "#414141")  +
  geom_sina(alpha = 0.15, color = median_color,
            stroke = 0.5,size = 6) +  scale_size(range=c(3,6)) + 
  ylab("z-score") + theme(legend.position = "none")
oh.save.pdf(p = p,outname = "figure1_rnaseq_hallepi_psr_regulon.pdf",outdir = "../figures")
