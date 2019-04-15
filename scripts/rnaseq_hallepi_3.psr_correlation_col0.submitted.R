library(ohchibi)
library(DESeq2)
library(paletteer)

#Set random seed
set.seed(130816)


source('plotting_parameters_hallepi.R')

Dat_rnaseq <- readRDS(file = "../cleandata/dat_rnaseq_hallepi.RDS")
Dat <- Dat_rnaseq$Dat_rnaseq_halle


#Subset dataset
Dat <- Dat %>% subset.Dataset(Fraction == "Root" & Genotype == "Col_0",
                              drop = T,clean = T)

# Create deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,colData = Dat$Map,
                              design = ~Phosphate)

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


#Calculate the mean per gene so we dont overinflate the points
melted_psr <- dcast(data = melted_vsd,formula = Gene ~ Phosphate,
                    fun.aggregate = mean,value.var = "value") %>%
  melt
colnames(melted_psr)[2] <- "Phosphate"



#Load phenotypic data 
df <- read.table("../rawdata/phenotype.csv",header = T,sep = ",",
                 quote = "",comment.char = "")

df$Phosphate <- df$Soil %>% 
  gsub(pattern = "^[a-c]-",replacement = "",perl = T) %>%
  gsub(pattern = "\\+",replacement = "_") %>%
  factor(levels = c("low","medium","high","low_Pi"))

df$PhosphateNum <- rep(0,nrow(df))
df$PhosphateNum[which(df$Phosphate == "medium")] <- "15"
df$PhosphateNum[which(df$Phosphate == "high")] <- "45"
df$PhosphateNum <- df$PhosphateNum %>% as.numeric

#Readjust some levels
df$Genotype <- df$Genotype %>% gsub(pattern = "c",replacement = "C") %>%
  gsub(pattern = "phr/phl",replacement = "phr1_phl1") %>%
  gsub(pattern = "-",replacement = "_") %>%
  gsub(pattern = "phf",replacement = "phf1") %>% 
  factor(levels = c("Col_0","phf1","phr1_phl1"))

df <- df %>% subset(Genotype == "Col_0") %>% droplevels

df$ZPiperG <- df$PiperG %>% scale

melted_pi <- melt(data = df,id.vars = "Phosphate",measure.vars = "ZPiperG")
melted_psr <- data.frame(Phosphate = melted_psr$Phosphate,
                         variable = rep("Regulon",nrow(melted_psr)),
                         value = melted_psr$value
)

merged <- rbind(melted_psr,melted_pi)

#Load palette for phosphate gradient
palette_pi_soil <- paletteer_d(package = "awtools",
                               palette = "mpalette",direction = -1)[1:4]
names(palette_pi_soil) <- c("low","medium","high","low_Pi")

p_raw <- chibi.boxplot(Map = merged,facet_formula = "variable",
                       x_val = "Phosphate",y_val = "value",col_val = "Phosphate",
                       mpalette = palette_pi_soil,
                       stroke_point = 0,median_colored_as_points = T,
                       size_point = 0,size_median = size_median) +
  theme_hallepi_boxplot + 
  geom_sina(alpha = 0.15, color = median_color,stroke = 0.5,size = 6) +  scale_size(range=c(3,6)) + 
  ylab("z-score") + theme(legend.position = "none") +
  geom_vline(xintercept = c(1.5,2.5,3.5),size = size_vline, color = color_vline)
oh.save.pdf(p = p_raw,outname = "figure1_rnaseq_hallepi_col0_regulon_pi.pdf",outdir = "../figures")

#CReate figure alone 
temp <- merged
p_alone <- temp %>% subset(variable == "Regulon") %>% droplevels %>%
  chibi.boxplot(Map = .,
                x_val = "Phosphate",y_val = "value",col_val = "Phosphate",
                mpalette = palette_pi_soil,
                stroke_point = 0,median_colored_as_points = T,
                size_point = 0,size_median = size_median) +
  theme_hallepi_boxplot + 
  geom_sina(alpha = 0.15, color = median_color,stroke = 0.5,size = 6) +  scale_size(range=c(3,6)) + 
  ylab("z-score") + theme(legend.position = "none") +
  geom_vline(xintercept = c(1.5,2.5,3.5),size = size_vline, color = color_vline)
saveRDS(object = p_alone,file = "../cleandata/fig1_zscore_193_col0.RDS")
#Correlation 
df_cor <- dcast(data = merged,fun.aggregate = mean,
                value.var = "value",Phosphate ~ variable) 

#Taken from 
#https://stackoverflow.com/questions/7549694/adding-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
mdf <- df_cor
colnames(mdf) <- c("Phosphate","y","x")

p_cor <- ggplot(df_cor, aes( ZPiperG, Regulon , fill = Phosphate)) + 
  geom_vline(xintercept = 0,size = 0.6, color = "#D9D9D9",linetype = "longdash") + 
  geom_hline(yintercept = 0,size = 0.6, color = "#D9D9D9",linetype = "longdash") +
  geom_smooth(data = df_cor, aes(ZPiperG,Regulon),inherit.aes = F,method = "lm",
              color = "#414141",alpha = 0.2) +
  geom_point(shape = 21,size = 25) + 
  theme_ohchibi(size_axis_text.x = size_axis_text.y,
                size_axis_text.y = size_axis_text.y,
                size_axis_title.x = size_axis_title.y,
                size_axis_title.y = size_axis_title.y,font_family = "Arial") + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),legend.position = "none") + 
  scale_fill_manual(values = palette_pi_soil) + ylab("z-score (PSR Regulon)") +
  xlab(label = "z-score (Pi Accumulation)") 

#geom_text(x = 0.3, y = 0.3, label = lm_eqn(mdf), parse = TRUE,family = "Arial",size = 20)
oh.save.pdf(p = p_cor,outname = "figure1_correlation_rna_pi_col0.pdf",outdir = "../figures/")

##Keep the information of the correlation
m1 <- lm(df_cor$Regulon ~ df_cor$ZPiperG)
m1_cor <- cor.test(df_cor$Regulon,df_cor$ZPiperG)
capture.output(file = "../figures/figure1_correlation_rna_pi_col0_summary.doc",
               append = F,print(summary(m1)))
capture.output(file = "../figures/figure1_correlation_rna_pi_col0_summary.doc",
               append = T,print((m1_cor)))
