library(ohchibi)
library(palettesPM)
library(paletteer)
library(corrplot)
library(egg)
library(emmeans)


source('plotting_parameters_hallepi.R')

df <- read.table("../rawdata/phenotype.csv",header = T,sep = ",",
                 quote = "",comment.char = "")

df$Phosphate <- df$Soil %>% 
  gsub(pattern = "^[a-c]-",replacement = "",perl = T) %>%
  factor(levels = c("low","medium","high","low+Pi"))

df$PhosphateNum <- rep(0,nrow(df))
df$PhosphateNum[which(df$Phosphate == "medium")] <- "15"
df$PhosphateNum[which(df$Phosphate == "high")] <- "45"
df$PhosphateNum <- df$PhosphateNum %>% as.numeric

#Readjust some levels
df$Genotype <- df$Genotype %>% gsub(pattern = "c",replacement = "C") %>%
  gsub(pattern = "phr/phl",replacement = "phr1/phl1") %>%
  gsub(pattern = "phf",replacement = "phf1") %>% 
  factor(levels = c("Col-0","phf1","phr1/phl1"))


colnames(df)[which(colnames(df) == "plot" )] <- "Plot"
df$Plot <- df$Plot %>% paste0("Plot",.) %>% factor

m1 <- aov(formula = PiperG ~ Genotype + Phosphate,data = df)
emmeans(object = m1,specs = "Genotype") %>% CLD

paleta <- palettesPM::pm.colors.genotypes()
p <- chibi.boxplot(Map = df,x_val = "Genotype",y_val = "PiperG",col_val = "Genotype",
              ,alpha_point = 1,style = "mix",mpalette = paleta,
              size_point = 20,size_median = size_median,median_colored_as_points = TRUE) +
  theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none")+ scale_y_continuous(limits = c(0,75))

oh.save.pdf(p = p,outname = "supfigure1_phf1.pdf",outdir =  "../figures/",height = 20,width = 20)
