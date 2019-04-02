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

### Plot for Phosphate
#Test the low,medium high in anova context categorical
df %>% subset(Phosphate != "low+Pi" & Genotype == "Col-0" ) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ Phosphate ) %>% summary

df %>% subset(Phosphate != "low+Pi" & Genotype == "phf1" ) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ Phosphate ) %>% summary

df %>% subset(Phosphate != "low+Pi" & Genotype == "phr1/phl1" & PiperG != 0) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ Phosphate ) %>% summary

#Test the low,medium high in anova context numeric
df %>% subset(Phosphate != "low+Pi" & Genotype == "Col-0" ) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ PhosphateNum ) %>% summary

df %>% subset(Phosphate != "low+Pi" & Genotype == "phf1" ) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ PhosphateNum ) %>% summary

df %>% subset(Phosphate != "low+Pi" & Genotype == "phr1/phl1" & PiperG != 0) %>%
  droplevels %>% 
  aov(formula = (PiperG) ~ PhosphateNum ) %>% summary

##Do a paired test  comparing low versus low+Pi
df_temp <- df %>% subset((Phosphate == "low+Pi"  | Phosphate == "low")& 
                           Genotype == "Col-0") %>% droplevels %>%
  acast(data = .,value.var = "PiperG",formula = Phosphate ~ Plot)
m1_col0 <- t.test(x = df_temp[1,],y = df_temp[2,],paired = TRUE)

df_temp <- df %>% subset((Phosphate == "low+Pi"  | Phosphate == "low")& 
                           Genotype == "phf1") %>% droplevels %>%
  acast(data = .,value.var = "PiperG",formula = Phosphate ~ Plot)
m1_phf1 <- t.test(x = df_temp[1,],y = df_temp[2,],paired = TRUE)

#Create segment df
dfseg <- data.frame(Genotype = c("Col-0","phf1"),
                    x = c(1,1),xend = c(4,4),
                    y = c(74,74), yend = c(74,74)
)
dfsig <- data.frame(Genotype = c("Col-0","phf1"),
                    x = c(2.5,2.5),
                    y = c(75,75), label = c("*","*")
)


#Create figure
p <- chibi.boxplot(Map = df,x_val = "Phosphate",y_val = "PiperG",col_val = "Phosphate",
                   facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                   size_point = size_point,size_median = size_median,median_colored_as_points = TRUE) +
  theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none")
p <- p +geom_segment(aes(x = x, y = y, xend = xend, yend = yend), 
                     data = dfseg,size = 2,color = median_color) +
  geom_text(mapping = aes(x = x , y = y, label = label),data = dfsig,
            size = 30,color = median_color,family = "Arial")
oh.save.pdf(p = p,outname = "figure1_phenotype_shootpi.pdf",outdir = "../figures")


#Create a split version of this figure
p_wt <- df %>% subset(Genotype == "Col-0") %>% droplevels %>% 
  chibi.boxplot(Map = .,x_val = "Phosphate",y_val = "PiperG",col_val = "Phosphate",
                facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                size_point = size_point,size_median = size_median,median_colored_as_points = TRUE) +
  theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none") + scale_y_continuous(limits = c(0,75))

p_single <- df %>% subset(Genotype == "phf1") %>% droplevels %>% 
  chibi.boxplot(Map = .,x_val = "Phosphate",y_val = "PiperG",col_val = "Phosphate",
                facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                size_point = size_point,size_median = size_median,median_colored_as_points = TRUE) +
  theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none")+ scale_y_continuous(limits = c(0,75))

p_double <- df %>% subset(Genotype == "phr1/phl1") %>% droplevels %>% 
  chibi.boxplot(Map = .,x_val = "Phosphate",y_val = "PiperG",col_val = "Phosphate",
                facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                size_point = size_point,size_median = size_median,median_colored_as_points = TRUE) +
  theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none")+ scale_y_continuous(limits = c(0,75))

#Read the rnaseq figure
p_rna <- readRDS(file = "../cleandata/fig1_zscore_193_col0.RDS")

#Create composition
composition <- egg::ggarrange(nrow = 2,ncol = 3,byrow = T,
                              p_wt,p_single,p_double,
                              p_rna)
oh.save.pdf(p = composition,outname = "figure1_phenotype_shootpi_composition.pdf",outdir = "../figures",width = 22,height = 25)



#Plot for weight
#Test inside each genotype
lala <- aov(formula = normWeight ~ Genotype +Phosphate +Genotype:Phosphate,data = df) %>%
  emmeans(pairwise ~ Genotype:Phosphate) %>% CLD
lala <- NULL
for(geno in df$Genotype %>% levels){
  lal <- df %>% subset(Genotype == geno) %>% droplevels %>%
    aov(formula = normWeight ~ Phosphate,data = .) %>%
    emmeans(specs = "Phosphate") %>% CLD
  lal$Genotype <- rep(geno,nrow(lal))
  lala <- rbind(lala,lal)
  
}
p <- chibi.boxplot(Map = df,x_val = "Phosphate",y_val = "normWeight",col_val = "Phosphate",
                   facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                   size_point = size_point,size_median = size_median,median_colored_as_points = TRUE,
) + theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  theme(legend.position = "none")+
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) 
oh.save.pdf(p = p,outname = "figure1_phenotype_shootweight.pdf",outdir = "../figures")


## Load shoot size information 
df_shoot <- read.table("../rawdata/Shoot_area_Halle.csv",header = T,sep = "\t",
                       quote = "",comment.char = "")
df_shoot$Genotype <- df_shoot$Genotype %>% 
  gsub(pattern = "col-0",replacement = "Col-0") %>%
  gsub(pattern = "phf",replacement = "phf1") %>% 
  gsub(pattern = "phr/phl",replacement = "phr1/phl1") %>%
  factor(levels = c("Col-0","phf1","phr1/phl1"))
df_shoot$Phosphate <- df_shoot$Phosphate %>% as.character %>%
  gsub(pattern = "low\\+pi",replacement = "low\\+Pi") %>%
  factor(levels = c("low","medium","high","low+Pi"))
df_shoot$Plot <- df_shoot$plot %>% paste0("Plot",.) %>% factor


df_shoot %>% aov(formula = TotSurfArea_cm2~Genotype+Phosphate + Genotype:Phosphate,data = .) %>%
  summary
lala <- aov(formula = TotSurfArea_cm2 ~ Genotype +Phosphate +Genotype:Phosphate,data = df_shoot) %>%
  emmeans(pairwise ~ Genotype:Phosphate) %>% CLD

p <- chibi.boxplot(Map = df_shoot,x_val = "Phosphate",y_val = "TotSurfArea_cm2",col_val = "Phosphate",
                   facet_formula = "Genotype",mpalette = palette_pi_soil,alpha_point = 1,
                   size_point = size_point,size_median = size_median,median_colored_as_points = TRUE,
) + theme_hallepi_boxplot + 
  ylab(label = "Phosphate/Gram") + xlab(label = "") + 
  geom_vline(xintercept = c(1.5,2.5,3.5), size = size_vline , color = color_vline) +
  theme(legend.position = "none")
oh.save.pdf(p = p,outname = "figure1_phenotype_shootarea.pdf",outdir = "../figures")


### Correlation
df$Uid <- paste(df$Genotype,df$Phosphate,df$Plot,df$plantnum,sep = "_")
df_shoot$Uid <- paste(df_shoot$Genotype,df_shoot$Phosphate,df_shoot$Plot,df_shoot$plant,sep ="_")

merged <- merge(df,df_shoot, by = "Uid") 
Tab <- merged[,c(6,10,17,18,19)] %>% as.matrix
rownames(Tab) <- merged$Uid
Tab_cor <- cor(Tab)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
cairo_pdf(filename = "../figures/figure1_correlation_phenotypes.pdf",
          width = 20,height = 12,family = "Arial")
corrplot(Tab_cor, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = 3,cl.cex = 2.5,tl.cex = 2.5,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()

#Write Supplemental Table 1
df <- merged[,c(1,2,4,5,6,7,10,11,13,17,18,19)]
df <- df[,c(1,2,4,9,3,8,5,6,7,10,11,12)]
colnames(df) <- c("UId","Pot","Plot","Plant","Genotype","PiSoil","PiperGram","Weight","NormWeight","TotalLength_cm",
                  "TotalProjArea_cm2","TotalSurfArea_cm2")
write.table(x = df,file = "../cleandata/sup_table_1.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)
