library(ohchibi)
library(palettesPM)
library(paletteer)
library(emmeans)
library(extrafont)
loadfonts(device = "pdf")
library(egg)

#Set random seed
set.seed(130816)



setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')
source('plotting_parameters_hallepi.R')


#Read bacterial  dat
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$Rarefied


#Calculate Shannon and Richness
Dat_rar$Map$Shannon <- vegan::diversity(x = Dat_rar$Tab, index = "shannon", MARGIN = 2 )
Dat_rar$Map$Richness <- colSums(Dat_rar$Tab > 0)
Dat_rar$Map$LDepth <- log(Dat_rar$Map$Depth)



Dat_rar$Map %>%
  write.table(x = .,file = "../cleandata/sup_table_3_fungi.csv",
              append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,subset = Genotype!="No_Plant",drop = T,clean = T) 

###Shannon### anova framework 
m1 <- aov(formula = Shannon ~ Fraction + Genotype + Phosphate + Fraction:Genotype ,
          data = Dat_sub$Map)
p_anova <- chibi.anova(m1) + scale_fill_manual(values = palette_variance) + 
  theme(legend.position = "none")

Res_em <- m1 %>% emmeans(pairwise ~Fraction:Genotype) %$% emmeans %>% CLD %>% as.data.frame 


Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"

#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>%
  strsplit(split = "") %>%
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>%
    unlist %>% toupper  %>% paste0(collapse = "")
}


#Rest across all 

#Define palette for fraction
paleta <- palettesPM::pm.colors.fractions()
paleta <- paleta[c(10,7,9)]

p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "Genotype",x_val = "Fraction",
                   y_val = "Shannon",col_val = "Fraction",
                   size_median = size_median,style = "mix",
                   size_point = 0,alpha_point = 0,stroke_point = 0,
                   median_colored_as_points = T,mpalette = paleta,
                   size_axis_text.x = size_axis_text.x,
                   size_axis_text.y = size_axis_text.y,
                   size_axis_title.x = size_axis_title.x,
                   size_axis_title.y = size_axis_title.y,
                   size_legend_text = size_legend_text,
                   strip_text_size = strip_text_size,
                   size_title_text = size_title_text,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index") + 
  # geom_vline(xintercept = c(1.5,2.5),size = size_vline,color = color_vline) +
  geom_sina(aes(fill = Fraction), size = size_sina,pch = 21,alpha = alpha_sina)
p <- p + geom_text(data = Res_em,
                   aes(x = Fraction,y = 4.3,label = Letters),
                   inherit.aes = F,size = size_anova,family ="Arial",color = median_color) 
p <- p + theme(legend.position = "none",
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5))

composition <- egg::ggarrange(p_anova,p,widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure2_fungi_asv_alphadiversity.pdf",outdir = "../figures")

#Make supplemental figures splitting by plot

p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "Fraction",x_val = "Genotype",
                   y_val = "Shannon",col_val = "Phosphate",
                   size_median = size_median,style = "mix",
                   size_point = 8,alpha_point = 0.8,stroke_point = 0,
                   median_colored_as_points = T,mpalette = palette_pi_soil,
                   size_axis_text.x = size_axis_text.x,
                   size_axis_text.y = size_axis_text.y,
                   size_axis_title.x = size_axis_title.x,
                   size_axis_title.y = size_axis_title.y,
                   size_legend_text = size_legend_text,
                   strip_text_size = strip_text_size,
                   size_title_text = size_title_text,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index") + 
  # geom_vline(xintercept = c(1.5,2.5),size = size_vline,
  #           color = color_vline) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5))

oh.save.pdf(p = p,outname = "supfigure2_fungi_asv_alphadiversity_extended.pdf",outdir = "../figures")

rm(list=ls())
dev.off()
