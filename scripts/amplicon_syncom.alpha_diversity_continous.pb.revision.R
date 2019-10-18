library(ohchibi)
library(emmeans)
library(paletteer)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf")


palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("typebyTissue","Genotype","condition",
                             "typebyTissue:condition","typebyTissue:Genotype","Residual")


set.seed(130816)

source(file = "plotting_parameters_hallepi.R")

Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_rar <- Dat_ori$RelativeAbundance

#Palette

colores_phosphate <- c("black","#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")

#Calculate Shannon and Richness
Dat_rar$Map$Shannon <- vegan::diversity(x = Dat_rar$Tab, index = "shannon", MARGIN = 2 )
Dat_rar$Map$Richness <- colSums(Dat_rar$Tab > 0)

#Readjust levels
Dat_rar$Map$typebyTissue <- Dat_rar$Map$typebyTissue %>% 
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "Agar_NP",replacement = "AgarNoPlant") %>% 
  gsub(pattern = "Agar_SC",replacement = "AgarPlant") %>%
  gsub(pattern = "Root_SC",replacement = "Root") %>% 
  gsub(pattern = "Shoot_SC",replacement = "Shoot") %>% 
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))

Mapa <- Dat_rar$Map %>% subset(typebyTissue != "Inoculum") %>%
  droplevels


###Continous version of the analysis
Mapa$gradient <- Mapa$condition %>% 
  gsub(pattern = " micromolar.*",replacement = "") %>% 
  as.numeric
Mapa$LPi <- (Mapa$gradient +1) %>% log10
Mapa$LPi2 <- Mapa$LPi^2

#Test lineal model per fraction
Res <- NULL
for( fraction in Mapa$typebyTissue %>% unique){
  Mapa_sub <- Mapa %>% subset(typebyTissue == fraction) %>% droplevels
  m1 <- lm(formula = Shannon ~ LPi + Rep,data = Mapa_sub)
  mcor <- cor.test(Mapa_sub$Shannon,Mapa_sub$LPi)
  df <- summary(m1)  %$% coefficients %>% as.data.frame
  df <- df[2,]
  label <- paste(round(m1$coefficients[1],2),round(m1$coefficients[2]  ,3),round(summary(m1)$r.squared,3))
  df <- data.frame(typebyTissue = fraction,estimate = df$Estimate,LM_p.value =round (df$`Pr(>|t|)`,6),
                   cor = mcor$estimate,Cor_p.value  = mcor$p.value,Label = label)
  Res <- rbind(Res,df)
  
}
rownames(Res) <- NULL
Res$Label <- Res$Label %>% as.character

#Lineal
p <- ggplot(data = Mapa,aes(LPi,Shannon)) +
  stat_smooth(method = "lm",
              size = 4,color = "black" )+ 
  facet_grid(.~typebyTissue,space = "free",scales = "free")  +
  geom_point(aes(fill = condition),shape = 21,size = 10 ) + 
  scale_fill_manual(values = colores_phosphate[2:7])+
  theme_ohchibi() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(family = "Arial",face = "bold",size = 30)
  ) + xlab(label = "Phosphate") + 
  ylab(label = "Shannon Diversity Index") +
  geom_text(data = Res,aes(x = 1,y = 3.1,label = LM_p.value),size = 6)
oh.save.pdf(p = p,outname = "alpha_diversity_syncom_continous.pdf",outdir = "../figures/")
