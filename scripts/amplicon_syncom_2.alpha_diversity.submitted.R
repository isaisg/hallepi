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
#colores_phosphate <- RColorBrewer::brewer.pal(n=6,"Oranges")
#colores_phosphate <- c("black",colores_phosphate)
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

#Write supplemental table 9
Mapa <- Dat_rar$Map
rownames(Mapa) <- NULL
Mapa <- Mapa[,-c(1,2,10:12)]
Mapa <- Mapa[,c(2,7,8,4,3,5,1,6,9,10)] 
Mapa <- Mapa[,-3]
write.table(x = Mapa,file = "../cleandata/sup_table_9.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)

m1 <- Dat_rar$Map %>% subset(typebyTissue != "Inoculum") %>%
  subset(typebyTissue != "AgarNoPlant") %>%
  droplevels %>% aov(formula = Shannon ~ typebyTissue + condition  + typebyTissue:condition + Rep)
anova(m1)

p_perm <- chibi.anova(m1)+ scale_fill_manual(values = palette_variance) + 
  theme(legend.position = "none")


#df <- Dat_rar$Map %>% subset(typebyTissue != "Inoculum") %>% droplevels
#m1 <- lm(formula = Shannon  ~ typebyTissue + Pi  +  typebyTissue:Pi + Rep,data = df)
#anova(m1)

##Compute statistics
Res_em <- Dat_rar$Map %>% aov(formula = Shannon ~ typebyTissue  + Rep)%>% 
  emmeans(.,specs = "typebyTissue") %>%
  CLD %>% as.data.frame

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

pal_frac <- palettesPM::pm.colors.fractions()
pal_frac <- match(Dat_rar$Map$typebyTissue %>% levels,names(pal_frac)) %>%
  pal_frac[.]

p <- chibi.boxplot(Map = Dat_rar$Map,x_val = "typebyTissue",
                   y_val = "Shannon",col_val = "typebyTissue",
                   size_median = size_median,style = "mix",
                   size_point = 0,alpha_point = 0,stroke_point = 0,
                   median_colored_as_points = T,mpalette = pal_frac,
                   size_axis_text.x = 30,
                   size_axis_text.y = size_axis_text.y,
                   size_axis_title.x = size_axis_title.x,
                   size_axis_title.y = size_axis_title.y,
                   size_legend_text = size_legend_text,
                   strip_text_size = strip_text_size,
                   size_title_text = size_title_text,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), size = 0.8, color = "#D9D9D9") +
  ylab(label = "Shannon Diversity Index") 
p <- p + geom_sina(aes(fill = typebyTissue), size = 8.5,pch = 21,alpha = 0.3)

p <- p + geom_text(data = Res_em,
                   aes(x = typebyTissue,y = 3.8,label = Letters),
                   inherit.aes = F,size = size_anova,family ="Arial",color = "#414141")  +
  theme(legend.position = "none")

composition <- egg::ggarrange(p_perm,p,widths = c(0.05,1))
oh.save.pdf(p = p,outname = "figure5_syncom_alphadiversity_all.pdf",outdir = "../figures")

### Figure of phosphate  ###
#Test inside each fraction

temp_ua <- Dat_rar$Map %>% subset(typebyTissue == "AgarNoPlant") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition") %>%
  CLD %>% as.data.frame
temp_ua$typebyTissue <- rep("AgarNoPlant",nrow(temp_ua))
Res_ua <- Dat_rar$Map %>% subset(typebyTissue == "AgarNoPlant") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition") %>%as.data.frame
Res_ua$typebyTissue <- rep("AgarNoPlant",nrow(Res_ua))

temp_pa <- Dat_rar$Map %>% subset(typebyTissue == "AgarPlant") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition") %>%
  CLD %>% as.data.frame
temp_pa$typebyTissue <- rep("AgarPlant",nrow(temp_pa))
Res_pa <- Dat_rar$Map %>% subset(typebyTissue == "AgarPlant") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition")%>%as.data.frame
Res_pa$typebyTissue <- rep("AgarPlant",nrow(Res_pa))


temp_root <- Dat_rar$Map %>% subset(typebyTissue == "Root") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition") %>%
  CLD %>% as.data.frame
temp_root$typebyTissue <- rep("Root",nrow(temp_root))
Res_root <- Dat_rar$Map %>% subset(typebyTissue == "Root") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition")%>%as.data.frame
Res_root$typebyTissue <- rep("Root",nrow(Res_root))


temp_shoot <- Dat_rar$Map %>% subset(typebyTissue == "Shoot") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition") %>%
  CLD %>% as.data.frame
temp_shoot$typebyTissue <- rep("Shoot",nrow(temp_shoot))
Res_shoot <- Dat_rar$Map %>% subset(typebyTissue == "Shoot") %>% droplevels %>%
  aov(formula = Shannon ~ condition + Rep,data = .) %>% emmeans(.,specs = "condition")%>%as.data.frame
Res_shoot$typebyTissue <- rep("Shoot",nrow(Res_shoot))

#Bind all the results
Res_em <- rbind(temp_ua,temp_pa,temp_root,temp_shoot)

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

Res_em$typebyTissue <- Res_em$typebyTissue %>%
  factor(levels = c("AgarNoPlant","AgarPlant","Root","Shoot"))

Map <- Dat_rar$Map %>% subset(typebyTissue != "Inoculum") %>% droplevels

Res_pval <- Res_em
Res_em <- rbind(Res_ua,Res_pa,Res_root,Res_shoot)
colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
#colores_phosphate <- c("#fff7f3","#fde0dd","#fa9fb5","#dd3497","#7a0177","#49006a")

p <- ggplot(data = Res_em,mapping = aes(condition,emmean,fill = condition)) + 
  facet_grid(.~typebyTissue) +
  geom_pointrange(aes(y = emmean, ymin = lower.CL, ymax= upper.CL,color = condition), 
                  ,size = 3,alpha = 1) +
  geom_point(shape = 21,aes(color = condition),alpha = 1)  +
  geom_line(aes(y = emmean,group = 1),size =5,alpha = 1) +
  theme_ohchibi() + 
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = colores_phosphate)+ scale_fill_manual(values = colores_phosphate) +
  geom_text(data = Res_pval,mapping = aes(x = condition,y = 2.95,label = Letters),inherit.aes = F)
oh.save.pdf(p = p,outname = "figure6_syncom_alphadiversity_pigradient.pdf",outdir = "../figures")


#Write this dataset
write.table(x = Dat_rar$Map,file = "../data_figures/data_Fig5B_Fig5E_census.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)

