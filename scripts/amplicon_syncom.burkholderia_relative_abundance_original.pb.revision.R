library(ohchibi)
library(RColorBrewer)

set.seed(130816)
source(file = "plotting_parameters_hallepi.R")

Dat <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat <- Dat$RelativeAbundance
Dat <- Dat %>% 
  subset.Dataset(subset = typebyTissue != "Inoculum_Inoculum",drop = T,clean = T)


toremove <- Dat$Tax$ID %>% grep(pattern = "OTU",value = T)
Dat <- remove_taxons(Dat = Dat,taxons = toremove)
Dat <- clean(Dat)


Dat <- collapse_by_taxonomy.Dataset(Dat = Dat,level = 7)


colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")

name <- rownames(Dat$Tab) %>% grep(pattern = "g__Burkholderia",value = T)


#Plot by conditions
Mapa <-Dat$Map
Tab<-Dat$Tab
Mapa[, name] <- as.vector(Tab[name, ])
colnames(Mapa)[14] <- "USeq"
p <-  ggplot(data = Mapa,mapping = aes(condition,USeq,color = condition)) + facet_grid(.~typebyTissue) +
  stat_summary(fun.y = "mean", size = 5, geom = "point",
               position=position_nudge(x = 0.1, y = 0)) +
  stat_summary(fun.data = "mean_cl_normal", size = 2, 
               geom = "errorbar", width = 0, fun.args = list(conf.int = 0.9),
               position=position_nudge(x = 0.1, y = 0)) +
  scale_color_manual(values = colores_phosphate) + theme_ohchibi() +
  stat_summary(aes(y = USeq), fun.y=mean,
               colour="black", geom="line",group = "condition",size = 2) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),axis.text.x = element_blank())  +
  ylab(label = "Relative Abundance") +
  theme(legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 20)
        )

#File name
filename <- paste("figure_7a_burkholderia_genus_relative_abundance.pdf")
oh.save.pdf(p = p,outname = filename,outdir = "../figures/")


#Write the file

write.table(x = Mapa,file = "../data_figures/data_Fig6A.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


