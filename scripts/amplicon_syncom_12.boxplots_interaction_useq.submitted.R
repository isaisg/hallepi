library(ohchibi)
library(RColorBrewer)

set.seed(130816)
source(file = "plotting_parameters_hallepi.R")


Dat <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat <- Dat$RelativeAbundance
Dat <- Dat %>% 
  subset.Dataset(subset = typebyTissue != "Inoculum_Inoculum",drop = T,clean = T)
Dat$Tab <- Dat$Tab[grep(pattern = "OTU",x = rownames(Dat$Tab),invert = T),]

chosen <- c("Sequence_1","Sequence_16","Sequence_30","Sequence_32",
            "Sequence_39","Sequence_67","Sequence_72","Sequence_87",
            "Sequence_97")
Dat$Tab <- Dat$Tab[which(rownames(Dat$Tab)%in%chosen),]
#colores_phosphate <- brewer.pal(n=6,"Oranges")
colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")

j <- 2
name<-rownames(Dat$Tab)[j]
#Plot by conditions
Mapa <-Dat$Map
Tab<-Dat$Tab
Mapa[, name] <- as.vector(Tab[name, ])


write.table(x = Dat$Tab,file = "../data_figures/data_S7BC.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


for(j in 1:length(rownames(Dat$Tab))){
  name<-rownames(Dat$Tab)[j]
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
    theme(legend.position = "none")
  
  #File name
  filename <- paste("figure_7a_condition_typebyTissue_interaction.",name,"pdf",sep="")
  oh.save.pdf(p = p,outname = filename,outdir = "../figures/")
}
