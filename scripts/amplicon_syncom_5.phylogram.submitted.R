library(ohchibi)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf") 


set.seed(130816)
## Abundance
Dat <- readRDS(file = "../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_rab <- Dat$RelativeAbundance
#Collapse taxonomy
Dat_phyla <- Dat_rab %>% collapse_by_taxonomy.Dataset(Dat = .,level = 3)

#Split Tab and Map to create another structure
Tab  <- Dat_phyla$Tab
Map <- Dat_phyla$Map

rownames(Tab) <- Tab %>% rownames %>%
  strsplit(split = "\\;") %>% unlist %>%
  grep(pattern = "p__",value = T) %>% 
  gsub(pattern = "p__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")

mphyla <- palettesPM::pm.names.phyla()

Tab_pr <- Tab[which((rownames(Tab) %in% mphyla )),]
Other <- Tab[which(!(rownames(Tab) %in% mphyla )),] 
Tab <- rbind(Tab_pr,Other)

#Create new phyla dataset
Dat_phyla <- create_dataset(Tab = Tab,Map = Map)

Dat_phyla$Map$typebyTissue <- Dat_phyla$Map$typebyTissue %>%
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  factor(levels = c("Inoculum","Agar_NP","Agar_SC","Root_SC","Shoot_SC"))
#Call phylogram
res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "typebyTissue",
                       size_ticks_x = 0,size_strip_text = 35,size_axis_text = 25,
                       legend_proportion_size = 4,size_legend_text = 30,
                       size_axis_title = 0,font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_phyla() +
  theme(aspect.ratio = 10) + theme(legend.position = "none")
oh.save.pdf(p = p, outdir = "../figures/",outname = "figure4_syncom_phylogram_fraction.pdf")



#Write this dataset
write.table(x = res$p_mean$data,file = "../data_figures/data_Fig4B.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



