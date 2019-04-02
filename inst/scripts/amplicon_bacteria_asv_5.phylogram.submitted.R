library(ohchibi)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf")


source('plotting_parameters_hallepi.R')

#Load dataset
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
Dat_rar <- Dat$Rarefied

#Remove bulksoil
Dat_rar <- Dat_rar %>% 
  subset.Dataset(subset = Fraction != "BulkSoil",drop = T,clean = T)

#Collapse taxonomy
Dat_phyla <- Dat_rar %>% collapse_by_taxonomy.Dataset(Dat = .,level = 3)

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
Other <- Tab[which(!(rownames(Tab) %in% mphyla )),] %>%
  colSums
Tab <- rbind(Tab_pr,Other)

#Create new phyla dataset
Dat_phyla <- create_dataset(Tab = Tab,Map = Map)


#Call phylogram
res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "Fraction",
                       size_ticks_x = 0,size_strip_text = 35,size_axis_text = 25,
                       legend_proportion_size = legend_proportion_size,size_legend_text = 30,
                       size_axis_title = 0,font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_phyla() +
  theme(aspect.ratio = 10) + theme(legend.position = "none")
oh.save.pdf(p = p, outdir = "../figures/",outname = "figure2_bacteria_asv_phylogram_fraction.pdf")

res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "Fraction+Genotype",
                       size_ticks_x = 0,size_strip_text = 22,size_axis_text = 25,
                       ,size_axis_title = 0,
                       legend_proportion_size = legend_proportion_size,size_legend_text = 30,
                       font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_phyla()+
  theme(aspect.ratio = 10) + theme(legend.position = "none")
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_bacteria_asv_phylogram_fraction_by_genotype.pdf")




res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "Fraction+Phosphate",
                       size_ticks_x = 0,size_strip_text = 22,size_axis_text = 25,
                       ,size_axis_title = 0,
                       legend_proportion_size = 3,size_legend_text = 30,
                       font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_phyla()+
  theme(aspect.ratio = 10) + theme(legend.position = "none")
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_bacteria_asv_phylogram_fraction_by_phosphate.pdf")



res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "Fraction+Genotype+Phosphate",
                       size_ticks_x = 0,size_strip_text = 8,size_axis_text = 25,
                       ,size_axis_title = 0,
                       legend_proportion_size = 0,size_legend_text = 0,
                       font_family = "Arial")
p <- res$p_mean +  scale_fill_phyla()+
  theme(aspect.ratio = 20) + theme(legend.position = "none")
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_bacteria_asv_phylogram_fraction_by_genotype_by_phosphate.pdf")


