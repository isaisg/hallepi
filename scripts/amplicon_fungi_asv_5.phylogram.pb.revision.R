library(ohchibi)
library(extrafont)
loadfonts(device = "pdf")

#Done
legend_proportion_size <- 4

#Load dataset
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$Rarefied

#Remove bulksoil
Dat_rar <- Dat_rar %>% 
  subset.Dataset(subset = Fraction != "BulkSoil",drop = T,clean = T)

#Collapse taxonomy
Dat_order <- Dat_rar %>% collapse_by_taxonomy.Dataset(Dat = .,level = 5)

#Split Tab and Map to create another structure
Tab  <- Dat_order$Tab
Map <- Dat_order$Map



torders <- Tab %>% rownames %>%
  strsplit(split = "\\;") %>% unlist %>%
  grep(pattern = "o__",value = T) %>% 
  gsub(pattern = "o__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")

tphylum <- Tab %>% rownames %>%
  strsplit(split = "\\;") %>% unlist %>%
  grep(pattern = "p__",value = T) %>% 
  gsub(pattern = "p__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")


df_ins <- data.frame(Phylum = tphylum, Order = torders)
rownames(Tab) <- df_ins$Order

df_ab <- Tab %>% rowSums %>%
  data.frame(Order = names(.), Freq = .,row.names = NULL)
df_ab <- merge(df_ab,df_ins,by = "Order")
df_ab <- with(df_ab,order(-Freq)) %>% df_ab[.,]

chosen_orders <- df_ab$Order[1:19] %>% as.character
Other <- Tab[which(!(rownames(Tab) %in% chosen_orders)),] %>% colSums
Tab_chosen <- Tab[which((rownames(Tab) %in% chosen_orders)),]
Tab <- rbind(Tab_chosen,Other)


df_other <- data.frame(Phylum = "Other", Order = "Other",
                       Freq =which(!(df_ab$Order %in% chosen_orders)) %>%
                         df_ab[.,] %$% Freq %>% sum)
df_ab <- which(df_ab$Order %in% chosen_orders) %>% df_ab[.,] %>%
  rbind(.,df_other)
df_ab$Phylum <- df_ab$Phylum %>% factor(levels = c("Ascomycota","Basidiomycota",
                                                   "Mortierellomycota","Fungi_unclassified","Other"))
df_ab <- with(df_ab,order(Phylum,Order)) %>%
  df_ab[.,] %>% droplevels

Tab <- match(df_ab$Order,rownames(Tab)) %>% Tab[.,]

#Create new phyla dataset
Dat_order <- create_dataset(Tab = Tab,Map = Map)

# paleta <- readRDS(file = "../rawdata/palette_fungi_orders.RDS")
# noms <- paleta %>% names
# mords <- Dat_order$Tab %>% rownames
# 
# conserved <- intersect(noms,mords)
# paleta_kept <- names(paleta) %in% conserved %>%
#   paleta[.]
# 
# paleta_kept <- names(paleta) %in% conserved %>%
#   paleta[.]
# 
# paleta_missing <- paleta[!(names(paleta) %in% conserved) ]
# names(paleta_missing) <- mords[!(mords %in% conserved)]
# 
# paleta <- c(paleta_kept,paleta_missing)
# paleta <- match(rownames(Tab),names(paleta)) %>%
#   paleta[.]
# saveRDS(object = paleta,file = "../rawdata/palette_fungi_orders.revisited.RDS")
paleta <- readRDS(file = "../rawdata/hallepi_palette_fungi_orders.RDS")



#Call phylogram
res <- chibi.phylogram(Tab = Dat_order$Tab,Map = Dat_order$Map,facet_formula = "Fraction",
                       size_ticks_x = 0,size_strip_text = 35,size_axis_text = 25,
                       legend_proportion_size = 1,size_legend_text = 30,
                       size_axis_title = 0,font_family = "Arial",size_title_text = 35)

#Save figure for legends
p_leg <- res$p_mean + 
  scale_fill_manual(values = paleta) 
  oh.save.pdf(p = p_leg, outdir = "../figures/",
              outname = "figure2_fungi_asvs_phylogram_fraction.legend.pdf")


p <- res$p_mean + 
  scale_fill_manual(values = paleta) +
  theme(aspect.ratio = 10) + theme(legend.position = "none") + 
  scale_y_continuous(c(0,0.25,0.5,0.75,1))
oh.save.pdf(p = p, outdir = "../figures/",outname = "figure2_fungi_asvs_phylogram_fraction.pdf")

write.table(x = res$p_mean$data,file = "../data_figures/data_Fig2F_census.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)





res <- chibi.phylogram(Tab = Dat_order$Tab,Map = Dat_order$Map,facet_formula = "Fraction+Genotype",
                       size_ticks_x = 0,size_strip_text = 20,size_axis_text = 25,
                       size_axis_title = 0,
                       legend_proportion_size = 3,size_legend_text = 30,
                       font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_manual(values = paleta) +
  theme(aspect.ratio = 10) + theme(legend.position = "none") + 
  scale_y_continuous(c(0,0.25,0.5,0.75,1))
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_fungi_asvs_phylogram_fraction_by_genotype.pdf")



res <- chibi.phylogram(Tab = Dat_order$Tab,Map = Dat_order$Map,facet_formula = "Fraction+Phosphate",
                       size_ticks_x = 0,size_strip_text = 18,size_axis_text = 25,
                       size_axis_title = 0,
                       legend_proportion_size = 3,size_legend_text = 30,
                       font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_manual(values = paleta) +
  theme(aspect.ratio = 10) + theme(legend.position = "none") + 
  scale_y_continuous(c(0,0.25,0.5,0.75,1))
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_fungi_asvs_phylogram_fraction_by_phosphate.pdf")



res <- chibi.phylogram(Tab = Dat_order$Tab,Map = Dat_order$Map,facet_formula = "Fraction+Genotype+Phosphate",
                       size_ticks_x = 0,size_strip_text = 8,size_axis_text = 25,
                       size_axis_title = 0,
                       legend_proportion_size = 0,size_legend_text = 0,
                       font_family = "Arial")
p <- res$p_mean + scale_fill_manual(values = paleta) +
  theme(aspect.ratio = 20) + theme(legend.position = "none") + 
  scale_y_continuous(c(0,0.25,0.5,0.75,1))
oh.save.pdf(p = p,outdir = "../figures/",outname = "figure2_fungi_asvs_phylogram_fraction_by_genotype_by_phosphate.pdf")


write.table(x = res$p_mean$data,file = "../data_figures/data_S2F.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


