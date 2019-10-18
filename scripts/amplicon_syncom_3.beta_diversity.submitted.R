library(ohchibi)
library(palettesPM)
library(paletteer)

size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4


Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_rar <- Dat_ori$RelativeAbundance

colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
#colores_phosphate <- c("#fff7f3","#fde0dd","#fa9fb5","#dd3497","#7a0177","#49006a")


Dat_rar$Map$typebyTissue <- Dat_rar$Map$typebyTissue %>% as.character %>%
  gsub(pattern = "Agar_NP",replacement = "AgarNoPlant") %>%
  gsub(pattern = "Agar_SC",replacement = "AgarPlant") %>%
  gsub(pattern = "Root_SC",replacement = "Root") %>%
  gsub(pattern = "Shoot_SC",replacement = "Shoot") %>%
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))

pal_frac <- palettesPM::pm.colors.fractions()
pal_frac <- match(Dat_rar$Map$typebyTissue %>% levels,names(pal_frac)) %>%
  pal_frac[.]

#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")


mpco <- oh.pco(Tab = t(Dat_rar$Tab),Map = Dat_rar$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "ID_Matrix")


palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("typebyTissue","Genotype","condition",
                             "typebyTissue:condition","typebyTissue:Genotype","Residual")

#Permanova
Tab_bray <- distfun(t(Dat_rar$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition + typebyTissue:condition ,
                      data = Dat_rar$Map,strata = Dat_rar$Map$Rep,permutations = 10000)

capture.output(file = "../figures/figure5_syncom_betadiversity_summary.doc",
               append = F,print(mypermanova))



p_perm <- chibi.permanova(mypermanova = mypermanova,  
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") +
  theme(legend.position = "none")


p <- chibi.pco(list_ohpco = mpco,col_val = "typebyTissue",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_fraction() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure5_syncom_betadiversiy_pco1_pco2_sup_composition.pdf",outdir = "../figures/")

### Main Figure ###
Dat_rar <- Dat_rar %>% subset.Dataset(subset = typebyTissue != "Inoculum",
                                      drop = T,clean = T) 

#Permanova
Tab_bray <- distfun(t(Dat_rar$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition + typebyTissue:condition ,
                      data = Dat_rar$Map,strata = Dat_rar$Map$Rep,permutations = 10000)

capture.output(file = "../figures/figure5_syncom_betadiversity_summary.doc",
               append = T,print(mypermanova))

#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S5D_dist.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_rar$Map,file = "../data_figures/data_S5D_metadata.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4) 
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") +
  theme(legend.position = "none")


mpco <- oh.pco(Tab = t(Dat_rar$Tab),Map = Dat_rar$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "ID_Matrix")


p <- chibi.pco(list_ohpco = mpco,col_val = "typebyTissue",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,
               size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_fraction() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure5_syncom_betadiversiy_pco1_pco2_main_composition.pdf",outdir = "../figures/")

#Perform cap

mcap <- oh.cap(Tab = Dat_rar$Tab,Map = Dat_rar$Map,formula = "typebyTissue + condition + Condition(Rep)",
               distfun = distfun,perms = 10000,sqrt = T)



#Write the  dataset
write.table(x = mcap$Map_cap,file = "../data_figures/data_Fig5C.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



p <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
               mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
               size_legend_text = size_legend_text,
               size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure5_syncom_betadiversiy_cap1_cap2_main_composition.pdf",outdir = "../figures/")

### Subset Root
Dat_sub <- Dat_rar %>% subset.Dataset(subset = typebyTissue == "Root",drop = T,clean = T)

mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 10000,sqrt = T)


#Write the  dataset
write.table(x = mcap$Map_cap,file = "../data_figures/data_Fig5F.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)




#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 10000)
capture.output(file = "../figures/figure6_syncom_betadiversity_root_summary.doc",
               append = F,print(mypermanova))
p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4) 
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") +
  theme(legend.position = "none")


#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S6C_dist.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_sub$Map,file = "../data_figures/data_S6C_metadata.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)




p <- chibi.cap(list_ohpco = mcap,col_val = "condition",
               mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
               size_legend_text = size_legend_text,
               size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_manual(values = colores_phosphate) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure6_syncom_betadiversiy_cap1_cap2_root_composition.pdf",outdir = "../figures/")






### Subset Agar
Dat_sub <- Dat_rar %>% subset.Dataset(subset = typebyTissue == "AgarPlant",drop = T,clean = T)

mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 10000,sqrt = T)

#Write the  dataset
write.table(x = mcap$Map_cap,file = "../data_figures/data_S6B.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 10000)
capture.output(file = "../figures/figure6_syncom_betadiversity_agar_summary.doc",
               append = F,print(mypermanova))
p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4) 
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") +
  theme(legend.position = "none")


#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S6E_dist.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_sub$Map,file = "../data_figures/data_S6E_metadata.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)




p <- chibi.cap(list_ohpco = mcap,col_val = "condition",
               mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
               size_legend_text = size_legend_text,
               size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_manual(values = colores_phosphate) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure6_syncom_betadiversiy_cap1_cap2_agar_composition.pdf",outdir = "../figures/")



### Subset Shoot
Dat_sub <- Dat_rar %>% subset.Dataset(subset = typebyTissue == "Shoot",drop = T,clean = T)

mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 10000,sqrt = T)

#Write the  dataset
write.table(x = mcap$Map_cap,file = "../data_figures/data_S6A.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 10000)
capture.output(file = "../figures/figure6_syncom_betadiversity_shoot_summary.doc",
               append = F,print(mypermanova))
p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,size_axis_title = 15,
                          size_axis_text = 30,size_title_text = size_legend_text,legend_proportion_size = 4) 
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") +
  theme(legend.position = "none")

#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S6D_dist.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_sub$Map,file = "../data_figures/data_S6D_metadata.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)




#Root phosphate

p <- chibi.cap(list_ohpco = mcap,col_val = "condition",
               mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
               size_legend_text = size_legend_text,
               size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_manual(values = colores_phosphate) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure6_syncom_betadiversiy_cap1_cap2_shoot_composition.pdf",outdir = "../figures/")

