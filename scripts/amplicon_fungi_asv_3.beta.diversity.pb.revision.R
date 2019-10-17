library(ohchibi)
library(palettesPM)
library(paletteer)
library(extrafont)
loadfonts(device = "pdf")
library(egg)
#Set random seed
set.seed(130816)


size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4

setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')

source('plotting_parameters_hallepi.R')

#Read bacterial asvs dat
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$Rarefied



#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")


###########################
###########################
######### All #############
###########################
###########################

#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T)
Dat_sub <- subset.Dataset(x = Dat_rar,Fraction!="BulkSoil",drop = T,clean = T)


mpco <- oh.pco(Tab = Dat_sub$Tab %>% t,Map = Dat_sub$Map,ndim = 4,
               eig = T,distfun = distfun,id_var = "DADA2_Header")
pco_frac <- chibi.pco(list_ohpco = mpco,col_val = "Fraction",
                      mypch = 21,size = 20,alpha=1,
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  scale_fill_fraction() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        
  )
pco_plot <- chibi.pco(list_ohpco = mpco,col_val = "Plot",
                      mypch = 21,size = 20,alpha=1,
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )


pco_phosphate<- chibi.pco(list_ohpco = mpco,col_val = "Phosphate",
                          mypch = 21,size = 20,alpha=1,
                          size_legend_text = size_legend_text,
                          size_axis_title = size_axis_title,size_axis_line = 2,
                          size_title_text = size_legend_title,
                          font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

pco_genotype<- chibi.pco(list_ohpco = mpco,col_val = "Genotype",
                         mypch = 21,size = 20,alpha=1,
                         size_legend_text = size_legend_text,
                         size_axis_title = size_axis_title,size_axis_line = 2,
                         size_title_text = size_legend_title,
                         font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

egg::ggarrange(pco_frac,pco_plot,
               pco_phosphate,pco_genotype,nrow = 2,ncol = 2)


#Control for the plot effect observed
mycap <-oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "Fraction + Genotype + Phosphate  + Condition(Plot) ",
               distfun = distfun,perms = 10000,sqrt = T)


mycap$perm_anova_terms
capture.output(file = "../figures/figure2_fungi_asvs_samples_soilrootshoot_summary.doc",
               append = F,print(mycap$perm_anova_terms))

#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  Fraction + Genotype  + Phosphate + Fraction:Genotype +
                        Fraction:Phosphate + Genotype:Phosphate + Fraction:Genotype:Phosphate + Plot,
                      data = Dat_sub$Map,permutations = 10000)

mypermanova
capture.output(file = "../figures/figure2_fungi_asvs_samples_soilrootshoot_summary.doc",
               append = T,print(mypermanova))

#Build permanova
p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,
                          size_axis_title = 35,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")
oh.save.pdf(p = p_perm,outname = "figure2_fungi_asvs_beta_diversity_maineffects_chibipermanova.pdf",outdir = "../figures/")

#Plot cap 
p <- chibi.cap(list_ohpco = mycap,col_val = "Fraction",comp_a = "CAP1",comp_b = "CAP2",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_fraction() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
oh.save.pdf(p = p,outname = "figure2_fungi_asvs_beta_diversity_maineffects_cap1_cap2.pdf",outdir = "../figures/")

p <- p+ theme(legend.position = "none")
p_perm <- p_perm + theme(legend.position = "none")


composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))

oh.save.pdf(p = composition,outname = "figure2_fungi_asvs_beta_diversity_maineffects_composition.pdf",outdir = "../figures/")

#Check the identity of the soil outliers
mycap$Map_cap %>% subset(Fraction == "Soil" & CAP1 < 0) 

###########################
###########################
######### Root #############
###########################
###########################
#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root",drop = T, clean =T)

mpco <- oh.pco(Tab = Dat_sub$Tab %>% t,Map = Dat_sub$Map,ndim = 4,
               eig = T,distfun = distfun,id_var = "DADA2_Header")
pco_frac <- chibi.pco(list_ohpco = mpco,col_val = "Fraction",
                      mypch = 21,size = 20,alpha=1,
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  scale_fill_fraction() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        
  )
pco_plot <- chibi.pco(list_ohpco = mpco,col_val = "Plot",
                      mypch = 21,size = 20,alpha=1,
                      size_legend_text = size_legend_text,
                      size_axis_title = size_axis_title,size_axis_line = 2,
                      size_title_text = size_legend_title,
                      font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

#egg::ggarrange(pco_frac,pco_plot,nrow = 1)

pco_phosphate<- chibi.pco(list_ohpco = mpco,col_val = "Phosphate",
                          mypch = 21,size = 20,alpha=1,
                          size_legend_text = size_legend_text,
                          size_axis_title = size_axis_title,size_axis_line = 2,
                          size_title_text = size_legend_title,
                          font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

pco_genotype<- chibi.pco(list_ohpco = mpco,col_val = "Genotype",
                         mypch = 21,size = 20,alpha=1,
                         size_legend_text = size_legend_text,
                         size_axis_title = size_axis_title,size_axis_line = 2,
                         size_title_text = size_legend_title,
                         font_family = "Arial",legend_proportion_size = 0.5,lines_zero = T)  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

egg::ggarrange(pco_frac,pco_plot,
               pco_phosphate,pco_genotype,nrow = 2,ncol = 2)



#Control for the plot effect
mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Phosphate  + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)


mycap$perm_anova_terms
capture.output(file = "../figures/figure3_fungi_asvs_samples_root_summary.doc",
               append = F,print(mycap$perm_anova_terms))


#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/figure3_fungi_asvs_samples_root_summary.doc",
               append = T,print(mypermanova))


p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,
                          size_axis_title = 35,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")

#Plot cap 
p <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
               shape_val = "Genotype",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_color_manual(values = palette_pi_soil) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dfpval <- data.frame (label = paste0("R2=",round((p_perm$data[2,1])/100,3),
                                     " p-value=",format.pval(p_perm$data[2,2])))



p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  ) +
  geom_text(dfpval,mapping = aes(1.1,1.5,label =label),size = 10, color = "black",fontface = "bold") 

dfpval <- data.frame (label = paste0("R2=",round((p_perm$data[1,1])/100,3),
                                     " p-value=",format.pval(p_perm$data[1,2])))


p_phos <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_manual(values = palette_pi_soil) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        
  )+
  geom_text(dfpval,mapping = aes(1.1,1.5,label =label),size = 10, color = "black",fontface = "bold") 


composition <- egg::ggarrange(p_geno,p_phos,nrow = 1)

oh.save.pdf(composition,outname = "figure_3_legendcaps_twopanels_root_fungi.pdf",
            outdir = "../figures",width = 20,height = 10)



p <- p+ theme(legend.position = "none")
p_perm <- p_perm + theme(legend.position = "none")
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure3_fungi_asvs_beta_diversity_onlyroot_cap1_cap2_composition.pdf",outdir = "../figures/")


###Test the genotype effect inside each phosphate level ####
#low#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root" & Phosphate == "low",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_root_low_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_root_low_genotypeeffect.doc",
               append = T,print(mypermanova))


#medium#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root" & Phosphate == "medium",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_root_medium_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_root_medium_genotypeeffect.doc",
               append = T,print(mypermanova))

#high#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root" & Phosphate == "high",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_root_high_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_root_high_genotypeeffect.doc",
               append = T,print(mypermanova))

#low+Pi#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root" & Phosphate == "low+Pi",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=0.9,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_root_low_Pi_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)


Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_root_low_Pi_genotypeeffect.doc",
               append = T,print(mypermanova))

###########################
###########################
######### Soil #############
###########################
###########################
#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Soil",drop = T, clean =T)

#Control for the plot effect
mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Phosphate  + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)


mycap$perm_anova_terms
capture.output(file = "../figures/figure3_fungi_asvs_samples_soil_summary.doc",
               append = F,print(mycap$perm_anova_terms))


Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/figure3_fungi_asvs_samples_soil_summary.doc",
               append = T,print(mypermanova))

p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,
                          size_axis_title = 35,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")

#Plot cap 
p <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
               shape_val = "Genotype",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_color_manual(values = palette_pi_soil) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df_temp <- mypermanova$aov.tab %>% as.data.frame
dfpval <- data.frame (label = paste0("R2=",round((df_temp[1,5]),3),
                                     " p-value=",format.pval(df_temp[1,6])))


p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )+
  geom_text(dfpval,mapping = aes(-2,2,label =label),size = 10, color = "black",fontface = "bold") 

dfpval <- data.frame (label = paste0("R2=",round((df_temp[2,5]),3),
                                     " p-value=",format.pval(df_temp[2,6])))


p_phos <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_manual(values = palette_pi_soil) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        
  )+
  geom_text(dfpval,mapping = aes(-2,2,label =label),size = 10, color = "black",fontface = "bold") 



composition <- egg::ggarrange(p_geno,p_phos,nrow = 1)

oh.save.pdf(composition,outname = "figure_3_legendcaps_twopanels_soil_fungi.pdf",outdir = "../figures",
            width = 20,height = 10)



#oh.save.pdf(p = p,outname = "figure3_bacteria_asvs_beta_diversity_onlyroot_cap1_cap2.pdf",outdir = "../figures/")
p <- p+ theme(legend.position = "none")
p_perm <- p_perm + theme(legend.position = "none")
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure3_fungi_asvs_beta_diversity_onlysoil_cap1_cap2_composition.pdf",outdir = "../figures/")


###########################
###########################
######### Shoot #############
###########################
###########################
#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot",drop = T, clean =T)

#Control for the plot effect
mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Phosphate  + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)


mycap$perm_anova_terms
capture.output(file = "../figures/figure3_fungi_asvs_samples_shoot_summary.doc",
               append = F,print(mycap$perm_anova_terms))


Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/figure3_fungi_asvs_samples_shoot_summary.doc",
               append = T,print(mypermanova))

p_perm <- chibi.permanova(mypermanova = mypermanova,
                          size_legend_text = size_legend_text,
                          size_axis_title = 35,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")

#Plot cap 
p <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
               shape_val = "Genotype",
               mypch = 21,size = 20,alpha=1,
               size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
               size_title_text = size_legend_title,
               font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_color_manual(values = palette_pi_soil) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

df_temp <- mypermanova$aov.tab %>% as.data.frame
dfpval <- data.frame (label = paste0("R2=",round((df_temp[1,5]),3),
                                     " p-value=",format.pval(df_temp[1,6])))


p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )+
  geom_text(dfpval,mapping = aes(0,2,label =label),size = 10, color = "black",fontface = "bold") 



dfpval <- data.frame (label = paste0("R2=",round((df_temp[2,5]),3),
                                     " p-value=",format.pval(df_temp[2,6])))

p_phos <- chibi.cap(list_ohpco = mycap,col_val = "Phosphate",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_manual(values = palette_pi_soil) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        
  )+
  geom_text(dfpval,mapping = aes(0,2,label =label),size = 10, color = "black",fontface = "bold") 



composition <- egg::ggarrange(p_geno,p_phos,nrow = 1)

oh.save.pdf(composition,
            outname = "figure_3_legendcaps_twopanels_shoot_fungi.pdf",
            outdir = "../figures",width = 20,height = 10)



#oh.save.pdf(p = p,outname = "figure3_bacteria_asvs_beta_diversity_onlyroot_cap1_cap2.pdf",outdir = "../figures/")
p <- p+ theme(legend.position = "none")
p_perm <- p_perm + theme(legend.position = "none")
composition <- egg::ggarrange(p_perm,p, widths = c(0.05,1))
oh.save.pdf(p = composition,outname = "figure3_fungi_asvs_beta_diversity_onlyshoot_cap1_cap2_composition.pdf",
            outdir = "../figures/")

#low#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot" & Phosphate == "low",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=0.9,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_shoot_low_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_shoot_low_genotypeeffect.doc",
               append = T,print(mypermanova))


#medium#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot" & Phosphate == "medium",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=0.9,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_shoot_medium_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_shoot_medium_genotypeeffect.doc",
               append = T,print(mypermanova))

#high#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot" & Phosphate == "high",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=0.9,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_shoot_high_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)

Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_shoot_high_genotypeeffect.doc",
               append = T,print(mypermanova))

#low+Pi#
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot" & Phosphate == "low+Pi",drop = T, clean =T)

mycap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
                formula = "Genotype + Condition(Plot) ",
                distfun = distfun,perms = 10000,sqrt = T)

p_geno <- chibi.cap(list_ohpco = mycap,col_val = "Genotype",comp_a = "CAP1",comp_b = "CAP2",
                    mypch = 21,size = 20,alpha=0.9,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) + 
  scale_fill_genotype() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
  )

oh.save.pdf(p_geno,outname = "supfigure3_fungi_permanova_shoot_low_Pi_genotypeeffect.pdf",
            outdir = "../figures",width = 20,height = 20)


Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova <- adonis(Tab_bray ~   Genotype + Plot,
                      data = Dat_sub$Map,permutations = 10000)
mypermanova
capture.output(file = "../figures/supfigure3_fungi_permanova_shoot_low_Pi_genotypeeffect.doc",
               append = T,print(mypermanova))


