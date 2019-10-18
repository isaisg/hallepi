library(ohchibi)
library(palettesPM)
library(paletteer)

set.seed(130816)

size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4


Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_piDO_useq97.RDS")
Dat_rar <- Dat_ori$RelativeAbundance
Dat_rar <- Dat_rar %>% subset.Dataset(subset = SynCom == "Full",drop = T,clean = T)


colores_phosphate <- c("#ffffe5","#41ab5d","#005a32")
names(colores_phosphate) <- c("0uM","50uM","1000uM")

colores_geno <- c(palettesPM::pm.colors.fractions(),pm.colors.genotypes())
#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
mpco <- oh.pco(Tab = t(Dat_rar$Tab),Map = Dat_rar$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "ID_Matrix")


p_frac <- chibi.pco(list_ohpco = mpco,col_val = "Fraction",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_fraction() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p_pi <- chibi.pco(list_ohpco = mpco,col_val = "Pi",
                  mypch = 21,size = 20,alpha=1,
                  size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                  size_title_text = size_legend_title,
                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_manual(values = colores_phosphate) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


p_geno<- chibi.pco(list_ohpco = mpco,col_val = "Genotype",
                   mypch = 21,size = 20,alpha=1,
                   size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                   size_title_text = size_legend_title,
                   font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_manual(values = colores_geno) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
composition <- egg::ggarrange(p_frac,p_pi,p_geno, nrow = 1)


#### Check the genoype effect in the full experiment ###
Dat_sub <- Dat_rar %>% subset.Dataset(subset = Fraction == "Root",drop = T,clean = T)

mpco <- oh.pco(Tab = t(Dat_sub$Tab),Map = Dat_sub$Map,ndim = 3,eig = T,
               distfun = distfun,id_var = "ID_Matrix")


p_geno <- chibi.pco(list_ohpco = mpco,col_val = "Genotype",
                    mypch = 21,size = 20,alpha=1,
                    size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_manual(values = colores_geno) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


p_pi <- chibi.pco(list_ohpco = mpco,col_val = "Pi",
                  mypch = 21,size = 20,alpha=1,
                  size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                  size_title_text = size_legend_title,
                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T)  +
  scale_fill_manual(values = colores_phosphate) + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  Genotype + Pi + Genotype:Pi  ,
                      data = Dat_sub$Map,permutations = 5000)
mypermanova

### CAP constraining for the Pi effect
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,formula = "Genotype  + Condition(Pi)",
               distfun = distfun,perms = 5000,sqrt = T)

df_pv <- mypermanova$aov.tab %>% as.data.frame
dfpval <- data.frame(label = paste0("R2 = ",round(df_pv$R2[1],3),
                                    "\n p-value = ",base::format.pval(pv = df_pv$`Pr(>F)`[1],digits = 2))
)

#Write data fore figure S8
write.table(x = mcap$Map_cap,file = "../data_figures/data_S9AB.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



p_geno <- chibi.cap(list_ohpco = mcap,col_val = "Genotype",
                    mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                    size_legend_text = size_legend_text,
                    size_axis_title = size_axis_title,size_axis_line = 2,
                    size_title_text = size_legend_title,
                    font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_manual(values = colores_geno) +
  geom_text(dfpval,mapping = aes(1.3,-2.5,label =label),size = 10, color = "black",fontface = "bold") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,formula = "Pi  + Condition(Genotype)",
               distfun = distfun,perms = 5000,sqrt = T)


dfpval <- data.frame(label = paste0("R2 = ",round(df_pv$R2[2],3),
                                    "\n p-value = ",base::format.pval(pv = df_pv$`Pr(>F)`[2],digits = 2))
)

p_pi <- chibi.cap(list_ohpco = mcap,col_val = "Pi",
                  mypch = 21,size = 20,alpha=1,comp_a = "CAP1",comp_b = "CAP2",
                  size_legend_text = size_legend_text,
                  size_axis_title = size_axis_title,size_axis_line = 2,
                  size_title_text = size_legend_title,
                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_manual(values = colores_phosphate) + 
  geom_text(dfpval,mapping = aes(-0.75,1.8,label =label),size = 10, color = "black",fontface = "bold") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

composition <- egg::ggarrange(p_geno,p_pi,nrow = 1)
oh.save.pdf(p = composition,outname = "pido_genotype_pi_effect.pdf",
            outdir = "../figures/",width = 20,height = 10)

rm(list=ls())
dev.off()