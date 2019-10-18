library(ohchibi)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf")
library(egg)

#Set random seed
set.seed(130816)



size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4


source('plotting_parameters_hallepi.R')

#Read bacterial otus dat
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
Dat_rar <- Dat$Rarefied



#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")


#Load the edaphic factors dataset
dat_edaphic <- read.table(file = "../rawdata/dat_edaphic_factors_hallepi.tsv",
                          header = T,sep = "\t")
dat_edaphic$Plot <- dat_edaphic$Plot %>% paste0("Plot",.)

#Remove P Silt and Clay
dat_edaphic <- which(!(dat_edaphic$Parameter %in% c("P.a","P.t","Silt","Clay"))) %>%
  dat_edaphic[.,] %>% droplevels

#Create a matrix and explore the correlated nature of the edaphic factor
Tab_ed <- acast(data = dat_edaphic,Plot ~ Parameter,value.var = "Value") %>%
  scale

#Tab_ed

#Write dataset
write.table(x = Tab_ed,file = "../data_figures/data_S10A.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



#Create correlation plot
mggcor <- chibi.ggcor(Tab = Tab_ed,size_r_display = 5)

#Due to the existence of the correlation structure analyze deal with it using pca
mpca <- Tab_ed %>% prcomp(x = .,center = F,scale. = F,retx = T)
msum <- summary(mpca) 
mdf <- msum$importance[2,] %>% data.frame(PC = names(.),Var = .,row.names = NULL)
mdf$Var <- mdf$Var*100
mdf$PC <- factor(mdf$PC,levels = mdf$PC %>% as.character)
p_pca <- ggplot(data = mdf,aes(PC,Var)) +
  geom_bar(stat = "identity") + 
  theme_ohchibi() +
  theme(
    axis.text.x = element_text(angle = 90,family = "Arial",vjust = 0.5,hjust = 1)
  ) +
  xlab(label = NULL) + ylab(label = "Variance Explained (%)") +
  scale_y_continuous(breaks = seq(0,100,5))

#Write dataset
write.table(x = mdf,file = "../data_figures/data_S10B.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



#Now see the contribution of each one to the edaphic factors
melted <- mpca$rotation %>% melt

#Subset onlyu the three first components
chosen <- c("PC1","PC2","PC3")
melted <- which(melted$Var2 %in% chosen) %>%
  melted[.,] %>% droplevels
melted$a_value <- abs(melted$value)
melted$Direction <- "Down"
melted$Direction[which(melted$value > 0)] <- "Up"

p_cont <- ggplot(data = melted,aes(Var1,a_value)) +
  geom_bar(stat = "identity",aes(fill = Direction)) +
  scale_fill_manual(values = c("blue","red")) + 
  facet_grid(.~Var2,space = "free",scales = "free") +
  theme_ohchibi() +
  theme(
    axis.text.x = element_text(angle = 90,family = "Arial",vjust = 0.5,hjust = 1),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold")
  ) +
  xlab(label = NULL) + ylab(label = "Contribution PC") +
  scale_y_continuous(breaks = seq(0,1,0.1))

composition <- egg::ggarrange(mggcor$p,p_pca,p_cont,nrow = 1)
oh.save.pdf(p = composition,outname = "edaphic_composition_red.pdf",outdir = "../figures/",
            width = 30,height = 10)

##Write the table

write.table(x = melted,file = "../data_figures/data_S10C.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



#Merge with the dat object
df_ed <- mpca$x %>% as.data.frame
df_ed$Plot <- rownames(df_ed)
Dat_rar$Map <- merge(Dat_rar$Map,df_ed, by = "Plot",all.x = TRUE) 
rownames(Dat_rar$Map) <- Dat_rar$Map$DADA2_Header
Dat_rar$Map <- match(colnames(Dat_rar$Tab),rownames(Dat_rar$Map)) %>%
  Dat_rar$Map[.,]




#Build permanova
paleta <- palette_variance
paleta_pc <- c("blue","purple","green")
names(paleta_pc) <- c("PC1","PC2","PC3")
palette_variance <- c(paleta,paleta_pc)


###########################
###########################
######### Root #############
###########################
###########################
#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root",drop = T, clean =T)



#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova_o <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + Plot,
                        data = Dat_sub$Map,permutations = 10000)
mypermanova_o
capture.output(file = "../figures/edaphic_summary_bacteria_summary.doc",
               append = T,print(mypermanova_o))


mypermanova_e <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + PC1 + PC2 + PC3 + Plot,
                        data = Dat_sub$Map,permutations = 10000)
mypermanova_e
capture.output(file = "../figures/edaphic_summary_bacteria_summary.doc",
               append = T,print(mypermanova_e))

#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S10D_dist_bacteria.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_sub$Map,file = "../data_figures/data_S10D_metadata_bacteria.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



p_perm <- chibi.permanova(mypermanova = mypermanova_o,
                          size_legend_text = size_legend_text,
                          size_axis_title = 15,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot",
)
p_perm_o <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model") +
  theme(
    legend.position = "none"
  )


p_perm <- chibi.permanova(mypermanova = mypermanova_e,
                          size_legend_text = size_legend_text,
                          size_axis_title = 15,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm_e <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")


composition <- egg::ggarrange(p_perm_o,p_perm_e,nrow = 1)
oh.save.pdf(p = composition,outname = "edaphic_permanova_root_bacteria.pdf",outdir = "../figures/")




### Repeat for the fungi dataset
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$Rarefied

Dat_rar$Map <- merge(Dat_rar$Map,df_ed, by = "Plot",all.x = TRUE) 
rownames(Dat_rar$Map) <- Dat_rar$Map$DADA2_Header
Dat_rar$Map <- match(colnames(Dat_rar$Tab),rownames(Dat_rar$Map)) %>%
  Dat_rar$Map[.,]



###########################
###########################
######### Root #############
###########################
###########################
#Remove bulksoil samples
Dat_sub <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root",drop = T, clean =T)



#Permanova
Tab_bray <- distfun(t(Dat_sub$Tab))

mypermanova_o <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + Plot,
                        data = Dat_sub$Map,permutations = 10000)
mypermanova_o
capture.output(file = "../figures/edaphic_summary_fungi_summary.doc",
               append = T,print(mypermanova_o))


mypermanova_e <- adonis(Tab_bray ~   Genotype  + Phosphate  + Genotype:Phosphate + PC1 + PC2 + PC3 + Plot,
                        data = Dat_sub$Map,permutations = 10000)
mypermanova_e
capture.output(file = "../figures/edaphic_summary_fungi_summary.doc",
               append = T,print(mypermanova_e))


p_perm <- chibi.permanova(mypermanova = mypermanova_o,
                          size_legend_text = size_legend_text,
                          size_axis_title = 15,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot",
)
p_perm_o <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model") +
  theme(
    legend.position = "none"
  )


p_perm <- chibi.permanova(mypermanova = mypermanova_e,
                          size_legend_text = size_legend_text,
                          size_axis_title = 15,
                          size_axis_text = 30,
                          size_title_text = size_legend_text,
                          legend_proportion_size = 4,
                          terms_exclude_plot = "Plot"
)
p_perm_e <- p_perm$p + 
  scale_fill_manual(values = palette_variance) +
  xlab(label = "Term Model")


composition <- egg::ggarrange(p_perm_o,p_perm_e,nrow = 1)
oh.save.pdf(p = composition,outname = "edaphic_permanova_fungi_bacteria.pdf",outdir = "../figures/")


#Print data used to compute permanova
bray_tab <- Tab_bray %>% as.matrix %>%as.data.frame
write.table(x =bray_tab,file = "../data_figures/data_S10D_dist_fungi.csv",
            append = F,quote = F,sep = ",",row.names = T,col.names = T)


write.table(x =Dat_sub$Map,file = "../data_figures/data_S10D_metadata_fungi.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)

