library(ohchibi)
library(paletteer)
library(adephylo)
library(phytools)
library(egg)
library(extrafont)
loadfonts(device = "pdf") 



#Read enricment results
res_frac <- read.table(file = "../cleandata/df_dds_res_syncom_fraction_specific.tsv",header = T,sep = "\t")
res_int <- read.table(file = "../cleandata/df_dds_res_syncom_phosphateint.tsv",header = T,sep = "\t")

res <- rbind(res_frac,res_int)
res$Facet <- res$Facet %>% 
  factor(levels = c("Root vs Agar","Shoot vs Agar","Root vs Shoot","Phosphate Interaction"))
res$Significance <- rep("NotSignificant",nrow(res))
res$Significance[which(res$padj < 0.1)] <- "Significant"
res$Significance <- res$Significance %>% factor
res$log2FoldChange[which(res$log2FoldChange>5)] <- 5
res$log2FoldChange[which(res$log2FoldChange< -5)] <- -5

write.table(x = res,file = "../data_figures/data_Fig5D.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)



#Read phylogenetic tree
tree <- read.tree(file = "../cleandata/201_uncisolates_47markers.rooted.newick")
meta <- read.table(file = "../rawdata/metadata_97useq_unc_genomes_useq_freezer_mothur_taxonoid_genome_name.tsv",
                   header = T,sep = "\t",quote = "",comment.char = "")

suptab <- meta
suptab <- suptab[-3]
suptab$Useq <- suptab$Useq %>% gsub(pattern = "Sequence_",replacement = "USeq")
suptab$Num <- suptab$Useq %>% gsub(pattern = "USeq",replacement = "") %>% as.numeric
suptab <- with(suptab,order(Num)) %>% suptab[.,]
suptab <- suptab[,-6]

write.table(x = suptab,file = "../cleandata/sup_table_8.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)

ids <- res$Id %>% as.character %>% unique
meta <- which(meta$Useq %in% ids) %>% meta[.,] %>%
  droplevels
df_useq2taxon <- NULL
for(useq in levels(meta$Useq)){
  tids <- meta %>% subset(Useq == useq) %$% taxon_oid %>%
    na.omit %>% as.character %>% unique
  if(length(tids) ==0){
  }else if(length(tids) ==1){
    df_useq2taxon <- data.frame(USeq = useq,Representative = tids) %>%
      rbind(df_useq2taxon,.)
  }else{
    dfdist <- distRoot(tree,tids) %>% 
      data.frame(taxon_oid = names(.),Distance = .,row.names = NULL)
    dfdist <- with(dfdist,order(Distance)) %>% dfdist[.,]
    rep <- dfdist[1,1] %>% as.character
    df_useq2taxon <- data.frame(USeq = useq,Representative = rep) %>%
      rbind(df_useq2taxon,.)
  }
}
df_useq2taxon$Representative <- df_useq2taxon$Representative %>% as.character
df_useq2taxon$Representative[which(df_useq2taxon$USeq == "Sequence_86")] <- "2639762526"

#Subset the tree
todrop <- which(!(tree$tip.label %in% df_useq2taxon$Representative)) %>%  
  tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)
tree$tip.label <- match(tree$tip.label,df_useq2taxon$Representative) %>% 
  df_useq2taxon$USeq[.] %>% as.character
#ggtree(tree) + geom_text2(aes(subset=!isTip, label=node)) + geom_tiplab()
#tree <- ape::rotate(phy = tree,node = 108)
#tree <- ape::rotate(phy = tree,node = 111)


#Determine the orders of the tips 
tree <- tree %>% ladderize
#tree <- ape::rotate(phy = tree,node = 149)
#tree <- ape::rotate(phy = tree,node = 144)
#tree <- ape::rotate(phy = tree,node = 145)

is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]


#Customize the tree 
mtips <- tree$tip.label[ordered_tips] %>% rev


#Change some labels 
taxa_name <- match(df_useq2taxon$Representative,meta$taxon_oid) %>% 
  meta$Genome_Classification[.] %>% as.character
df_tree_tax <- taxa_name %>% strsplit(split = "\\;") %>% unlist %>%
  gsub(pattern = "[k|p|c|o|f|g]__",replacement = "") %>%
  matrix(data = .,ncol = 6,byrow = T) %>% as.data.frame
colnames(df_tree_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus")
df_tree_tax$Representative  <- df_useq2taxon$Representative
df_useq2taxon <- merge(df_useq2taxon,df_tree_tax, by = "Representative")
# 
# df_useq2taxon <- match(mtips,df_useq2taxon$USeq) %>% 
#   df_useq2taxon[.,]
df_useq2taxon <- match(tree$tip.label, df_useq2taxon$USeq) %>% 
  df_useq2taxon[.,]

mphyla <- palettesPM::pm.colors.phyla()
df_useq2taxon$PhylaColor <- match(df_useq2taxon$Phylum,names(mphyla)) %>%
  mphyla[.] %>% as.character

df_useq2taxon$Genus <- df_useq2taxon$Genus %>% as.character
#Sequence_36 is a leucobacter
df_useq2taxon$Genus[which(df_useq2taxon$USeq == "Sequence_36")] <- "Leucobacter"

#Create another set of labels to visualize
df_useq2taxon$NLabels <- df_useq2taxon$USeq %>% as.character %>% 
  gsub(pattern = "Sequence_",replacement = "USeq") %>%
  paste(.,df_useq2taxon$Genus,sep = " ")

tree_plot <- tree
tree_plot$tip.label <- df_useq2taxon$NLabels

plot(tree)
nodelabels()

cairo_pdf(filename = "../figures/figure5_tree_heatmap.pdf",onefile = FALSE, fallback_resolution = 1200,
          width = 5, height = 15, family = "Arial", antialias = "default",
          pointsize = 12)
plot(tree_plot,align.tip.label = T,show.tip.label = T,adj = 1,
     no.margin = TRUE,cex = 1)
tiplabels(pch = 21,col = "#414141",
          bg = df_useq2taxon$PhylaColor,cex = 1)
dev.off()

write.table(x = df_useq2taxon,file = "../cleandata/df_useq2taxon.tsv",append = F,
            quote = F,sep = "\t",row.names = F,col.names = T)

###### Heatmap #####
res <- which(res$Id %in% tree$tip.label) %>% res[.,] %>% droplevels
res$Id <- res$Id %>% factor(levels =(tree$tip.label[ordered_tips]))
#Reorder levels of Contrasts 
res$Contrast <- res$Contrast %>% factor(levels = 
                                          c("Root0Pi_vs_Agar0Pi","Root10Pi_vs_Agar10Pi","Root30Pi_vs_Agar30Pi","Root50Pi_vs_Agar50Pi",
                                            "Root100Pi_vs_Agar100Pi","Root1000Pi_vs_Agar1000Pi","Shoot0Pi_vs_Agar0Pi","Shoot10Pi_vs_Agar10Pi",
                                            "Shoot30Pi_vs_Agar30Pi","Shoot50Pi_vs_Agar50Pi","Shoot100Pi_vs_Agar100Pi","Shoot1000Pi_vs_Agar1000Pi",
                                            "Root0Pi_vs_Shoot0Pi","Root10Pi_vs_Shoot10Pi","Root30Pi_vs_Shoot30Pi","Root50Pi_vs_Shoot50Pi",
                                            "Root100Pi_vs_Shoot100Pi","Root1000Pi_vs_Shoot1000Pi","TissueRoot.PhosphateIntLow","TissueShoot.PhosphateIntLow"))

#Putative Palttes
#kovesi.diverging_bwr_55_98_c37
#kovesi.diverging_gwv_55_95_c39
#kovesi.diverging_gwr_55_95_c38

p_heatmap <- ggplot(data = res,mapping = aes(x = Contrast, y = Id,fill = log2FoldChangePlot)) + 
  geom_raster(aes(fill = log2FoldChange)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.85,height = 0.85) + 
  facet_grid(~Facet,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),na.value = "white") +
  scale_color_manual(values = c("#00000000","#414141"))+
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 10,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.x = 0,
                legend_proportion_size = 0,size_title_text = 0,size_legend_text = 0,
                size_panel_border = 0.1,size_lines_panel = 0) +
  theme(aspect.ratio = 0.5,axis.ticks = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "none")

p_heatmap_legend <- ggplot(data = res,mapping = aes(x = Contrast, y = Id,fill = log2FoldChangePlot)) + 
  geom_raster(aes(fill = log2FoldChange)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.85,height = 0.85) + 
  facet_grid(~Facet,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),na.value = "white") +
  scale_color_manual(values = c("#00000000","#414141")) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 10,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.x = 0,
                size_panel_border = 0.1,size_lines_panel = 0)
oh.save.pdf(p = p_heatmap_legend,outname = "figure_5_heatmap_legends.pdf",outdir = "../figures/")

#Comput ephylogenetic signal
Tab <- acast(data = res,formula =Id ~ Contrast ,value.var = "log2FoldChange") 
df_ps <- apply(X = Tab,MARGIN = 2,
               FUN = function(x)phylosig(tree = tree,x,method = "lambda",test = TRUE) %>% unlist) %>%
  melt
colnames(df_ps) <- c("Type","Contrast","value")
df_ps_p <- df_ps %>% subset(Type == "P") 
df_ps_p$value <- df_ps_p$value %>% p.adjust(method = "fdr")
df_ps_p$Type <- df_ps_p$Type %>% as.character %>% 
  gsub(pattern = "P",replacement = "padj")
df_ps <- rbind(df_ps,df_ps_p )

df_ps_p$Facet <- df_ps_p$Contrast %>% as.character
df_ps_p$Facet <- df_ps_p$Facet %>% gsub(pattern = "Root.*Agar.*",replacement = "Root vs Agar") %>%
  gsub(pattern = "Shoot.*Agar.*",replacement = "Shoot vs Agar") %>%
  gsub(pattern = "Root.*Shoot.*",replacement = "Root vs Shoot") %>%
  gsub(pattern = ".*Int.*",replacement = "Phosphate Interaction") %>%
  factor(levels = c("Root vs Agar","Shoot vs Agar","Root vs Shoot", "Phosphate Interaction"))


p_signal <- df_ps_p %>% ggplot(data = .,mapping = aes(x = Contrast, y = -log10(value))) + 
  geom_bar(stat = "identity") + facet_grid(~Facet,scales = "free",space = "free") +
  ylab(label = "-log10(padj)") +
  theme_ohchibi(size_axis_text.x = 0,size_axis_title.y =10,size_axis_title.x = 0,
                size_strip_text.x = 0,size_axis_text.y = 6,size_panel_border = 0.1) +
  theme(aspect.ratio = 0.5,axis.ticks.x = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "none") + scale_y_reverse() +
  geom_hline(yintercept = 1.30103,size = 0.3, color = "firebrick")



## Abundance
Dat <- readRDS(file = "../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat <- Dat$RelativeAbundance
ra_melted <- Dat$Tab %>% melt
colnames(ra_melted) <- c("USeq","ID_Matrix","value")
ra_melted <- which(ra_melted$USeq %in% unique(res$Id)) %>% ra_melted[.,] %>%
  droplevels
ra_melted <- ra_melted %>% merge(Dat$Map,by  = "ID_Matrix")
ra_melted <- dcast(data = ra_melted,formula = USeq ~ typebyTissue,value.var = "value",fun.aggregate = mean) %>%
  melt
ra_melted <- ra_melted %>% subset(variable == "Agar_SC" | variable == "Root_SC" |
                                    variable == "Shoot_SC") %>%
  droplevels

ra_melted$USeq <- ra_melted$USeq %>% factor(levels = res$Id %>% levels)


p_abundance <- ggplot(data = ra_melted,mapping = aes(x= USeq, y = value)) + 
  geom_bar(stat = "identity") +facet_grid(~variable) + 
  coord_flip() + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                               size_axis_title.y = 0,size_axis_text.x = 2.5,
                               size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 270),
        panel.spacing = unit(0.25, "lines"),panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",aspect.ratio = 0.5,
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9")) +
  scale_y_sqrt(breaks = c(0.01,0.05,0.1,0.2,0.3))  +
  geom_hline(yintercept = 0.01,color = "firebrick",size = 0.3)

#Create a heder with the descriptions of each column
colores_header <- c(RColorBrewer::brewer.pal(n=6,"Oranges"),"white","white")
df_header <- data.frame(Contrast = res$Contrast %>% levels,
                        Type = rep("Header",20),
                        Condition = res$Contrast %>% levels %>% gsub(pattern = "Root|Shoot|_vs_.*",replacement = ""))
df_header$Facet <- c(rep("RootvsAgar",6),rep("ShootvsAgar",6),rep("RootvsShoot",6),rep("Interaction",2)) %>%
  factor(levels = c("RootvsAgar","ShootvsAgar","RootvsShoot","Interaction"))
p_header <- ggplot(data = df_header,mapping = aes(Contrast, y = Type, fill = Condition)) +
  geom_tile() + facet_grid(~Facet,scales = "free",space = "free") + scale_fill_manual(values = colores_header) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,
                size_axis_title.x = 0,size_axis_title.y = 0,legend_proportion_size = 0,
                size_title_text = 0,size_legend_text = 0,
                size_strip_text.x = 0,size_panel_border = 0.1,size_lines_panel = 0)+
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        aspect.ratio = 2.5,legend.position = "none")

##Count the number of sequences per useq
meta$Count <- rep(1,nrow(meta))
df_numtaxons <- acast(data = meta,formula =  taxon_oid ~ Useq,value.var = "Count",fun.aggregate = sum) %>%
  melt %>% subset(value != 0) %$% Var2 %>% table %>% as.data.frame
colnames(df_numtaxons)[1] <- "USeq"
df_numtaxons <- match(levels(res$Id),df_numtaxons$USeq) %>% df_numtaxons[.,] %>% droplevels
df_numtaxons$USeq <- df_numtaxons$USeq %>% factor(levels = levels(res$Id))


p_num <- ggplot(df_numtaxons,aes(USeq,Freq)) + geom_bar(stat = "identity") +
  geom_bar(stat = "identity") +
  coord_flip() + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                               size_axis_title.y = 0,size_axis_text.x = 2.5,
                               size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 270),
        panel.spacing = unit(0.25, "lines"),
        legend.position = "none",aspect.ratio = 0.5,
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  scale_y_continuous(breaks = c(0,1,3,7,13))

#Create a bar for genera based on black and white
#Genus
df_numtaxons <- merge(df_numtaxons,df_useq2taxon[,c(2:8)], by = "USeq")
df_numtaxons <- match(levels(df_numtaxons$USeq),df_numtaxons$USeq) %>% df_numtaxons[.,]
uninames <- df_numtaxons$Genus %>% unique
df_barcolors <- data.frame(Genus = uninames, Color = c(rep(c("white","black"),16),"white"))
df_numtaxons_genus <- merge(df_numtaxons,df_barcolors, by = "Genus")

paleta_blanco_negro <- c("white","black")
names(paleta_blanco_negro) <- c("white","black")
#Now create a ggplot structure to plot the bar graph
p_colorbar <- ggplot(df_numtaxons_genus,aes(Kingdom,USeq, fill = Color,color = Color)) +
  geom_tile(size = 0)  + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                                       size_axis_title.y = 0,size_axis_text.x = 0,
                                       size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),axis.text.y = element_blank(),
        legend.position = "none",aspect.ratio = 0.5,
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_cartesian(expand = FALSE) + scale_fill_manual(values = paleta_blanco_negro) +
  scale_color_manual(values = paleta_blanco_negro)

p_colorbar_genus <- p_colorbar

#Family
uninames <- df_numtaxons$Family %>% unique
df_barcolors <- data.frame(Family = uninames, Color = c(rep(c("white","black"),12),"white"))
df_numtaxons_family <- merge(df_numtaxons,df_barcolors, by = "Family")

p_colorbar <- ggplot(df_numtaxons_family,aes(Kingdom,USeq, fill = Color,color = Color)) +
  geom_tile(size = 0)  + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                                       size_axis_title.y = 0,size_axis_text.x = 0,
                                       size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),axis.text.y = element_blank(),
        legend.position = "none",aspect.ratio = 0.5,
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_cartesian(expand = FALSE) + scale_fill_manual(values = paleta_blanco_negro) +
  scale_color_manual(values = paleta_blanco_negro)

p_colorbar_family <- p_colorbar



p_blank <- ggplot() + theme_minimal()
# composition <- egg::ggarrange(p_header,p_blank,p_blank,p_heatmap,p_abundance,p_blank,p_signal,
#                               ncol = 3,nrow = 3,
#                heights = c(0.1,1,0.1),padding = 0,
#                debug = FALSE,widths = c(1,0.3,0),draw = FALSE)

composition <- egg::ggarrange(p_header,p_blank,p_blank,
                              p_heatmap,p_abundance,p_colorbar_genus,
                              p_signal,p_blank,p_blank,
                              ncol = 3,nrow = 3,
                              heights = c(0.1,1,0.1),padding = 0,
                              debug = FALSE,widths = c(1,0.5,0.05),draw = FALSE)

oh.save.pdf(p = composition,outname = "figure5_syncom_heatmap_fractioneffect.pdf",
            height = 15,width = 15,outdir = "../figures/")


