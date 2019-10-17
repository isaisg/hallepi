library(ohchibi)

setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')
set.seed(130816)
pval_thres <- 0.1
res <- read.table(file = "../cleandata/df_dds_res_amplicon_fungi_asvs_fraction_asvlevel.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")

#Read phylum information
#Read palettes at the order level
paleta <- readRDS(file = "../rawdata/hallepi_palette_fungi_orders.RDS")
morder <- paleta %>% names
res$Order <- res$Order %>% as.character %>% gsub(pattern = " o__",replacement = "")
res$Order[which(!(res$Order %in% morder))] <- "Other"
res$Order <- res$Order %>% factor(levels = morder)
df_all <- NULL

#Create structures for enrichments
df_freq <- res %>% subset(Contrast == "Root-Soil" & 
                            padj < pval_thres & log2FoldChange >0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Enriched_RootvsSoil")
df_all <- rbind(df_all,df_freq)

#Create structures for enrichments
df_freq <- res %>% subset(Contrast == "Root-Soil" & 
                            padj < pval_thres & log2FoldChange <0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Depleted_RootvsSoil")
df_all <- rbind(df_all,df_freq)

df_freq <- res %>% subset(Contrast == "Shoot-Soil" & 
                            padj < pval_thres & log2FoldChange >0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Enriched_ShootvsSoil")
df_all <- rbind(df_all,df_freq)

#Create structures for enrichments
df_freq <- res %>% subset(Contrast == "Shoot-Soil" & 
                            padj < pval_thres & log2FoldChange <0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Depleted_ShootvsSoil")
df_all <- rbind(df_all,df_freq)

df_freq <- res %>% subset(Contrast == "Root-Shoot" & 
                            padj < pval_thres & log2FoldChange >0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Enriched_RootvsShoot")
df_all <- rbind(df_all,df_freq)

#Create structures for enrichments
df_freq <- res %>% subset(Contrast == "Root-Shoot" & 
                            padj <= pval_thres & log2FoldChange <0) %$%
  Order %>% table %>% as.data.frame
df_freq$Perc <- df_freq$Freq/sum(df_freq$Freq)
df_freq$Type <- rep("Depleted_RootvsShoot")
df_all <- rbind(df_all,df_freq)


colnames(df_all)[1] <- "Order"
df_all$Facet <- rep("RootvsSoil",nrow(df_all))
df_all$Facet[df_all$Type %>% grep(pattern = "Shoot.*Soil")] <- "ShootvsSoil"
df_all$Facet[df_all$Type %>% grep(pattern = "Root.*Shoot")] <- "RootvsShoot"
df_all$Facet <- df_all$Facet %>% 
  factor(levels = c("RootvsSoil","ShootvsSoil","RootvsShoot"))
df_all$Mode <- rep("Enriched",nrow(df_all))
df_all$Mode[df_all$Type %>% grep(pattern = "Depleted")] <- "Depleted"
df_all$Mode <- df_all$Mode %>% factor(levels = c("Enriched","Depleted"))

size_ticks_x = 0
size_strip_text = 22
size_axis_text = 25
size_axis_title = 0
font_family = "Arial"
size_title_text = 35
y_vjust=0.5
width_bars =1
size_ticks_x = 2.5
size_ticks_y =2.5
legend_proportion_size = 0
size_legend_text = 0
spacing_x = 0.4




p <- ggplot(data = df_all,mapping = aes(Mode,y = Freq,fill = Order)) + 
  geom_bar(stat = "identity")+ scale_fill_manual(values = paleta) +
  coord_cartesian(expand = FALSE) + 
  facet_grid(~Facet,scales = "free_x",space = "free_x") +
  theme(legend.position = "none",aspect.ratio = 10) +
  theme(
    axis.line = element_blank(),
    panel.background = element_rect(fill = "white",color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.spacing.x = unit(spacing_x, "lines"),
    axis.ticks.y =element_line(colour = "black",size = size_ticks_y),
    axis.ticks.x =element_line(colour = "black",size = size_ticks_x),
    axis.text.x =element_blank(),
    axis.text.y = element_text(family = font_family,face="plain",size=size_axis_text,colour="#414141",vjust = y_vjust),
    axis.title.x = element_text(family = font_family,face="plain",size = size_axis_title,colour = "#414141"),
    axis.title.y = element_text(family = font_family,face="plain",size=size_axis_title,colour="#414141"),
    strip.background = element_blank()) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100,120,140))
oh.save.pdf(p = p, outdir = "../figures/",outname = "figure2_fungi_otus_phylogram_fraction_enrichments_depletions_counts.pdf")
