library(ohchibi)
library(DESeq2)

set.seed(130816)

Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_raw <- Dat_ori$RawCounts

Dat_raw <- Dat_raw %>% 
  subset.Dataset(x = .,subset = Tissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(x = .,subset =  typebyTissue != "Agar_NP",drop = T,clean = T)

Dat_raw$Map$Phosphate <- Dat_raw$Map$condition %>% 
  gsub(pattern = " micromolar ",replacement = "") %>%
  factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))

Dat_raw$Map$group <- paste(Dat_raw$Map$Tissue,Dat_raw$Map$Phosphate,sep = "") %>%
  factor
Dat_raw$Map$Phosphate <- Dat_raw$Map$condition %>%
  gsub(pattern = " micromolar ",replacement = "") %>%
  factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))


# Full model
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                                   design = ~Rep + Tissue + Phosphate)
dds <- DESeq(dds)

res_ra <- results(dds,contrast = c("Tissue","Root","Agar")) %>%
  as.data.frame
res_sa <- results(dds,contrast = c("Tissue","Shoot","Agar")) %>%
  as.data.frame

#Add ids to aggregate
res_ra$Id <- rownames(res_ra)
res_sa$Id <- rownames(res_sa)

res_ra$Fraction <- rep("Root",nrow(res_ra))
res_sa$Fraction <- rep("Shoot",nrow(res_sa))

res <- rbind(res_ra,res_sa)
rownames(res) <- NULL

#Aggregate
df_lf <- aggregate(log2FoldChange ~ Id,res,mean)
df_pval <- aggregate(padj ~ Id,res,mean)

merged <- merge(df_lf,df_pval, by = "Id")

#Count
era <- res_ra %>% subset(padj < 0.1) %$% Id
esa <- res_sa %>% subset(padj < 0.1) %$% Id
all <- union(era,esa)
all <- all %>% grep(pattern = "OTU",value = T,invert = T)

#Load relative abundance
Dat <- readRDS(file = "../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat <- Dat$RelativeAbundance
melted <- Dat$Tab %>% melt
colnames(melted) <- c("Id","ID_Matrix","RA")
melted <- merge(melted,Dat$Map,by = "ID_Matrix")
df <- dcast(data = melted,formula = Id~typebyTissue,
      fun.aggregate = mean,fill = 0,value.var = "RA")

df_all <- match(all,df$Id) %>% df[.,] 
colSums(df_all[,2:ncol(df_all)] %>% as.matrix)

#Now ask the opposite question
chosen <- df %>% subset(Root_SC >= 0.01) %$% Id %>%
  as.character

res_ra_chosen <- match(chosen,rownames(res_ra)) %>%
  res_ra[.,] %>% droplevels
res_sa_chosen <- match(chosen,rownames(res_sa)) %>%
  res_sa[.,] %>% droplevels

chosen_ra <- res_ra_chosen %>% subset(padj < 0.1) %>% rownames
chosen_sa <- res_sa_chosen %>% subset(padj < 0.1) %>% rownames

which(!(chosen %in% chosen_ra)) %>% chosen[.]
which(!(chosen %in% chosen_sa)) %>% chosen[.]

#Make a figure to represent this
merged <- merge(merged,df, by = "Id") 
merged$RA_Plant <- (merged$Root_SC + merged$Shoot_SC) /2
merged <- merged$Id %>% grep(pattern = "OTU",invert = T) %>%
  merged[.,] %>% droplevels
merged$Significance <- rep("NoSignificant")
merged$Significance[which(merged$Id %in% all)] <- "Significant"
merged$Significance <- merged$Significance %>% factor
merged$padjtrans <- -log10(merged$padj)
merged$Type <- rep("Less",nrow(merged))
merged$Type[which(merged$RA_Plant >= 0.01)] <- "More"
merged$Type <- merged$Type %>% factor
merged$RA_Plant %>% max
p <- ggplot(data = merged,aes(padjtrans,log2FoldChange)) +
  geom_point(aes(shape = Type,fill = Significance,size = RA_Plant),stroke = 1.5) +
  scale_shape_manual(values = c(23,21)) +
  coord_flip() + 
  scale_x_sqrt() +
  theme_ohchibi(legend_proportion_size = 1) + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
  ) +
  scale_fill_manual(values = c("black","red")) +
  xlab(label = "-log10 q.value")
