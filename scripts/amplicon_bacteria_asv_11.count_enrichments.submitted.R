library(ohchibi)
library(palettesPM)

pval_thres <- 0.1
res <- read.table(file = "../cleandata/df_dds_res_amplicon_bacteria_asv_fraction_asvlevel.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")


#Read phylum information
mphylum <- palettesPM::pm.names.phyla()
res$Phylum <- res$Phylum %>% as.character
res$Phylum[which(!(res$Phylum %in% mphylum))] <- "Other"
res$Phylum <- res$Phylum %>% factor(levels = c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                               "Chloroflexi","Cyanobacteria","Firmicutes",
                                               "Gemmatimonadetes","Patescibacteria","Proteobacteria",
                                               "Verrucomicrobia","Other"))

ra <- res %>% subset(Contrast == "Root-Soil" & 
                 padj < pval_thres ) %$% ASV_Id %>% as.character
sa <- res %>% subset(Contrast == "Shoot-Soil" & 
                       padj < pval_thres ) %$% ASV_Id %>% as.character
all <- union(ra,sa) 
(all %>% length) / 3874

#load cumulative RA
Dat_asv <- readRDS(file = "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
melted_ab_all <- Dat_asv$RelativeAbundance$Tab %>% melt
colnames(melted_ab_all) <- c("Id","DADA2_Header","RA")
map_ab <- Dat_asv$RelativeAbundance$Map 
melted_ab_all <- merge(melted_ab_all,map_ab, by = "DADA2_Header")
#Compute average RE in fraction 
melted_ab <- dcast(data = melted_ab_all,formula = Id ~ Fraction,fun.aggregate = mean,value.var = "RA")


match(all,melted_ab$Id) %>% melted_ab[.,] %$% Root %>% sum
match(all,melted_ab$Id) %>% melted_ab[.,] %$% Shoot %>% sum
