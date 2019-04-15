library(ohchibi)
library(DESeq2)


Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_raw <- Dat_ori$RawCounts

Dat_raw <- Dat_raw %>% 
  subset.Dataset(x = .,subset = Tissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(x = .,subset =  typebyTissue != "Agar_NP",drop = T,clean = T)

Dat_raw$Map$Phosphate <- Dat_raw$Map$condition %>% 
  gsub(pattern = " micromolar ",replacement = "") %>%
  factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))

Dat_raw <- Dat_raw %>% subset.Dataset(subset = Phosphate != "50Pi",drop = T,clean = T) %>%
  subset.Dataset(subset = Phosphate != "30Pi",drop = T,clean = T) 

Dat_raw$Map$PhosphateInt <- Dat_raw$Map$Phosphate %>% as.character
Dat_raw$Map$PhosphateInt[which(Dat_raw$Map$PhosphateInt %in% c("100Pi","1000Pi"))] <- "High"
Dat_raw$Map$PhosphateInt[which(Dat_raw$Map$PhosphateInt %in% c("0Pi","10Pi"))] <- "Low"
Dat_raw$Map$PhosphateInt <- Dat_raw$Map$PhosphateInt %>% factor(levels = c("High","Low"))

#Create Deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design =~ Rep + Tissue + PhosphateInt + Tissue:PhosphateInt)
dds<-DESeq(object = dds)

df_root <- dds %>% results(name = "TissueRoot.PhosphateIntLow") %>%
  as.data.frame
df_root$Contrast <- rep("TissueRoot.PhosphateIntLow",nrow(df_root))
df_root$Id <- rownames(df_root)
df_root$Facet <- rep("Phosphate Interaction",nrow(df_root))
df_root$Specific <- rep("Root",nrow(df_root))

df_shoot <- dds %>% results(name = "TissueShoot.PhosphateIntLow") %>%
  as.data.frame
df_shoot$Contrast <- rep("TissueShoot.PhosphateIntLow",nrow(df_shoot))
df_shoot$Id <- rownames(df_shoot)
df_shoot$Facet <- rep("Phosphate Interaction",nrow(df_shoot))
df_shoot$Specific <- rep("Shoot",nrow(df_shoot))

res <- rbind(df_root,df_shoot)

res <- res$Id %>% grep(pattern = "OTU",invert = TRUE) %>%
  res[.,] %>% droplevels


write.table(x = res,file = "../cleandata/df_dds_res_syncom_phosphateint.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
