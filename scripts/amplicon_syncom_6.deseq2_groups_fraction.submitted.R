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

### Group model ###
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,design = ~ 0 + group + Rep)
dds <- DESeq(dds)

## Root vs Agar ###
df_ra0 <- dds %>% results(contrast = c("group","Root0Pi","Agar0Pi")) %>%
  as.data.frame
df_ra0$Contrast <- rep("Root0Pi_vs_Agar0Pi",nrow(df_ra0))
df_ra0$Id <- rownames(df_ra0)

df_ra10 <- dds %>% results(contrast = c("group","Root10Pi","Agar10Pi")) %>%
  as.data.frame
df_ra10$Contrast <- rep("Root10Pi_vs_Agar10Pi",nrow(df_ra10))
df_ra10$Id <- rownames(df_ra10)

df_ra30 <- dds %>% results(contrast = c("group","Root30Pi","Agar30Pi")) %>%
  as.data.frame
df_ra30$Contrast <- rep("Root30Pi_vs_Agar30Pi",nrow(df_ra30))
df_ra30$Id <- rownames(df_ra30)

df_ra50 <- dds %>% results(contrast = c("group","Root50Pi","Agar50Pi")) %>%
  as.data.frame
df_ra50$Contrast <- rep("Root50Pi_vs_Agar50Pi",nrow(df_ra50))
df_ra50$Id <- rownames(df_ra50)

df_ra100 <- dds %>% results(contrast = c("group","Root100Pi","Agar100Pi")) %>%
  as.data.frame
df_ra100$Contrast <- rep("Root100Pi_vs_Agar100Pi",nrow(df_ra100))
df_ra100$Id <- rownames(df_ra100)

df_ra1000 <- dds %>% results(contrast = c("group","Root1000Pi","Agar1000Pi")) %>%
  as.data.frame
df_ra1000$Contrast <- rep("Root1000Pi_vs_Agar1000Pi",nrow(df_ra1000))
df_ra1000$Id <- rownames(df_ra1000)


## Shoot vs Agar ###
df_sa0 <- dds %>% results(contrast = c("group","Shoot0Pi","Agar0Pi")) %>%
  as.data.frame
df_sa0$Contrast <- rep("Shoot0Pi_vs_Agar0Pi",nrow(df_sa0))
df_sa0$Id <- rownames(df_sa0)

df_sa10 <- dds %>% results(contrast = c("group","Shoot10Pi","Agar10Pi")) %>%
  as.data.frame
df_sa10$Contrast <- rep("Shoot10Pi_vs_Agar10Pi",nrow(df_sa10))
df_sa10$Id <- rownames(df_sa10)

df_sa30 <- dds %>% results(contrast = c("group","Shoot30Pi","Agar30Pi")) %>%
  as.data.frame
df_sa30$Contrast <- rep("Shoot30Pi_vs_Agar30Pi",nrow(df_sa30))
df_sa30$Id <- rownames(df_sa30)

df_sa50 <- dds %>% results(contrast = c("group","Shoot50Pi","Agar50Pi")) %>%
  as.data.frame
df_sa50$Contrast <- rep("Shoot50Pi_vs_Agar50Pi",nrow(df_sa50))
df_sa50$Id <- rownames(df_sa50)

df_sa100 <- dds %>% results(contrast = c("group","Shoot100Pi","Agar100Pi")) %>%
  as.data.frame
df_sa100$Contrast <- rep("Shoot100Pi_vs_Agar100Pi",nrow(df_sa100))
df_sa100$Id <- rownames(df_sa100)

df_sa1000 <- dds %>% results(contrast = c("group","Shoot1000Pi","Agar1000Pi")) %>%
  as.data.frame
df_sa1000$Contrast <- rep("Shoot1000Pi_vs_Agar1000Pi",nrow(df_sa1000))
df_sa1000$Id <- rownames(df_sa1000)


## Shoot vs Root ##
df_sr0 <- dds %>% results(contrast = c("group","Root0Pi","Shoot0Pi")) %>%
  as.data.frame
df_sr0$Contrast <- rep("Root0Pi_vs_Shoot0Pi",nrow(df_sr0))
df_sr0$Id <- rownames(df_sr0)

df_sr10 <- dds %>% results(contrast = c("group","Root10Pi","Shoot10Pi")) %>%
  as.data.frame
df_sr10$Contrast <- rep("Root10Pi_vs_Shoot10Pi",nrow(df_sr10))
df_sr10$Id <- rownames(df_sr10)

df_sr30 <- dds %>% results(contrast = c("group","Root30Pi","Shoot30Pi")) %>%
  as.data.frame
df_sr30$Contrast <- rep("Root30Pi_vs_Shoot30Pi",nrow(df_sr30))
df_sr30$Id <- rownames(df_sr30)

df_sr50 <- dds %>% results(contrast = c("group","Root50Pi","Shoot50Pi")) %>%
  as.data.frame
df_sr50$Contrast <- rep("Root50Pi_vs_Shoot50Pi",nrow(df_sr50))
df_sr50$Id <- rownames(df_sr50)

df_sr100 <- dds %>% results(contrast = c("group","Root100Pi","Shoot100Pi")) %>%
  as.data.frame
df_sr100$Contrast <- rep("Root100Pi_vs_Shoot100Pi",nrow(df_sr100))
df_sr100$Id <- rownames(df_sr100)

df_sr1000 <- dds %>% results(contrast = c("group","Root1000Pi","Shoot1000Pi")) %>%
  as.data.frame
df_sr1000$Contrast <- rep("Root1000Pi_vs_Shoot1000Pi",nrow(df_sr1000))
df_sr1000$Id <- rownames(df_sr1000)

#Merged results
res <- rbind(df_ra0,df_ra10,df_ra30,df_ra50,df_ra100,df_ra1000,
             df_sa0,df_sa10,df_sa30,df_sa50,df_sa100,df_sa1000,
             df_sr0,df_sr10,df_sr30,df_sr50,df_sr100,df_sr1000)

#res <- res$Id %>% grep(pattern = "OTU",invert = TRUE) %>%
# res[.,] %>% droplevels

res$Contrast <- res$Contrast %>%
  factor(levels =c("Root0Pi_vs_Agar0Pi","Root10Pi_vs_Agar10Pi","Root30Pi_vs_Agar30Pi",
                   "Root50Pi_vs_Agar50Pi","Root100Pi_vs_Agar100Pi","Root1000Pi_vs_Agar1000Pi",
                   "Shoot0Pi_vs_Agar0Pi","Shoot10Pi_vs_Agar10Pi","Shoot30Pi_vs_Agar30Pi",
                   "Shoot50Pi_vs_Agar50Pi","Shoot100Pi_vs_Agar100Pi","Shoot1000Pi_vs_Agar1000Pi",
                   "Root0Pi_vs_Shoot0Pi","Root10Pi_vs_Shoot10Pi","Root30Pi_vs_Shoot30Pi",
                   "Root50Pi_vs_Shoot50Pi","Root100Pi_vs_Shoot100Pi","Root1000Pi_vs_Shoot1000Pi") )
rownames(res) <- NULL

res$Facet <- res$Contrast %>% as.character
res$Facet <- res$Facet %>% gsub(pattern = "Root.*Agar.*",replacement = "Root vs Agar") %>%
  gsub(pattern = "Shoot.*Agar.*",replacement = "Shoot vs Agar") %>%
  gsub(pattern = "Root.*Shoot.*",replacement = "Root vs Shoot") %>%
  factor()
res$Specific <- res$Contrast %>% gsub(pattern = ".*_vs_[a-z|A-Z]+",replacement = "") %>%
  factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))

write.table(x = res,file = "../cleandata/df_dds_res_syncom_fraction_specific.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)



#### Now compute main effects #######
lista <- list(c("groupRoot0Pi","groupRoot10Pi","groupRoot30Pi","groupRoot50Pi","groupRoot100Pi","groupRoot1000Pi"),
              c("groupAgar0Pi","groupAgar10Pi","groupAgar30Pi","groupAgar50Pi","groupAgar100Pi","groupAgar1000Pi"))
res_global_root_agar <- results(object = dds,contrast = lista,listValues = c(1/6,-1/6)) %>%
  as.data.frame  

lista <- list(c("groupShoot0Pi","groupShoot10Pi","groupShoot30Pi","groupShoot50Pi","groupShoot100Pi","groupShoot1000Pi"),
              c("groupAgar0Pi","groupAgar10Pi","groupAgar30Pi","groupAgar50Pi","groupAgar100Pi","groupAgar1000Pi"))
res_global_shoot_agar <- results(object = dds,contrast = lista,listValues = c(1/6,-1/6)) %>%
  as.data.frame  


lista <- list(c("groupRoot0Pi","groupRoot10Pi","groupRoot30Pi","groupRoot50Pi","groupRoot100Pi","groupRoot1000Pi"),
              c("groupShoot0Pi","groupShoot10Pi","groupShoot30Pi","groupShoot50Pi","groupShoot100Pi","groupShoot1000Pi"))
res_global_root_shoot <- results(object = dds,contrast = lista,listValues = c(1/6,-1/6)) %>%
  as.data.frame  

#Add contrast
res_global_root_agar$Contrast <- "Global_Root_vs_Agar"
res_global_shoot_agar$Contrast <- "Global_Shoot_vs_Agar"
res_global_root_shoot$Contrast <- "Global_Root_vs_Shoot"

#Add Ids
res_global_root_agar$Id <- rownames(res_global_root_agar)
res_global_shoot_agar$Id <- rownames(res_global_shoot_agar)
res_global_root_shoot$Id <- rownames(res_global_root_shoot)

#Add Facet
res_global_root_agar$Facet <- rep("Root vs Agar",nrow(res_global_root_agar))
res_global_shoot_agar$Facet <- rep("Shoot vs Agar",nrow(res_global_shoot_agar))
res_global_root_shoot$Facet <- rep("Root vs Shoot",nrow(res_global_root_shoot))



res_global <- rbind(res_global_root_agar,res_global_shoot_agar,res_global_root_shoot)
rownames(res_global) <- NULL
write.table(x = res_global,file = "../cleandata/df_dds_res_syncom_fraction_global.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
