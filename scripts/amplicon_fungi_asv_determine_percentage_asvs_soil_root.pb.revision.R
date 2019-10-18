library(ohchibi)
library(palettesPM)
library(paletteer)
library(extrafont)
loadfonts(device = "pdf")
library(egg)
library(harrietr)
#Set random seed
set.seed(130816)


size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4




source('plotting_parameters_hallepi.R')

#Read bacterial otus dat
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$RawCounts

#Determine not problematic pots
Map <- Dat_rar$Map
Map$Count <- rep(1,nrow(Map))
df_freq <- acast(data = Map,formula = Pot~Fraction,fun.aggregate = sum,
                 value.var = "Count")

res <- apply(X = df_freq,MARGIN = 1,
             FUN = function(x)which(x>1)%>% length)
noms <- which(res == 1) %>% names
Map <- Map[which(Map$Pot %in% noms),] %>% droplevels
samples_remove <- Map$DADA2_Header %>% as.character
Dat_rar <- remove_samples(Dat = Dat_rar,samples = samples_remove,droplevels = T)


#Remove bulksoil samples
Dat_root <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root",drop = T, clean =T)

Dat_soil <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Soil",drop = T, clean =T)

Dat_shoot <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot",drop = T, clean =T)

#Determine share pots
pots_root <- Dat_root$Map$Pot %>% unique %>% as.character
pots_soil <- Dat_soil$Map$Pot %>% unique %>% as.character
pots_shoot <- Dat_shoot$Map$Pot %>% unique %>% as.character

int_root_shoot <- intersect(pots_root,pots_shoot)
int_all <- intersect(int_root_shoot,pots_soil)
#int_all <- intersect(pots_root,pots_soil)


## Add soil
Map_root <- Dat_root$Map
Map_shoot <- Dat_shoot$Map
Map_soil <- Dat_soil$Map

Map_root <- match(int_all,Map_root$Pot) %>%
  Map_root[.,]
Map_shoot <- match(int_all,Map_shoot$Pot) %>%
  Map_shoot[.,]
Map_soil <- match(int_all,Map_soil$Pot) %>%
  Map_soil[.,]

Tab_root <- Dat_root$Tab
Tab_shoot <- Dat_shoot$Tab
Tab_soil <- Dat_soil$Tab

Tab_root <- match(Map_root$DADA2_Header,colnames(Tab_root)) %>%
  Tab_root[,.]
Tab_shoot <- match(Map_shoot$DADA2_Header,colnames(Tab_shoot)) %>%
  Tab_shoot[,.]
Tab_soil <- match(Map_soil$DADA2_Header,colnames(Tab_soil)) %>%
  Tab_soil[,.]

#Define function to use
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

dist_root <- distfun(x = Tab_root %>% t)
dist_shoot <- distfun(x = Tab_shoot %>% t)
dist_soil <- distfun(x = Tab_soil %>% t)

#New version
melted_soil <- Tab_soil %>% melt
melted_root <- Tab_root %>% melt
melted_shoot <- Tab_shoot %>% melt

colnames(melted_soil) <- c("Id","DADA2_Header","value")
colnames(melted_root) <- c("Id","DADA2_Header","value")
colnames(melted_shoot) <- c("Id","DADA2_Header","value")


melted_soil <- merge(melted_soil,Dat$RawCounts$Map, by = "DADA2_Header")  %>% droplevels
melted_root <- merge(melted_root,Dat$RawCounts$Map, by = "DADA2_Header")  %>% droplevels
melted_shoot <- merge(melted_shoot,Dat$RawCounts$Map, by = "DADA2_Header") %>% droplevels


melted <- rbind(melted_soil,melted_root,melted_shoot)
allTab <- acast(data = melted,formula = Id~Pot+Fraction,value.var = "value")
allTab[allTab>1] <- 1


merged <- allTab %>% melt
colnames(merged) <- c("Id","PotbyFraction","value")
merged$Pot <- merged$PotbyFraction %>%
  gsub(pattern = "_.*",replacement = "") %>% factor
merged$Fraction <- merged$PotbyFraction %>%
  gsub(pattern = ".*_",replacement = "") %>% factor

#Loop over the pots
Res <- NULL
for(pot in merged$Pot %>% unique){
  df_sub <- merged %>% subset(Pot == pot) %>% droplevels %>%
    dcast(formula = Id~Fraction,value.var = "value")
  
  total_soil_root <- df_sub %>% subset(Root == 1 | Soil == 1) %>% nrow
  shared_soil_root <- df_sub %>% subset(Root == 1 &  Soil == 1) %>% nrow
  perc_share_soil_root <- (shared_soil_root/total_soil_root)*100
  
  total_soil_shoot <- df_sub %>% subset(Shoot == 1 | Soil == 1) %>% nrow
  shared_soil_shoot <- df_sub %>% subset(Shoot == 1 &  Soil == 1) %>% nrow
  perc_share_soil_shoot <- (shared_soil_shoot/total_soil_shoot)*100
  
  total_root_shoot <- df_sub %>% subset(Shoot == 1 | Root == 1) %>% nrow
  shared_root_shoot <- df_sub %>% subset(Shoot == 1 &  Root == 1) %>% nrow
  perc_share_root_shoot <- (shared_root_shoot/total_root_shoot)*100
  
  Res <-  data.frame( PercSoilRoot = perc_share_soil_root,
             PercSoilShoot = perc_share_soil_shoot,
             PercRootShoot = perc_share_root_shoot) %>%
    rbind(Res,.)
}
melted  <- Res %>% melt
Res$PercSoilRoot %>% mean

rm(list=ls())