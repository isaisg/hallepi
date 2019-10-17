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

setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')



source('plotting_parameters_hallepi.R')

#Read bacterial otus dat
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")
Dat_rar <- Dat$Rarefied

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
# int_all <- read.table(file = "../rawdata/pots_otus_fungi.txt") %$% x %>%
#   as.character
# int_all <- int_all[-6]


#Root versus Shoot ###
samples_remove <- Dat_root$Map$DADA2_Header[which(!(Dat_root$Map$Pot %in% int_root_shoot))] %>%
  as.character
Dat_root <- remove_samples(Dat = Dat_root,samples = samples_remove,droplevels = T)
samples_remove <- Dat_shoot$Map$DADA2_Header[which(!(Dat_shoot$Map$Pot %in% int_root_shoot))] %>%
  as.character
Dat_shoot <- remove_samples(Dat = Dat_shoot,samples = samples_remove,droplevels = T)

#Order
Map_root <- Dat_root$Map
Map_shoot <- Dat_shoot$Map
Map_root <- match(int_root_shoot,Map_root$Pot) %>%
  Map_root[.,]
Map_shoot <- match(int_root_shoot,Map_shoot$Pot) %>%
  Map_shoot[.,]
Tab_root <- Dat_root$Tab
Tab_shoot <- Dat_shoot$Tab

Tab_root <- match(Map_root$DADA2_Header,colnames(Tab_root)) %>%
  Tab_root[,.]
Tab_shoot <- match(Map_shoot$DADA2_Header,colnames(Tab_shoot)) %>%
  Tab_shoot[,.]

#Define distfun to use
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")
dist_root <- distfun(x = Tab_root %>% t)
dist_shoot <- distfun(x = Tab_shoot %>% t)

m1 <- mantel(xdis = dist_root,ydis = dist_shoot,
             permutations = 10000)

capture.output(file = "../figures/supfigure2_mantel_fungi_root_shoot.doc",
               append = F,print(m1))

df1 <- melt_dist(dist = dist_root %>% as.matrix) 
df2 <- melt_dist(dist = dist_shoot %>% as.matrix)

dfpval <- data.frame(r = m1$statistic %>% round(3) ,p = m1$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Root = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_root)
p <- ggplot(data = melted,aes(Root,Shoot)) + geom_point(shape = 21,,size = 5) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.9,0.3,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_fungi_root_shoot_big.pdf",outdir = "../figures/",width = 10,height = 10)



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

dist_root <- distfun(x = Tab_root %>% t)
dist_shoot <- distfun(x = Tab_shoot %>% t)
dist_soil <- distfun(x = Tab_soil %>% t)

m_root_shoot <- mantel(xdis = dist_root,ydis = dist_shoot,
                       permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_fungi_root_shoot_small.doc",
               append = F,print(m_root_shoot))

df1 <- melt_dist(dist = dist_root %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_root_shoot$statistic %>% round(3) ,p = m_root_shoot$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Root = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_root)
p <- ggplot(data = melted,aes(Root,Shoot)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.9,0.3,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_fungi_root_shoot_small.pdf",outdir = "../figures/",width = 10,height = 10)



m_root_soil <- mantel(xdis = dist_root,ydis = dist_soil,
                      permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_fungi_root_soil_small.doc",
               append = F,print(m_root_soil))
df1 <- melt_dist(dist = dist_soil %>% as.matrix) 

df2 <- melt_dist(dist = dist_root %>% as.matrix) 

dfpval <- data.frame(r = m_root_soil$statistic %>% round(3) ,p = m_root_soil$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Soil = df1$dist,Root = df2$dist)
mlm <- lm(dist_root~dist_soil)
p <- ggplot(data = melted,aes(Soil,Root)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.7,0.4,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_fungi_soil_root_small.pdf",outdir = "../figures/",width = 10,height = 10)


m_shoot_soil <- mantel(xdis = dist_shoot,ydis = dist_soil,
                       permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_fungi_shoot_soil_small.doc",
               append = F,print(m_shoot_soil))


df1 <- melt_dist(dist = dist_soil %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_shoot_soil$statistic %>% round(3) ,p = m_shoot_soil$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Soil = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_soil)
p <- ggplot(data = melted,aes(Soil,Shoot)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.65,0.4,label =label),size = 5)
oh.save.pdf(p = p,outname = "mantel_fungi_soil_shoot_small.pdf",outdir = "../figures/",width = 10,height = 10)


### Check the cloud between root and soil
m_root_soil <- mantel(xdis = dist_root,ydis = dist_soil,
                      permutations = 10000)

df1 <- melt_dist(dist = dist_soil %>% as.matrix) 

df2 <- melt_dist(dist = dist_root %>% as.matrix) 

dfpval <- data.frame(r = m_root_soil$statistic %>% round(3) ,p = m_root_soil$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(S1 = df1$iso1,S2 = df2$iso2,Soil = df1$dist,Root = df2$dist)
mlm <- lm(dist_root~dist_soil)


#Determine which ones are over 0.7 distance in soil and remove them 
df1 <- melt_dist(dist = dist_soil %>% as.matrix) 
temp <- df1 %>% subset(dist > 0.7)
temp$iso1 %>% table %>% sort
temp$iso2 %>% table %>% sort

toremove <- c("Plate3D3","Plate3B3","Plate3E5","Plate3B4")


#Remove bulksoil samples
Dat_root <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Root",drop = T, clean =T)

Dat_soil <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Soil",drop = T, clean =T)

Dat_shoot <- subset.Dataset(x = Dat_rar,Genotype!="No_Plant",drop = T,clean = T) %>%
  subset.Dataset(x = ., Fraction == "Shoot",drop = T, clean =T)

#Determine share pots
which(Dat_soil$Map$Pot %in% toremove) %>%
  Dat_soil$Map[.,]
Dat_soil <- AMOR::remove_samples(Dat = Dat_soil,samples = toremove,droplevels = T)


pots_root <- Dat_root$Map$Pot %>% unique %>% as.character
pots_soil <- Dat_soil$Map$Pot %>% unique %>% as.character
pots_shoot <- Dat_shoot$Map$Pot %>% unique %>% as.character

int_root_shoot <- intersect(pots_root,pots_shoot)
int_all <- intersect(int_root_shoot,pots_soil)

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

dist_root <- distfun(x = Tab_root %>% t)
dist_shoot <- distfun(x = Tab_shoot %>% t)
dist_soil <- distfun(x = Tab_soil %>% t)

m_root_shoot <- mantel(xdis = dist_root,ydis = dist_shoot,
                       permutations = 10000)


df1 <- melt_dist(dist = dist_root %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_root_shoot$statistic %>% round(3) ,p = m_root_shoot$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Root = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_root)
p <- ggplot(data = melted,aes(Root,Shoot)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.9,0.3,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_fungi_root_shoot_small.withoutscloud.pdf",outdir = "../figures/",width = 10,height = 10)





m_root_soil <- mantel(xdis = dist_root,ydis = dist_soil,
                      permutations = 10000)
df1 <- melt_dist(dist = dist_soil %>% as.matrix) 

df2 <- melt_dist(dist = dist_root %>% as.matrix) 

dfpval <- data.frame(r = m_root_soil$statistic %>% round(3) ,p = m_root_soil$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Soil = df1$dist,Root = df2$dist)
mlm <- lm(dist_root~dist_soil)
p <- ggplot(data = melted,aes(Soil,Root)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.7,0.4,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_fungi_soil_root_small.withoutscloud.pdf",
            outdir = "../figures/",width = 10,height = 10)



m_shoot_soil <- mantel(xdis = dist_shoot,ydis = dist_soil,
                       permutations = 10000)


df1 <- melt_dist(dist = dist_soil %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_shoot_soil$statistic %>% round(3) ,p = m_shoot_soil$signif %>% format.pval())
dfpval$label <- paste0("r= ",dfpval$r,"\npval=",dfpval$p)
melted <- data.frame(Soil = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_soil)
p <- ggplot(data = melted,aes(Soil,Shoot)) + geom_point(shape = 21,,size = 8) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.65,0.4,label =label),size = 5)
oh.save.pdf(p = p,outname = "mantel_fungi_soil_shoot_small.withoutscloud.pdf",outdir = "../figures/",width = 10,height = 10)

