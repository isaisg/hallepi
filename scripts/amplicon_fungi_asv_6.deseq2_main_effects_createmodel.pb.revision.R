library(ohchibi)
library(DESeq2)

set.seed(seed = 130816)
setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')

#Load dataset
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")

##### ASV Level #####
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",
                 drop = T,clean = T)  

#Explore the diversity of sequences
df_freq <- Dat_raw$Tab %>% rowSums %>% 
  data.frame(Id = names(.), Freq = .,row.names = NULL)
df_present <- Dat_raw$Tab %>% 
  apply(X = .,MARGIN = 1,FUN = function(x)which(x>0)%>% length)%>%
  data.frame(Id = names(.),PresentIn = .,row.names = NULL)
df_merged <- merge(df_freq,df_present,by = "Id")
df_merged$Average <- df_merged$Freq/df_merged$PresentIn

ggplot(df_merged,aes(PresentIn,Average)) + geom_point() +
  scale_y_continuous(limits = c(0,3000))


original_reads <- Dat_raw$Tab %>% sum
#Identify measurable ESVs
Dat_raw <- measurable_taxa(Dat = Dat_raw,min_samples_otu =10,min_reads_otu = 25,method = "absolute",clean = T)
filtered_reads <- Dat_raw$Tab %>% sum

(filtered_reads/original_reads)*100

df_freq <- Dat_raw$Tab %>% rowSums %>% 
  data.frame(Id = names(.), Freq = .,row.names = NULL)
df_present <- Dat_raw$Tab %>% 
  apply(X = .,MARGIN = 1,FUN = function(x)which(x>0)%>% length)%>%
  data.frame(Id = names(.),PresentIn = .,row.names = NULL)
df_merged <- merge(df_freq,df_present,by = "Id")
df_merged$Average <- df_merged$Freq/df_merged$PresentIn

ggplot(df_merged,aes(PresentIn,Average)) + geom_point() +
  scale_y_continuous(limits = c(0,500))


#Write the chosen esvs
Dat_raw$Tax$ID %>% as.character  %>% 
  write(x = .,file = "../cleandata/386_asvs.txt")


#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~ Fraction + Genotype + Phosphate)


dds <- DESeq(dds)
dds_otus <- dds


##### Genus Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant" ,drop = T,clean = T) 

original_reads <- Dat_raw$Tab %>% sum

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 7)

dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Fraction + Genotype + Phosphate)


dds <- DESeq(dds)
dds_genus <- dds

##### Family Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T) 

original_reads <- Dat_raw$Tab %>% sum

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 6)

dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Fraction + Genotype + Phosphate)

dds <- DESeq(dds)

dds_family <- dds

##### Order Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant" ,drop = T,clean = T) 

original_reads <- Dat_raw$Tab %>% sum

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 5)


dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Fraction + Genotype + Phosphate)

dds <- DESeq(dds)

dds_order <- dds

###### Class Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T) 

original_reads <- Dat_raw$Tab %>% sum

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 4)

dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Fraction + Genotype + Phosphate)

dds <- DESeq(dds)

dds_class <- dds



##### Phylum Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant" ,drop = T,clean = T) 

original_reads <- Dat_raw$Tab %>% sum

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 3)


dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Fraction + Genotype + Phosphate)

dds <- DESeq(dds)

dds_phylum <- dds
rm(dds)


dds_res <- list(ASV = dds_otus, Genus = dds_genus, Family = dds_family,
                Order = dds_order, Class = dds_class, Phylum = dds_phylum)
saveRDS(object = dds_res,file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_asvs_maineffects.RDS")
