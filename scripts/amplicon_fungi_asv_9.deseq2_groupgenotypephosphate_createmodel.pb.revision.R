library(ohchibi)
library(DESeq2)



set.seed(seed = 130816)
#Load dataset
Dat <- readRDS(file = "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS")

####################################
############ Root ##################

##### ASV Level #####
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

original_reads <- Dat_raw$Tab %>% sum

Dat_raw <- measurable_taxa(Dat = Dat_raw,min_samples_otu = 3,min_reads_otu = 12,
                           method = "absolute",clean = T)

filtered_reads <- Dat_raw$Tab %>% sum
(filtered_reads/original_reads)*100

Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor



#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_otu_root <- dds

###### Genus Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 7)

Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor

#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_genus_root  <- dds

###### Family Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)


#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 6)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor

#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_family_root  <- dds

###### Order Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 5)

Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_order_root  <- dds

###### Class Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 4)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor

#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_class_root  <- dds

###### Phylum Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Root",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 3)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor

#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_phylum_root  <- dds



####################################
############ Shoot #################

##### ASVS Level #####
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

original_reads <- Dat_raw$Tab %>% sum

Dat_raw <- measurable_taxa(Dat = Dat_raw,min_samples_otu = 2,min_reads_otu = 12,method = "absolute",clean = T)

filtered_reads <- Dat_raw$Tab %>% sum
(filtered_reads/original_reads)*100

Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_otu_shoot <- dds

###### Genus Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 7)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_genus_shoot  <- dds

###### Family Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 6)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_family_shoot  <- dds

###### Order Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 5)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_order_shoot  <- dds

###### Class Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 4)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_class_shoot  <- dds

###### Phylum Level #####
rm(dds)
Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction=="Shoot",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype!="No_Plant",drop = T,clean = T)  %>%
  #subset.Dataset(x = .,Genotype != "phr1/phl1",drop = T,clean = T) %>%
  subset.Dataset(x = . ,Phosphate == "low" | Phosphate == "low+Pi",drop = T,clean = T)

#Collapse to the genus level
Dat_raw <- AMOR::collapse_by_taxonomy.Dataset(Dat = Dat_raw,level = 3)
Dat_raw$Map$Phosphate <- Dat_raw$Map$Phosphate %>% gsub(pattern = "\\+",replacement = "_") %>%
  factor

Dat_raw$Map$Genotype <- Dat_raw$Map$Genotype %>% gsub(pattern = "\\/",replacement = "_") %>%
  gsub(pattern = "\\-",replacement = "_") %>% factor
#Create group variable
Dat_raw$Map$group <- paste(Dat_raw$Map$Genotype,Dat_raw$Map$Phosphate,sep = "_") %>%
  factor
#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~Plot + group)
dds <- DESeq(dds)
dds_phylum_shoot  <- dds



##Save results as lists
Root <- list(ASV = dds_otu_root,Genus = dds_genus_root,Family = dds_family_root,
             Order = dds_order_root,Class = dds_class_root,Phylum = dds_phylum_root)
Shoot <- list(ASV = dds_otu_shoot,Genus = dds_genus_shoot, Family = dds_family_shoot,
              Order = dds_order_shoot,Class = dds_class_shoot, Phylum = dds_phylum_shoot)



dds_res <- list(Root = Root, Shoot = Shoot)
saveRDS(object = dds_res,file = "../cleandata/dds_deseq2_hallepi_amplicon_fungi_asvs_groupsgenotypephosphate.RDS")
