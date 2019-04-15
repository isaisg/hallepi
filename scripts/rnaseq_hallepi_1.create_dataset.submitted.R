library(ohchibi)

#Set random seed
set.seed(130816)



#Read Tab and Map files
Tab <- read.table(file = "../rawdata/HALLE_mRNA_counts_filtered.txt",
                  header = T,row.names = 1,sep = "\t",quote = "",comment.char = "",check.names = F)
Map <- read.table(file = "../rawdata/RNASeqMetadata.txt",
                  header = T,sep = "\t",quote = "",comment.char = "")
Map$sample <- Map$sample %>% as.character() %>%
  gsub(pattern = "0+",replacement = "")

#Remove sample 43 clear outlier
Map <- Map %>% subset(sample != "HALLE43") %>% droplevels

#Create a dataframe with the length of the genes
df_length_genes <- data.frame(gene_id =rownames(Tab),
                              length_gene = Tab$Length)

Tab <- Tab[,-1]
ids <- colnames(Tab) %>% intersect(Map$sample)
Tab <- match(ids,colnames(Tab)) %>% Tab[,.]
Map <- match(ids,Map$sample) %>% Map[.,] %>% droplevels
rownames(Map) <- Map$sample

#Remove genes without counts
toremove <- Tab %>% rowSums %>% 
  data.frame(Gene = names(.), Freq = .,row.names = NULL) %>%
  subset(Freq == 0) %$% Gene %>% as.character

Tab <- which(!(rownames(Tab) %in% toremove)) %>% Tab[.,]

#Clean the Map and modify variables to keep consistnecy
colnames(Map)[4] <- "Fraction"
colnames(Map)[5] <- "Phosphate"
colnames(Map)[7] <- "Plot"

Map$Plot <- Map$Plot %>% paste0("Plot",.) %>% factor

Map$Fraction <- Map$Fraction %>% gsub(pattern = "root",replacement = "Root") %>%
  gsub(pattern = "shoot",replacement = "Shoot") %>%
  factor(levels = c("Root","Shoot"))
Map$Phosphate <- Map$Phosphate %>% 
  gsub(pattern = "[abc]-",replacement = "") %>%
  gsub(pattern = "\\+",replacement = "_") %>%
  factor(levels = c("low","medium","high","low_Pi"))
Map$Genotype <- Map$Genotype %>% gsub(pattern = "col",replacement = "Col") %>%
  gsub(pattern = "phf",replacement = "phf1") %>%
  gsub(pattern = "phr/phl",replacement = "phr1/phl1") %>%
  gsub(pattern = "\\/", replacement = "_") %>%
  gsub(pattern = "-",replacement = "_")
  factor(c("Col_0","phf1","phr1_phl1"))

Map$group <- paste(Map$Genotype,Map$Phosphate,sep = "_") %>%
  factor %>% relevel(ref = c("Col_0_low"))

Dat <- create_dataset(Tab = Tab ,Map = Map)



Dat_rnaseq_halle <- list(Dat_rnaseq_halle = Dat, 
                         df_length_genes = df_length_genes)

saveRDS(file = "../cleandata/dat_rnaseq_hallepi.RDS",object = Dat_rnaseq_halle)
