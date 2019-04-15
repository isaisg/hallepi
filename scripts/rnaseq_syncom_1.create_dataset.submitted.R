library(ohchibi)

#Set random seed
set.seed(130816)



#Read Tab and Map files
Tab <- read.table(file = "../rawdata/rnaseq_tab_pi_syncom.tsv",
                  header = T,row.names = 1,sep = "\t",quote = "",comment.char = "",check.names = F)
Map <- read.table(file = "../rawdata/rnaseq_map_pi_syncom.tsv",
                  header = T,row.names = 1,sep = "\t",quote = "",comment.char = "")
df_length_genes <- read.table(file = "../rawdata/rnaseq_df_length_genes.tsv",
                              header = T,sep = "\t",quote = "",comment.char = "")

Map$Syncom <- Map$Syncom %>% 
  gsub(pattern = "full",replacement = "Full") %>%
  factor(levels = c("NB","Full"))

Map$Conditions <- Map$Conditions %>% 
  gsub(pattern = "control",replacement = "1000Pi") %>%
  gsub(pattern = "-P",replacement = "50Pi") %>%
  gsub(pattern = "\\(|\\)",replacement = "") %>%
  factor(levels = c("50Pi","1000Pi"))

Map$group <- paste(Map$Syncom,Map$Conditions,sep = "_") %>%
  factor %>% relevel(ref = "NB_1000Pi")

Dat <- create_dataset(Tab = Tab ,Map = Map)



Dat_rnaseq_syncom <- list(Dat_rnaseq_syncom = Dat, 
                          df_length_genes = df_length_genes)
saveRDS(object = Dat_rnaseq_syncom,file = "../cleandata/dat_rnaseq_syncom.RDS")
