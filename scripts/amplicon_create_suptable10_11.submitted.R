library(ohchibi)



#Create supplemental tables for  the GLM in the syncom
#Read enricment results
res_frac <- read.table(file = "../cleandata/df_dds_res_syncom_fraction_specific.tsv",header = T,sep = "\t")

#Create modifications for the supplemental table
res_frac$Id <- res_frac$Id %>% gsub(pattern = "Sequence_",replacement = "USeq")
res_frac <- res_frac[,c(8,1:7,9:10)]
res_frac <- res_frac[,1:8]

res_int <- read.table(file = "../cleandata/df_dds_res_syncom_phosphateint.tsv",header = T,sep = "\t")
res_int$Id <- res_int$Id %>% gsub(pattern = "Sequence_",replacement = "USeq")
res_int <- res_int[,c(8,1:7,9:10)]
res_int <- res_int[,c(1:8,10)]
colnames(res_int)[9] <- "Fraction"


#Write both tables
write.table(x = res_frac,file = "../cleandata/sup_table_10.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)
write.table(x = res_int,file = "../cleandata/sup_table_11.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)
