library(ohchibi)

#Set random seed
set.seed(130816)




#Read both alpha diversities
df_bac <- read.table(file = "../cleandata/sup_table_3_bacteria.csv",header = T,
                     sep = "\t")

df_fungi <- read.table(file = "../cleandata/sup_table_3_fungi.csv",header = T,
                     sep = "\t")


colnames(df_bac)
df_bac <- df_bac[,-c(14,15)]
df_fungi <- df_fungi[,c(1:13,17,18)]

temp1 <- df_bac[,-c(14:15)]
temp2 <- df_fungi[,-c(14:15)]

df_all <- rbind(temp1,temp2) %>% unique

#Populate df_all with bacteria
df_all$Shannon_Bacteria <- rep(NA,nrow(df_all))
df_all$Shannon_Bacteria[match(df_bac$DADA2_Header,df_all$DADA2_Header)] <- df_bac$Shannon

df_all$Richness_Bacteria <- rep(NA,nrow(df_all))
df_all$Richness_Bacteria[match(df_bac$DADA2_Header,df_all$DADA2_Header)] <- df_bac$Richness

#Populate df_all with fungi
df_all$Shannon_Fungi <- rep(NA,nrow(df_all))
df_all$Shannon_Fungi[match(df_fungi$OTU_Header,df_all$OTU_Header)] <- df_fungi$Shannon


df_all$Richness_Fungi <- rep(NA,nrow(df_all))
df_all$Richness_Fungi[match(df_fungi$OTU_Header,df_all$OTU_Header)] <- df_fungi$Richness


write.table(x = df_all,file = "../cleandata/sup_table_3.csv",append = F,
            quote = F,sep = ",",row.names = F,col.names = T)
