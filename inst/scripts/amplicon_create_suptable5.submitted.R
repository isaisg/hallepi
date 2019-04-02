library(ohchibi)

res <- read.table(file = "../cleandata/df_dds_res_amplicon_fungi_otus_fraction_asvlevel.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")


#Write sup table 4
suptab <- res
suptab$Contrast <- suptab$Contrast %>% factor(levels = c("Root-Soil","Shoot-Soil","Root-Shoot"))
suptab$Num <- suptab$OTU_Id %>% gsub(pattern = "OTU_",replacement = "") %>% as.numeric
suptab <- with(suptab,order(Num,Contrast)) %>% suptab[.,]
suptab <- suptab[,-15]

#Load relative abundance information
Dat <- readRDS("../cleandata/dat_hallepi_amplicon_fungi_otus_soil.2019.RDS")
melted <- Dat$RelativeAbundance$Tab %>% melt
colnames(melted) <- c("OTU_Id","OTU_Header","RA")
melted <- merge(melted,Dat$RelativeAbundance$Map, by = "OTU_Header")

#Create cumulative relative abundance
df_ab <- dcast(data = melted,formula = OTU_Id~Fraction,
               fun.aggregate = mean,value.var = "RA")


suptab <- merge(suptab,df_ab, by = "OTU_Id")
colnames(suptab)[15] <- "CumulativeRA_BulkSoil"
colnames(suptab)[16] <- "CumulativeRA_Soil"
colnames(suptab)[17] <- "CumulativeRA_Root"
colnames(suptab)[18] <- "CumulativeRA_Shoot"

#Clean the taxonomy columns
suptab$Kingdom <- suptab$Kingdom %>% gsub(pattern = "k__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")
suptab$Phylum <- suptab$Phylum %>% gsub(pattern = "p__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Class <- suptab$Class %>% gsub(pattern = "c__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")
suptab$Order <- suptab$Order %>% gsub(pattern = "o__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Family <- suptab$Family %>% gsub(pattern = "f__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")
suptab$Genus <- suptab$Genus %>% gsub(pattern = "g__",replacement = "")%>%
  gsub(pattern = " ",replacement = "")


write.table(x = suptab,file = "../cleandata/sup_table_5.csv",append = F,
            quote = F,sep = ",",row.names = F,col.names = T)
