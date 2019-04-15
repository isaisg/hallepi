library(ohchibi)

res <- read.table(file = "../cleandata/df_dds_res_amplicon_bacteria_asv_fraction_asvlevel.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")


#Write sup table 4
suptab <- res
suptab$Contrast <- suptab$Contrast %>% factor(levels = c("Root-Soil","Shoot-Soil","Root-Shoot"))
suptab$Num <- suptab$ASV_Id %>% gsub(pattern = "ASV",replacement = "") %>% as.numeric
suptab <- with(suptab,order(Num,Contrast)) %>% suptab[.,]
suptab <- suptab[,-15]

#Load relative abundance information
Dat <- readRDS("../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS")
melted <- Dat$RelativeAbundance$Tab %>% melt
colnames(melted) <- c("ASV_Id","DADA2_Header","RA")
melted <- merge(melted,Dat$RelativeAbundance$Map, by = "DADA2_Header")

#Create cumulative relative abundance
df_ab <- dcast(data = melted,formula = ASV_Id~Fraction,fun.aggregate = mean,value.var = "RA")


suptab <- merge(suptab,df_ab, by = "ASV_Id")
colnames(suptab)[15] <- "CumulativeRA_BulkSoil"
colnames(suptab)[16] <- "CumulativeRA_Soil"
colnames(suptab)[17] <- "CumulativeRA_Root"
colnames(suptab)[18] <- "CumulativeRA_Shoot"

write.table(x = suptab,file = "../cleandata/sup_table_4.csv",append = F,
            quote = F,sep = ",",row.names = F,col.names = T)
