library(ohchibi)

df <- read.table(file = "../cleandata/sup_table_2.csv",header = T,sep = ",")

df_desc <- read.table(file = "../cleandata/mgene_description_20131231.txt",
                      header = F,sep = "\t",fill = NA,quote = "",comment.char = "")
df_desc$Gene <- df_desc$V1 %>% gsub(pattern = "\\..*",replacement = "")

Description <- match(df$Gene,df_desc$Gene) %>%
  df_desc$V2[.] %>% as.character
df$Description <- Description

df$PSR_Regulon <- df$PSR_Regulon %>% factor %>% relevel(ref = "Yes")
df <- with(df,order(PSR_Regulon)) %>% 
  df[.,]

df$Significance <- df$Significance %>% gsub(pattern = "NoSignificant",replacement = "NotSignificant")

write.table(x = df,file = "../cleandata/sup_table_2.annotated.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)
