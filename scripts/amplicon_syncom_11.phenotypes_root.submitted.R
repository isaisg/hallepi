library(ohchibi)
library(emmeans)
library(egg)
library(paletteer)



source('plotting_parameters_hallepi.R')

#Read the data
df <- read.table("../rawdata/phenotypes_syncom_pi.csv",
                 header = T,sep = ",",quote = "",comment.char = "")
df$Treatment <- df$Treatment %>% paste0(.,"Pi") %>%
  factor (levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))

df$BioRep <- df$BioRep %>% factor

paleta <- c("#D59F0F","#002B7A")


#Compute statistics for plot
df_res <- NULL
for(pi in levels(df$Treatment)){
  temp_df <- df %>% subset(Treatment == pi) %>%
    droplevels
  x <- temp_nb <- temp_df %>% subset(Bacteria == "NB") %$% TotalRootLength_cm
  y <- temp_nb <- temp_df %>% subset(Bacteria == "SC") %$% TotalRootLength_cm
  m1 <- t.test(x = x,y = y)
  temp <- data.frame(Treatment = pi,p.value = m1$p.value)
  df_res <- rbind(df_res,temp)
}
df_res$q.value <- p.adjust(p = df_res$p.value,method = "fdr")
df_res$q.value <- format.pval(pv = df_res$q.value,digits = 2)
### Create plot 

p_total <- chibi.boxplot(Map = df,x_val = "Bacteria",y_val = "TotalRootLength_cm",
              col_val = "Bacteria",facet_formula = "Treatment",
              mpalette = paleta,median_colored_as_points =T,size_point = 10) +
  theme(legend.position = "none",axis.title.x = element_blank()) +
  geom_text(data = df_res,mapping = aes(x = 1.5,y = 90,label = q.value),
            size = 5,inherit.aes = F)

#Compute statistics for plot
df_res <- NULL
for(pi in levels(df$Treatment)){
  temp_df <- df %>% subset(Treatment == pi) %>%
    droplevels
  x <- temp_nb <- temp_df %>% subset(Bacteria == "NB") %$% MainRootLength_cm
  y <- temp_nb <- temp_df %>% subset(Bacteria == "SC") %$% MainRootLength_cm
  m1 <- t.test(x = x,y = y)
  temp <- data.frame(Treatment = pi,p.value = m1$p.value)
  df_res <- rbind(df_res,temp)
}
df_res$q.value <- p.adjust(p = df_res$p.value,method = "fdr")
df_res$q.value <- format.pval(pv = df_res$q.value,digits = 2)
### Create plot 

p_main <- chibi.boxplot(Map = df,x_val = "Bacteria",y_val = "MainRootLength_cm",
                         col_val = "Bacteria",facet_formula = "Treatment",
                         mpalette = paleta,median_colored_as_points =T,size_point = 10) +
  theme(legend.position = "none",axis.title.x = element_blank()) +
  geom_text(data = df_res,mapping = aes(x = 1.5,y = 9,label = q.value),
            size = 5,inherit.aes = F)


oh.save.pdf(p = p_total,outname = "syncom_phenotypes_totalroot.pdf",
            outdir = "../figures")

oh.save.pdf(p = p_main,outname = "syncom_phenotypes_mainroot.pdf",
            outdir = "../figures")

