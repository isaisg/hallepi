library(ohchibi)
library(emmeans)

set.seed(130816)

source(file = "plotting_parameters_hallepi.R")


df <- read.table(file = "../rawdata/hallepi_syncom_dropout_shootdata.tsv",
                 header = T,sep = "\t")

df$Treatment <- df$Treatment %>% factor(levels = c("50uM","1000uM","50uM-1000uM"))
df$Syncom <- df$Syncom %>% 
  factor(levels = c("NB","Full","-Burk","-Ralst"))


#Do testing inside each Pi treatment
Res_em <- NULL
for(treat in levels(df$Treatment)) {
  temp <- df %>% subset(Treatment == treat) %>% droplevels
  m1 <- aov(data = temp,formula = TotalLength_cm ~ Syncom + Exp)
  m1_em <- emmeans(object = m1,specs = "Syncom") %>% CLD
  m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
  colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
  m1_em$Treatment <- rep(treat,nrow(m1_em))
  Res_em <- rbind(Res_em,m1_em)
  
}


#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>% 
  strsplit(split = "") %>% 
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop 
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
    unlist %>% toupper  %>% paste0(collapse = "")
}

Res_em$Treatment <- Res_em$Treatment %>% 
  factor(levels = c("50uM","1000uM","50uM-1000uM"))


#Create boxplot
p <- chibi.boxplot(Map = df,x_val = "Syncom",y_val = "TotalLength_cm",
                   style = "open",
                   facet_formula = "Treatment",size_median = 2,
                   size_point = 0,stroke_point = 0,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "red",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "TotSurfArea_cm2") + 
  xlab(label = "Syncom") +
  geom_sina(alpha = 0.15, color = "#414141",
            stroke = 0,size = 4) +
  geom_text(data = Res_em,
            aes(x = Syncom,y = 0.03,label = Letters),
            inherit.aes = F,size = 7.5,family ="Arial") +
  coord_cartesian(expand = TRUE) +
  theme(
    #axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )   
oh.save.pdf(p = p,outname = "shoot_dropout_syncoms.pdf",outdir = "../figures/")

#Write dataset
write.table(x = df,file = "../data_figures/data_S8.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)


