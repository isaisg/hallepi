library(ohchibi)
library(emmeans)
source('plotting_parameters_hallepi.R')
#Read the data
df <- read.table("../rawdata/PiDOcombined.csv",
                 header = T,sep = ",",quote = "",comment.char = "")

#Clean the headers
colnames(df) <- c("Tube","Syncom","Treatment","PiConcentration","Rep")

df <- df %>% subset(Rep != "a") %>% droplevels
df$Treatment <- df$Treatment %>% gsub(pattern = ">",replacement = "-") %>%
  gsub(pattern = "--",replacement = "_")
#df <- df %>% subset(Treatment != "50_1000") %>% droplevels
df$Treatment <- df$Treatment %>% factor(levels = c("50","1000","50_1000"))

df$Syncom <- df$Syncom %>% gsub(pattern = "\\(|\\)",replacement = "") %>%
  factor(levels = c("NB","full","-Burk","-Ralst"))
#Take into account missing normalization when creating the metadata

Res_em <- NULL
for(treat in levels(df$Treatment)) {
  temp <- df %>% subset(Treatment == treat) %>% droplevels
  temp$logPi <- log(temp$PiConcentration)
  m1 <- aov(data = temp,formula = logPi ~ Syncom + Rep)
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

#Remove outliers for visualization
df_sub <- df %>% subset(PiConcentration < 0.03) %>% droplevels 
Res_em$Treatment <- Res_em$Treatment %>% factor(levels = c("50","1000","50_1000"))

#Create boxplot
p <- chibi.boxplot(Map = df_sub,x_val = "Syncom",y_val = "PiConcentration",
                   style = "open",
                   facet_formula = "Treatment",size_median = 2,
                   size_point = 0,stroke_point = 0,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "red",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "Phosphate") + 
  xlab(label = "Syncom") +
  geom_sina(alpha = 0.15, color = "#414141",
            stroke = 0,size = 4) +
  geom_text(data = Res_em,
            aes(x = Syncom,y = 0.03,label = Letters),
            inherit.aes = F,size = 7.5,family ="Arial") +
  coord_cartesian(expand = TRUE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())   +
  geom_vline(xintercept = c(1.5,2.5,3.5),size = 0.6,color = "#D9D9D9")
oh.save.pdf(p = p,outname = "pi_dropout_syncoms.pdf",outdir = "../figures/")
