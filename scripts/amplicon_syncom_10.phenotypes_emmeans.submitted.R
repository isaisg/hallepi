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

#Loop over each concentratrion
Res_em <- NULL
Res_pval <- NULL
for(cn in levels(df$Treatment)){
  temp <- df %>% subset(Treatment == cn) %>% droplevels
  m1 <- aov(data = temp,formula = TotalShootLength_cm~Bacteria+BioRep)
  m1_em <- m1 %>% emmeans(specs = "Bacteria")
  temp_pval <- pairs(m1_em,reverse = T) %>% as.data.frame 
  temp_em <- m1_em %>% as.data.frame
  temp_pval$Treatment <- rep(cn,nrow(temp_pval))
  temp_em$Treatment <- rep(cn,nrow(temp_em))
  Res_em <- rbind(Res_em,temp_em)
  Res_pval <- rbind(Res_pval,temp_pval)
}

### Create plot 
paleta <- c("#D59F0F","#002B7A")

Res_em$Treatment <- Res_em$Treatment %>% factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))
Res_pval$Treatment <- Res_pval$Treatment  %>% factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi"))
Res_pval$p.adj <- Res_pval$p.value %>% p.adjust(method = "BY")
Res_pval$Fig <- rep("NS",nrow(Res_pval))
Res_pval$Fig[Res_pval$p.adj <= 0.05] <- "*"
Res_pval$Fig[Res_pval$p.adj <= 0.01] <- "**"
Res_pval$Fig[Res_pval$p.adj <= 0.001] <- "***"


p <- ggplot(data = Res_em,mapping = aes(Treatment,emmean,fill = Bacteria)) +
  geom_pointrange(aes(y = emmean, ymin = lower.CL, ymax= upper.CL,color = Bacteria), 
                  ,size = 3) +
  geom_point(shape = 21,size = 3,aes(color = Bacteria))  +
  geom_line(aes(y = emmean,color = Bacteria,group = Bacteria),size =5) +
  theme_ohchibi() + 
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  scale_color_manual(values = paleta)+ scale_fill_manual(values = paleta) +
  geom_text(data = Res_pval,mapping = aes(x = Treatment,y = 20,label = Fig),inherit.aes = F) 
oh.save.pdf(p = p,outname = "figure6_syncom_shootphenotypes_anova_emmean.pdf",outdir = "../figures")


