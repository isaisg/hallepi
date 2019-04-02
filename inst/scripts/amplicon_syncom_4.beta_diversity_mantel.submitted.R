library(ohchibi)
library(harrietr)

size_legend_text <- 45
size_axis_title <- 35
size_legend_title <- 55
legend_proportion_size <- 4


Dat_ori <- readRDS("../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS")
Dat_rar <- Dat_ori$RelativeAbundance

colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
#colores_phosphate <- c("#fff7f3","#fde0dd","#fa9fb5","#dd3497","#7a0177","#49006a")


Dat_rar$Map$typebyTissue <- Dat_rar$Map$typebyTissue %>% as.character %>%
  gsub(pattern = "Agar_NP",replacement = "AgarNoPlant") %>%
  gsub(pattern = "Agar_SC",replacement = "AgarPlant") %>%
  gsub(pattern = "Root_SC",replacement = "Root") %>%
  gsub(pattern = "Shoot_SC",replacement = "Shoot") %>%
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))

#Compute distances
#Define Bray Curtis as the dissimilarity
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

Dat_rar_root <- Dat_rar %>% subset.Dataset(typebyTissue == "Root",drop = T,clean = T)
Dat_rar_shoot <- Dat_rar %>% subset.Dataset(typebyTissue == "Shoot",drop = T,clean = T)
Dat_rar_agar <- Dat_rar %>% subset.Dataset(typebyTissue == "AgarPlant",drop = T,clean = T)

#Calculate the distances 
dist_root <- distfun(x = Dat_rar_root$Tab %>% t)
dist_shoot <- distfun(x = Dat_rar_shoot$Tab %>% t)
dist_agar <- distfun(x = Dat_rar_agar$Tab %>% t)

#Now here we want to create and average
melted_root <- dist_root %>% as.matrix %>% melt
melted_shoot <- dist_shoot %>% as.matrix %>% melt
melted_agar <- dist_agar %>% as.matrix %>% melt

#Associate the metadata to each melted object
Mapa <- Dat_rar$Map
melted_root$D1 <- melted_root$Var1 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character
melted_root$D2 <- melted_root$Var2 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character


melted_shoot$D1 <- melted_shoot$Var1 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character
melted_shoot$D2 <- melted_shoot$Var2 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character


melted_agar$D1 <- melted_agar$Var1 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character
melted_agar$D2 <- melted_agar$Var2 %>% match(Mapa$ID_Matrix) %>% Mapa$TTC[.] %>% as.character

#Clean names
melted_root$D1 <- melted_root$D1 %>% gsub(pattern = "Root_SC_",replacement = "")
melted_root$D2 <- melted_root$D2 %>% gsub(pattern = "Root_SC_",replacement = "")

melted_shoot$D1 <- melted_shoot$D1 %>% gsub(pattern = "Shoot_SC_",replacement = "")
melted_shoot$D2 <- melted_shoot$D2 %>% gsub(pattern = "Shoot_SC_",replacement = "")

melted_agar$D1 <- melted_agar$D1 %>% gsub(pattern = "Agar_SC_",replacement = "")
melted_agar$D2 <- melted_agar$D2 %>% gsub(pattern = "Agar_SC_",replacement = "")



#Create summary per combination 
dist_root <- acast(data = melted_root,formula = D1~D2,fun.aggregate = mean,value.var = "value") %>% as.dist

dist_shoot <- acast(data = melted_shoot,formula = D1~D2,fun.aggregate = mean,value.var = "value") %>% as.dist

dist_agar <- acast(data = melted_agar,formula = D1~D2,fun.aggregate = mean,value.var = "value") %>% as.dist



m_root_shoot <- mantel(xdis = dist_root,ydis = dist_shoot,
                       permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_root_shoot_syncom.doc",
               append = F,print(m_root_shoot))


df1 <- melt_dist(dist = dist_root %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_root_shoot$statistic %>% round(3) ,p = m_root_shoot$signif %>% format.pval(digits = 2))
dfpval$label <- paste0(dfpval$r,"\n",dfpval$p)
melted <- data.frame(Root = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_root)
p <- ggplot(data = melted,aes(Root,Shoot)) + geom_point(shape = 21,,size = 10) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.25,0.2,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_bacteria_root_shoot_syncom.pdf",outdir = "../figures/",width = 10,height = 10)



m_root_soil <- mantel(xdis = dist_root,ydis = dist_agar,
                      permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_root_agar_syncom.doc",
               append = F,print(m_root_soil))

df1 <- melt_dist(dist = dist_agar %>% as.matrix) 

df2 <- melt_dist(dist = dist_root %>% as.matrix) 

dfpval <- data.frame(r = m_root_soil$statistic %>% round(3) ,p = m_root_soil$signif %>% format.pval(digits = 2))
dfpval$label <- paste0(dfpval$r,"\n",dfpval$p)
melted <- data.frame(Agar = df1$dist,Root = df2$dist)
mlm <- lm(dist_root~dist_agar)
p <- ggplot(data = melted,aes(Agar,Root)) + geom_point(shape = 21,,size = 10) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.3,0.2,label =label),size = 5)

oh.save.pdf(p = p,outname = "mantel_bacteria_agar_root_syncom.pdf",outdir = "../figures/",width = 10,height = 10)





m_shoot_soil <- mantel(xdis = dist_shoot,ydis = dist_agar,
                       permutations = 10000)
capture.output(file = "../figures/supfigure2_mantel_shoot_agar_syncom.doc",
               append = F,print(m_shoot_soil))

df1 <- melt_dist(dist = dist_agar %>% as.matrix) 

df2 <- melt_dist(dist = dist_shoot %>% as.matrix) 
dfpval <- data.frame(r = m_shoot_soil$statistic %>% round(3) ,p = m_shoot_soil$signif %>% format.pval(digits = 2))
dfpval$label <- paste0(dfpval$r,"\n",dfpval$p)
melted <- data.frame(Agar = df1$dist,Shoot = df2$dist)
mlm <- lm(dist_shoot~dist_agar)
p <- ggplot(data = melted,aes(Agar,Shoot)) + geom_point(shape = 21,,size = 10) + 
  theme_ohchibi() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  geom_abline(intercept = mlm$coefficients[1],slope = mlm$coefficients[2],color = "red",size = 1) +
  geom_text(dfpval,mapping = aes(0.25,0.25,label =label),size = 5)
oh.save.pdf(p = p,outname = "mantel_bacteria_agar_shoot_syncom.pdf",outdir = "../figures/",width = 10,height = 10)
