  library(ohchibi)
  library(DESeq2)
  #Set random seed
  set.seed(130816)
  
  
  
  Dat_raw <- readRDS(file = "../cleandata/dat_rnaseq_syncom.RDS")
  Dat <- Dat_raw$Dat_rnaseq_syncom
  
  ### Create DEseq object ###
  dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,colData = Dat$Map,
                              design = ~ Experiment + group)
  
  vsd <- vst(object = dds,blind = FALSE)
  
  #Melt object and z-scale it
  melted_vsd <- vsd %>% assay %>% t  %>% scale %>% t %>% melt
  
  colnames(melted_vsd) <- c("Gene","Sample_Id","value")
  
  melted_vsd <- merge(melted_vsd,Dat$Map, by = "Sample_Id")
  
  #Low the psr genes
  regulon_psr <- read.table("../rawdata/regulon_psr.csv",header = F) %$%
  V1 %>% as.character
  
  melted_vsd <- which(melted_vsd$Gene %in% regulon_psr) %>%
  melted_vsd[.,] %>% droplevels
  
  
  #Calculate the mean per gene so we dont overinflate the points
  melted_psr <- dcast(data = melted_vsd,formula = Gene ~ group,fun.aggregate = mean,value.var = "value") %>%
  melt
  colnames(melted_psr)[2] <- "group"
  
  melted_psr_end <- melted_vsd[,c(7:8,10)] %>% unique %>% merge(melted_psr, . , by = "group",all.x = TRUE)
  
  #Read model results
  lista <- readRDS("../cleandata/list_psr_regulon_significance_insyncomrnaseq.RDS")
  axenic <- lista$axenic
  full <- lista$full
  
  df_nb <- melted_psr_end %>% subset(Syncom == "NB") %>% droplevels
  df_full <- melted_psr_end %>% subset(Syncom == "Full") %>% droplevels
  
  df_nb$Significance <- rep("NoSignificant",nrow(df_nb))
  df_nb$Significance[which(df_nb$Gene %in% axenic)] <- "Significant"
  
  df_full$Significance <- rep("NoSignificant",nrow(df_full))
  df_full$Significance[which(df_full$Gene %in% full)] <- "Significant"
  
  melted_psr_end <- rbind(df_nb,df_full)
  melted_psr_end$Significance <- melted_psr_end$Significance %>% factor
  
  melted_psr_end <- na.omit(melted_psr_end)
  
  
  p <- ggplot(data = melted_psr_end,aes(Conditions,value)) +
  facet_grid(.~Syncom) +
  geom_boxplot(outlier.shape = NA,outlier.size = NA) + 
  geom_sina(aes(color = Significance,group = Conditions),size = 7,alpha = 0.4,fill = NA) +
  #geom_point(position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1),
  #           aes(color = Significance),shape = 21,size = 7) +
  theme_ohchibi() +
  geom_vline(xintercept = c(1.5),size = 0.5, color = "#D9D9D9") +
  ylab("z-score") + theme(legend.position = "none",
                          axis.title.x = element_blank()) +
  scale_size(range=c(3,6))  +
  scale_color_manual(values = c("black","#E55000"))
  
  
  write.table(x = melted_psr,file = "../data_figures/data_S5B.csv",
              append = F,quote = F,sep = ",",row.names = F,col.names = T)
  
  
  
  
  oh.save.pdf(p = p,outname = "figure1_extra_rnaseq_hallepi_psr_regulon_syncom.pdf",outdir = "../figures")
  
