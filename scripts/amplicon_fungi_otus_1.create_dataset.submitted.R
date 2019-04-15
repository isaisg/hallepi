library(ohchibi)
set.seed(130816)

#Create cleandata and figures folders
dir.create(path = "../cleandata")
dir.create(path = "../figures")


#Read otutable
Tab <- read.table(file = "../rawdata/otu_table_fungi.txt",header = T,sep = "\t",row.names = 1,quote = "",check.names = F)

#Read metadata
Map <- read.table("../rawdata/Soil_Pi_gradient_metadata.2019.csv.withotuheader.tsv",header=T,sep="\t")
Map$DADA2_Header <- factor(paste0(Map$Plate,"_",Map$Well))
Map$Genotype <- Map$Genotype %>% 
  gsub(pattern = "No_plant",replacement = "No_Plant") %>%
  gsub(pattern = "col-0",replacement = "Col-0") %>%
  gsub(pattern = "phf",replacement = "phf1") %>%
  gsub(pattern = "phr/phl",replacement = "phr1/phl1")

Map$Phosphate <- Map$Soil %>% as.character %>%
  gsub(pattern = "-.*",replacement = "")
Map$Phosphate[which(is.na(Map$Phosphate))] <- "Blank"


Map$fraction <- Map$fraction %>% as.character
Map$Genotype <- Map$Genotype %>% as.character
Map$plot <- Map$plot %>% as.character
Map$Soil <- Map$Soil %>% as.character

Map$Soil[which(Map$Sample == "B")] <- "Blank"
Map$plot[which(Map$Sample == "B")] <- "Blank"
Map$fraction[which(Map$Sample == "B")] <- "Blank"
Map$Genotype[which(Map$Sample == "B")] <- "Blank"
Map$fraction[which(Map$Genotype == "No_Plant")] <- "BulkSoil"

Map$Genotype <- Map$Genotype %>%
  factor(levels = c("Blank","No_Plant","Col-0","phf1","phr1/phl1"))

Map$Phosphate <- Map$Phosphate %>% 
  factor(levels = c("Blank","low","medium","high","low+Pi"))

Map$fraction <- Map$fraction %>%
  factor(levels = c("Blank","BulkSoil","Soil","Root","Shoot"))


#Rename the  rownames in the Metadata
rownames(Map) <- Map$OTU_Header
colnames(Map)[8] <- "Plot"
colnames(Map)[9] <- "Fraction"

#Create the pair variable
Map$Pot <- Map$Sample %>% as.character %>%
  gsub(pattern = "L|R|S",replacement = "") %>%
  paste0("Pot",.)

colnames(Tab) <- colnames(Tab) %>% toupper
Map$Fungi_Tab <- paste0(Map$Plate,Map$Well) %>% gsub(pattern = "Plate",replacement = "")

Tab <- Tab[,match(Map$Fungi_Tab,colnames(Tab))]
colnames(Tab) <- Map$OTU_Header

#Read the taxonomy
Tax <- read.table(file = "../rawdata/taxonomy_fungi_otus.txt",header = F,sep = "\t")
Tax$qiime <- paste("Root; k__",Tax$V3,"; p__",Tax$V6,"; c__",Tax$V9,"; o__",
                   Tax$V12,"; f__",Tax$V15,"; g__",Tax$V18,sep = "")
myTax <- data.frame(ID=Tax$V1,Taxonomy=Tax$qiime)
rownames(myTax) <- Tax$V1

#Reorder the structures to create the Dat object
rownames(Map) <- Map$OTU_Header
myTax<-myTax[match(rownames(Tab),rownames(myTax)),]


Dat<-create_dataset(Tab = Tab,Map = Map,Tax = myTax)

#Compute depth
Dat$Map$Depth <- colSums(Dat$Tab)

#Add Plate label to the Plate column in the Metadata
#Add Plate label to the Plate column in the Metadata
Dat$Map$Plate <- Dat$Map$Plate %>% factor
Dat$Map$Plot <- factor(paste("Plot",Dat$Map$Plot,sep=""),
                       levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
                                ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13",
                                "Plot14","Plot15","Plot16","Plot17","Plot18"))
Dat$Map$LDepth <- log(Dat$Map$Depth)

##Plot usable reads
p <- chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",
                   style = "open",size_axis_text.x = 30,size_axis_text.y = 30,
                   size_point = 10,size_median = 5) + scale_y_log10() + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),size = 0.3, color = "#D9D9D9") +
  ylab(label = "log Depth")
#oh.ggsave.svg(outdir = "../figures/",ggobject = p,outname = "amplicon_fungi.usable_reads_fraction.svg")


#Subset Root Shoot and Soil and color by genotype
Dat$Map %>% subset(Fraction == "Root" | Fraction == "Shoot" |
                     Fraction == "Soil") %>% droplevels %>% 
  chibi.boxplot(Map = .,x_val = "Genotype",y_val = "LDepth",
                col_val = "Phosphate",facet_formula = "Fraction",
                style = "open",size_axis_text.x = 30,size_axis_text.y = 30,
                size_point = 10,size_median = 5,median_color = "#414141") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),size = 0.3, color = "#D9D9D9") +
  ylab(label = "log10 Depth")


#Subset the Blank Samples
Dat<-subset.Dataset(x = Dat,subset = Fraction!="Blank",drop = T,clean = T)

#Remove sapmles with less than 3000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000 ])
Dat_raw <- clean(Dat_raw)
chibi.boxplot(Map = Dat_raw$Map,x_val = "Genotype",y_val = "LDepth",
              col_val = "Phosphate",facet_formula = "Fraction",
              style = "open",size_axis_text.x = 30,size_axis_text.y = 30,
              size_point = 10,size_median = 5,median_color = "#414141") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),size = 0.3, color = "#D9D9D9")


### Rarefaction
set.seed(3111)
Dat_rar <- rarefaction(x = Dat_raw,sample =3000)
Dat_rar <- clean(Dat_rar)


###Relative abundance 
Dat_rab <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000 ])
Dat_rab<-clean(Dat_rab)
Dat_rab$Tab<-scale(x = Dat_rab$Tab,center = F,scale = colSums(Dat_rab$Tab))

##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_amplicon<-list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_rab)
filename <- "../cleandata/dat_hallepi_amplicon_fungi_otus_soil.2019.RDS"
saveRDS(object = Dat_amplicon,file = filename)

#Create a dataframe of taxonomy
Tax <- Dat_raw$Tax
df_Taxonomy <- Tax$Taxonomy %>% as.character %>% 
  strsplit(split = "\\;") %>% unlist %>%
  matrix(ncol = 7,byrow = T) %>% 
  as.data.frame
colnames(df_Taxonomy) <- c("Root","Kingdom","Phylum","Class","Order","Family","Genus")
df_Taxonomy$Id <- Tax$ID
df_Taxonomy <- df_Taxonomy[,-1]
write.table(x = df_Taxonomy,file = "../rawdata/df_taxonomy_fungi_otus.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
