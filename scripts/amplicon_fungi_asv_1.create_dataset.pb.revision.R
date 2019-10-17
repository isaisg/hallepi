library(ohchibi)
setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts/')
set.seed(130816)

#Create cleandata and figures folders
dir.create(path = "../cleandata")
dir.create(path = "../figures")


#Read otutable
Tab <- read.table(file = "../rawdata/hallepi_wildsoil_fungi_dada2_res_structure.tsv",
                  header = T,sep = "\t",row.names = 1,quote = "",check.names = F)

#Read metadata
Map <- read.table("../rawdata/hallepi_metadata_soil_pi_gradient.tsv",header=T,sep="\t")


#Adjust the factors
Map$Genotype <- Map$Genotype %>%
  factor(levels = c("Blank","No_Plant","Col-0","phf1","phr1/phl1"))

Map$Phosphate <- Map$Phosphate %>% 
  factor(levels = c("Blank","low","medium","high","low+Pi"))

Map$Fraction <- Map$Fraction %>%
  factor(levels = c("Blank","BulkSoil","Soil","Root","Shoot"))

Map$Plot <- Map$Plot %>% factor(
  levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
           ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13",
           "Plot14","Plot15","Plot16","Plot17","Plot18")
)

#Rename the  rownames in the Metadata
rownames(Map) <- Map$DADA2_Header

rownames(Tab) <- rownames(Tab) %>% toupper

Tab <- Tab %>% t
colnames(Tab) <- colnames(Tab) %>% 
  gsub(pattern = "^P",replacement = "Plate")

#Intersect usable
chosen <- intersect(Map$DADA2_Header,colnames(Tab))
Map <- which(Map$DADA2_Header %in% chosen) %>%
  Map[.,] %>% droplevels

Tab <- Tab[,match(Map$DADA2_Header,colnames(Tab))]
colnames(Tab) <- Map$DADA2_Header


#Handle the taxtonomic object derived from mothur
Tax<-read.table(file = "../rawdata/hallepi_wildsoil_fungi_seqsformothur.2019.wang.taxonomy",
                header = F,sep = "\t",quote = "",comment.char = "")
df_tax <- Tax$V2 %>% as.character %>% 
  gsub(pattern = "\\([0-9]+\\)",replacement = "") %>% 
  strsplit(x = .,split = "\\;") %>% unlist %>%
  matrix(data = .,nrow = nrow(Tax),ncol = 7,byrow = T) %>% 
  as.data.frame
colnames(df_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

Tax$Taxonomy <- paste0("Root; k__",df_tax$Kingdom,"; p__",
                       df_tax$Phylum,"; c__",df_tax$Class,"; o__",df_tax$Order,
                       "; f__",df_tax$Family,"; g__",df_tax$Genus)
mTax <- data.frame(ID = Tax$V1,Taxonomy = Tax$Taxonomy)
rownames(mTax) <- mTax$ID
mTax <- match(rownames(Tab),rownames(mTax)) %>%
  mTax[.,]

#Create dataset object 
Dat <- create_dataset(Tab = Tab,Map = Map,Tax = mTax)

#Select only the fungi ASVs
contam_otus <- Dat$Tax$Taxonomy %>% grep(pattern = "k__Fungi",invert = T) %>% 
  Dat$Tax$ID[.] %>% as.character

#contam_otus <- Dat$Tax$Taxonomy %>% grep(pattern = "p__Fungi_unclassified") %>%
#  Dat$Tax$ID[.] %>% as.character  %>%
#  c(contam_otus,.)


df_tax$Id <- Tax$V1
df_tax <- which(!(df_tax$Id %in% contam_otus)) %>%
  df_tax[.,] %>% droplevels

df_tax$Root <- "Root"
df_tax <- df_tax[,c(8:9,1:7)]

#### Quantify how many reads are being discarded by removing contaminants
melted <- Dat$Tab %>% melt
df_sums <- aggregate(value ~ Var1, data = melted,FUN = sum)
sum_removed <- which(df_sums$Var1 %in% contam_otus) %>%
  df_sums$value[.] %>% sum 

(sum_removed / (df_sums$value %>% sum))*100

# Filter Dataset
Dat_filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat_filter <- clean(Dat = Dat_filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]


Dat <- Dat_filter


#Compute depth
Dat$Map$Depth <- colSums(Dat$Tab)
Dat$Map$LDepth <- log(Dat$Map$Depth)



##Plot usable reads
chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",
              style = "open",size_axis_text.x = 30,size_axis_text.y = 30,
              size_point = 10,size_median = 5) + scale_y_log10() + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),size = 0.3, color = "#D9D9D9") +
  ylab(label = "log Depth")


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
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000])
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
Dat_amplicon<-list(RawCounts=Dat_raw,
                   Rarefied=Dat_rar,
                   RelativeAbundance=Dat_rab,
                   df_tax = df_tax)
filename <- "../cleandata/dat_hallepi_amplicon_fungi_asvs_soil.2019.RDS"
saveRDS(object = Dat_amplicon,file = filename)
rm(list=ls())
dev.off()
