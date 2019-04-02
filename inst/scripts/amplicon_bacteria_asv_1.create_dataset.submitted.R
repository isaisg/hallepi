library(ohchibi)

#Create cleandata and figures folders
dir.create(path = "../cleandata")
dir.create(path = "../figures")


#Read files
Tab <- read.table(file = "../rawdata/dada2_structure_maxee0_nopool.tsv",
                  header = T,sep = "\t",check.names = F,row.names = 1)
rownames(Tab) <- Tab %>% rownames %>% 
  gsub(pattern = ".fastq.filt.fastq.gz",replacement = "")
Tab <- t(Tab)

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
rownames(Map) <- Map$DADA2_Header
colnames(Map)[8] <- "Plot"
colnames(Map)[9] <- "Fraction"

#Create the pair variable
Map$Pot <- Map$Sample %>% as.character %>%
  gsub(pattern = "L|R|S",replacement = "") %>%
  paste0("Pot",.)

# Match tab and map
Map <- Map[ (row.names(Map) %in% colnames(Tab)), ]
Tab <- Tab[ ,colnames(Tab) %in% row.names(Map)]
Map <- Map[colnames(Tab),]

#Handle the taxtonomic object derived from dada2
Tax<-read.table(file = "../rawdata/dada2_structure_maxee0_nopool.taxaclassification.tsv",
                header = T,sep = "\t",quote = "",comment.char = "")
Tax$Taxonomy <- paste0("Root; k__",Tax$Kingdom,"; p__",
                       Tax$Phylum,"; c__",Tax$Class,"; o__",Tax$Order,
                       "; f__",Tax$Family,"; g__",Tax$Genus)
mTax <- data.frame(ID = Tax$ASV_Id,Taxonomy = Tax$Taxonomy)
rownames(mTax) <- mTax$ID
mTax <- match(rownames(Tab),rownames(mTax)) %>%
  mTax[.,]

#Create dataset object 
Dat <- create_dataset(Tab = Tab,Map = Map,Tax = mTax)

#Create vector of contaminants to remove from the ASV mapping
contam_otus <- c(grep(pattern = "Chloroplast", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "Mitochondria", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__NA", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "k__Eukaryota", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 #Remove the  ones without phylum assignment
                 grep(pattern = "p__NA", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE)
                 
)
contam_otus <- row.names(Dat$Tax)[contam_otus]

# Filter Dataset
Dat_filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat_filter <- clean(Dat = Dat_filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]


Dat <- Dat_filter

#Compute depth
Dat$Map$Depth <- colSums(Dat$Tab)

#Add Plate label to the Plate column in the Metadata
Dat$Map$Plate <- Dat$Map$Plate %>% factor
Dat$Map$Plot <- factor(paste("Plot",Dat$Map$Plot,sep=""),
                       levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
                                ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13",
                                "Plot14","Plot15","Plot16","Plot17","Plot18"))
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
  ylab(label = "log Depth")


#Subset the Blank Samples
Dat <- subset.Dataset(x = Dat,subset = Fraction!="Blank",
                      drop = T,clean = T)


#Remove samples with less than 3000 reads
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
filename <- "../cleandata/dat_hallepi_amplicon_bacteria_ASV_soil.2019.RDS"
saveRDS(object = Dat_amplicon,file = filename)
