library(ohchibi)

#Create cleandata and figures folders
dir.create(path = "../cleandata")
dir.create(path = "../figures")


#####Create Dataset for Bacterial OTUs#####

#Create cleandata folders
outfolder<-paste0("../cleandata/bacteria_otus")
dir.create(outfolder)
#Read otutable
Dat <- read.am("../rawdata/otu_table_bacteria.txt",format = "qiime",taxonomy="taxonomy")

#Create contaminants to remove later
contam_otus <- c(grep(pattern = "chloroplast", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "mitochondri", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "oomycete", as.character(Dat$Tax$Taxonomy),ignore.case = TRUE),
                 which(Dat$Tax$Taxonomy == "Root; k__Bacteria"),
                 which(Dat$Tax$Taxonomy == "Root"))
contam_otus <- row.names(Dat$Tax)[contam_otus]

# Filter Dataset
Dat_filter <- remove_taxons(Dat = Dat, taxons = contam_otus)
Dat_filter <- clean(Dat = Dat_filter,verbose = TRUE)
contam_otus <- Dat$Tax[ contam_otus, ]

# Write resutls of filtering

filename <- paste(outfolder,"/tax_removed_97.txt",sep = "")
write.table(contam_otus,file = filename, col.names = NA, quote = FALSE, sep = "\t")

# Write filtered table
filename <- paste(outfolder,"/otu_table_filter_97.txt",sep = "")
write.qiime(Tab = Dat_filter$Tab, file = filename)
filename <- paste(outfolder,"/tax_filtered_97.txt",sep = "")
write.table(Dat_filter$Tax,file = filename)

# Clean
rm(Dat,Dat.filter,contam_otus,filename)

#Create folder of figures 
outfigures<-"../figures/bacteria_otus"
dir.create(path = outfigures)

#Reread the corresponding files to create final Dataset
filename <- paste(outfolder,"/otu_table_filter_97.txt",sep = "")
Tab <- read.am(filename,format = "qiime",taxonomy = FALSE)$Tab
filename <- paste(outfolder,"/tax_filtered_97.txt",sep = "")
Tax <- read.table(file = filename)
Map <- read.table("../rawdata/metadata_phosphate_gradient_hallesoil.tsv",header=T,sep="\t")

#Rename the  rownames in the Metadata
rownames(Map) <- Map$OTU_Header
# Match tab and map
Map <- Map[ (row.names(Map) %in% colnames(Tab)), ]
Tab <- Tab[ ,colnames(Tab) %in% row.names(Map)]
Map <- Map[colnames(Tab),]
Dat <- create_dataset(Tab = Tab,Map = Map, Tax = Tax)
Dat$Map$Depth <- colSums(Dat$Tab)

#Add Plate label to the Plate column in the Metadata
Dat$Map$Plate<-factor(paste("Plate",Dat$Map$Plate,sep=""))
Dat$Map$Plot<-factor(paste("Plot",Dat$Map$Plot,sep=""),levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
                                                                ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13","Plot14","Plot15","Plot16","Plot17","Plot18"))

##Plot usable reads
p<-chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",col_val = "Genotype")

oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot.svg")

Dat<-subset.Dataset(x = Dat,subset = Fraction!="Blank",drop = T,clean = T)
##Adjust the levels
Dat$Map$Condition<-factor(Dat$Map$Condition,levels=c("low","medium","high","low+Pi"))
Dat$Map$Genotype<-factor(Dat$Map$Genotype,levels=c("No_plant","Col-0","phf1","phr1/phl1"))
Dat$Map$Fraction<-factor(Dat$Map$Fraction,levels=c("BulkSoil","Soil","Root","Shoot"))

#Remove sapmles with less than 3000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000 ])
Dat_raw <- clean(Dat_raw)
p<-chibi.boxplot(Map = Dat_raw$Map,x_val = "Fraction",y_val = "Depth",col_val = "Genotype")
oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot_filtered3000.svg")

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
filename <- paste(outfolder,"dat_hallepi_bacteria_otus.RDS",sep = "/")
saveRDS(object = Dat_amplicon,file = filename)
rm(list=ls())

#####Create Dataset for Fungi OTUs#####
outfolder<-paste0("../cleandata/fungi_otus")
dir.create(outfolder)
outfigures<-"../figures/fungi_otus"
dir.create(outfigures)

#Read otutable
Tab <- read.table(file = "../rawdata/otu_table_fungi.txt",header = T,sep = "\t",row.names = 1,quote = "",check.names = F)
Map<-read.table(file = "../rawdata/metadata_phosphate_gradient_hallesoil.tsv",header = T,sep = "\t")

colnames(Tab)<-toupper(colnames(Tab))
Map$Fungi_Tab<-paste0(Map$Plate,Map$PosWell,Map$PosNum)

Tab<-Tab[,match(Map$Fungi_Tab,colnames(Tab))]
colnames(Tab)<-Map$OTU_Header

#Read the taxonomy
Tax<-read.table(file = "../rawdata/taxonomy_fungi_otus.txt",header = F,sep = "\t")
Tax$qiime<-paste("Root; k__",Tax$V3,"; p__",Tax$V6,"; c__",Tax$V9,"; o__",Tax$V12,"; f__",Tax$V15,"; g__",Tax$V18,sep = "")
myTax<-data.frame(ID=Tax$V1,Taxonomy=Tax$qiime)
rownames(myTax)<-Tax$V1

#Reorder the structures to create the Dat object
rownames(Map)<-Map$OTU_Header
myTax<-myTax[match(rownames(Tab),rownames(myTax)),]

Dat<-create_dataset(Tab = Tab,Map = Map,Tax = myTax)

#Compute depth
Dat$Map$Depth <- colSums(Dat$Tab)

#Add Plate label to the Plate column in the Metadata
Dat$Map$Plate<-factor(paste("Plate",Dat$Map$Plate,sep=""))
Dat$Map$Plot<-factor(paste("Plot",Dat$Map$Plot,sep=""),levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
                                                                ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13","Plot14","Plot15","Plot16","Plot17","Plot18"))


##Plot usable reads
p<-chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",col_val = "Genotype")
oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot_fungi.svg")

#Remove sapmles with less than 3000 reads
Dat <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000 ])
Dat<-subset.Dataset(x = Dat,subset = Fraction!="Blank",drop = T,clean = T)
Dat <- clean(Dat)
p<-chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",col_val = "Genotype")
oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot_fungi_filtered3000.svg")



##Adjust the levels
Dat$Map$Condition<-factor(Dat$Map$Condition,levels=c("low","medium","high","low+Pi"))
Dat$Map$Genotype<-factor(Dat$Map$Genotype,levels=c("No_plant","Col-0","phf1","phr1/phl1"))
Dat$Map$Fraction<-factor(Dat$Map$Fraction,levels=c("BulkSoil","Soil","Root","Shoot"))


### Rarefaction
set.seed(3111)
Dat_rar <- rarefaction(x = Dat,sample = 3000)
Dat_rar <- clean(Dat_rar)

#Relative abundance
Dat_rab<-Dat
Dat_rab$Tab<-scale(x = Dat_rab$Tab,center = F,scale = colSums(Dat_rab$Tab))
Dat_rab<-clean(Dat_rab)

Dat_raw<-Dat
##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_fungi<-list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_rab)
filename <- paste(outfolder,"dat_hallepi_fungi_otus.RDS",sep = "/")
saveRDS(object = Dat_fungi,file = filename)
rm(list=ls())

#####Create Dataset for Bacterial ESV#####

outfolder<-paste0("../cleandata/bacteria_esv")
dir.create(outfolder)
outfigures<-"../figures/bacteria_esv"
dir.create(outfigures)

#Read dada2 seq object
Tab<-readRDS(file = "../rawdata/seqtab_final_bacteria_esv.rds")
ids<-data.frame(Id=paste0("ESV",1:ncol(Tab)),Seq=colnames(Tab))
colnames(Tab)<-ids$Id
rows_new<-gsub(pattern = ".fastq.filt.fastq.gz",replacement = "",x = rownames(Tab))
rownames(Tab)<-rows_new
Tab<-t(Tab)
#Read metadata
Map <- read.table("../rawdata/metadata_phosphate_gradient_hallesoil.tsv",header=T,sep="\t")
Map$DADA2_Header<-factor(paste0("Plate",Map$Plate,"_",Map$PosWell,Map$PosNum))

#Rename the  rownames in the Metadata
rownames(Map) <- Map$DADA2_Header

# Match tab and map
Map <- Map[ (row.names(Map) %in% colnames(Tab)), ]
Tab <- Tab[ ,colnames(Tab) %in% row.names(Map)]
Map <- Map[colnames(Tab),]

#Handle the taxtonomic object derived from dada2
Tax<-readRDS(file = "../rawdata/tax_final_bacteria_esv.rds")
Tax<-as.data.frame(x = Tax)
Tax$Taxonomy<-factor(paste0("Root;k__",Tax$Kingdom,";p__",Tax$Phylum,";c__",Tax$Class,
                            ";o__",Tax$Order,";f__",Tax$Family,";g__",Tax$Genus,";s__",Tax$Species))

Tax<-droplevels(subset(Tax,Kingdom=="Bacteria"))

aTax<-data.frame(ID=ids$Id[match(rownames(Tax),ids$Seq)],Taxonomy=Tax$Taxonomy,Seq=ids$Seq[match(rownames(Tax),ids$Seq)])
rownames(aTax)<-aTax$ID

Tab<-Tab[which(rownames(Tab)%in%aTax$ID),]

Dat <- create_dataset(Tab = Tab,Map = Map,Tax = aTax)
Dat$Map$Depth <- colSums(Dat$Tab)

#Add Plate label to the Plate column in the Metadata
Dat$Map$Plate<-factor(paste("Plate",Dat$Map$Plate,sep=""))
Dat$Map$Plot<-factor(paste("Plot",Dat$Map$Plot,sep=""),levels=c("PlotBlank","Plot1","Plot2","Plot3","Plot4","Plot5","Plot6","Plot7"
                                                                ,"Plot8","Plot9","Plot10","Plot11","Plot12","Plot13","Plot14","Plot15","Plot16","Plot17","Plot18"))

##Plot usable reads
p<-chibi.boxplot(Map = Dat$Map,x_val = "Fraction",y_val = "Depth",col_val = "Genotype")
oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot_bacteria_esv.svg")

#Remove sapmles with less than 5000 reads
Dat <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 3000 ])
#Remove the Blank
Dat<-subset.Dataset(x = Dat,subset = Fraction!="Blank",drop = T,clean = T)
Dat <- clean(Dat)

oh.ggsave.svg(outdir = outfigures,ggobject = p,outname = "usable_reads_boxplot_bacteria_esv_filtered3000.svg")

##Adjust the levels
Dat$Map$Condition<-factor(Dat$Map$Condition,levels=c("low","medium","high","low+Pi"))
Dat$Map$Genotype<-factor(Dat$Map$Genotype,levels=c("No_plant","Col-0","phf1","phr1/phl1"))
Dat$Map$Fraction<-factor(Dat$Map$Fraction,levels=c("BulkSoil","Soil","Root","Shoot"))


#Rarefaction
set.seed(3111)
Dat_rar <- rarefaction(x = Dat,sample = 3000)
Dat_rar <- clean(Dat_rar)

#Relative abundance
Dat_rab<-Dat
Dat_rab<-clean(Dat_rab)
Dat_rab$Tab<-scale(x = Dat_rab$Tab,center = F,scale = colSums(Dat_rab$Tab))

Dat_raw<-Dat
##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_esv<-list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_rab)
filename <- paste(outfolder,"dat_hallepi_bacteria_esv.RDS",sep = "/")
saveRDS(object = Dat_esv,file = filename)


