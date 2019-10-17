library(ohchibi)
library(gridExtra)

set.seed(130816)
setwd('/home/isai/Documents/results/hallepi/revision_plosbiology/scripts')

#Read the data
Tab_usearch <- read.table(file = "../rawdata/hallepi_syncom_pido_useq.tsv",
                          header = T,sep = "\t",row.names = 1,
                          quote = "",check.names = F,comment.char = "")

Dat_otus <-read.am(file = "../rawdata/hallepi_syncom_pido_otus.txt",
                   format = "qiime",taxonomy = "taxonomy")
Tab_otus <- Dat_otus$Tab
Tax_otus <- Dat_otus$Tax

#Compute proportion of Mapped to strain 
usearch_counts <- data.frame((colSums(Tab_usearch)/sum(Tab_usearch))*100)
colnames(usearch_counts) <- "Global_abundance"
usearch_counts$Strain <- rownames(usearch_counts)
percentage_mapped_to_strains <- sum(usearch_counts[grep(pattern = "Sequence_",x = usearch_counts$Strain),][,1])
percentage_mapped_to_contaminants <- sum(usearch_counts[grep(pattern = "contaminant",x = usearch_counts$Strain),][,1])
percentage_unmapped <- sum(usearch_counts[grep(pattern = "Unmapped",x = usearch_counts$Strain),][,1])

###Create a table with the mapping data results
df<-data.frame(Category=c("Reads mapped to strains","Reads mapped to contaminants","Reads Unmapped"),
               Percentage=c(percentage_mapped_to_strains,percentage_mapped_to_contaminants,percentage_unmapped))
grid.table(df,rows=NULL)
dev.off()

#Remove the ones that hit to contaminants
Tab_usearch <- Tab_usearch[,-grep(pattern = "contaminant",colnames(Tab_usearch))]

#Remove the unmapped
Tab_usearch <- Tab_usearch[,-which(colnames(Tab_usearch)=="Unmapped")]



#Read the Metadata
Map <- read.table(file = "../rawdata/hallepi_syncom_pido_counts_16s_metadata.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")

Map <- which(!(is.na(Map$sample))) %>%
  Map[.,] %>% droplevels

df_trans <- data.frame (Number = paste0("P",1:96),
                        Well = c("A1","B1","C1","D1","E1","F1","G1","H1",
                                 "A2","B2","C2","D2","E2","F2","G2","H2",
                                 "A3","B3","C3","D3","E3","F3","G3","H3",
                                 "A4","B4","C4","D4","E4","F4","G4","H4",
                                 "A5","B5","C5","D5","E5","F5","G5","H5",
                                 "A6","B6","C6","D6","E6","F6","G6","H6",
                                 "A7","B7","C7","D7","E7","F7","G7","H7",
                                 "A8","B8","C8","D8","E8","F8","G8","H8",
                                 "A9","B9","C9","D9","E9","F9","G9","H9",
                                 "A10","B10","C10","D10","E10","F10","G10","H10",
                                 "A11","B11","C11","D11","E11","F11","G11","H11",
                                 "A12","B12","C12","D12","E12","F12","G12","H12")
)

Map$Number <- Map$position %>% gsub(pattern = "Plate.*",replacement = "")
Map <- merge(Map,df_trans, by = "Number")
Map$Plate <- Map$position %>% gsub(pattern = "^P[0-9]+",replacement = "")
Map$OTU_Header <- paste0(Map$position,"PiDOR1")
Map$ID_Matrix <-  paste0(Map$Plate,"PiDOR1",Map$Well)



#Use the Map to put both tables in the same order
Tab_usearch <- match(Map$ID_Matrix,rownames(Tab_usearch)) %>% 
  Tab_usearch[.,]

Tab_otus <- match(Map$OTU_Header,colnames(Tab_otus)) %>% 
  Tab_otus[,.] %>% t
Tab_otus[is.na(Tab_otus)] <- 0
rownames(Tab_otus) <- Map$ID_Matrix
Tab_merged <- cbind(Tab_usearch,Tab_otus)

#Add the taxonomic information here
Tax_usearch <- read.table(file = "../rawdata/strains_4stresses_august2017_manuallycurated.fasta_unique_sequences.ng.wang.taxonomy",
                          header=F,sep="\t")

#Modify the structure to contain the proper qiime format needed by AMOR
temp_new<-NULL
for(element in Tax_usearch$V2){
  temp_vec<-unlist(strsplit(x = element,split = ";"))
  qiime_format<-paste("Root; k__",temp_vec[1],"; p__",temp_vec[2],"; c__",temp_vec[3],"; o__",temp_vec[4],"; f__",temp_vec[5],"; g__",temp_vec[6],sep="")
  temp_new<-c(temp_new,qiime_format)
}

Tax_usearch <- data.frame(ID=Tax_usearch$V1,Taxonomy=temp_new)
rownames(Tax_usearch) <- Tax_usearch$ID

#Merge taxonomies
Tax_merged<-rbind(Tax_usearch,Tax_otus)

Tax_merged <- match(colnames(Tab_merged),rownames(Tax_merged)) %>%
  Tax_merged[.,] %>% droplevels

#Create dataset
rownames(Map) <- Map$ID_Matrix
Dat <- create_dataset(Tab = t(Tab_merged),Map = Map,Tax = Tax_merged)
Dat$Map$Pi  <- Dat$Map$Pi %>% factor(levels = c("0uM","50uM","1000uM"))
Dat$Map$Genotype <- Dat$Map$Genotype %>% as.character %>% gsub(pattern = "Col",replacement = "Col-0") %>%
  gsub(pattern = "phf",replacement = "phf1") %>%
  gsub(pattern = "p/p",replacement = "phr1/phl1") 
Dat$Map$Genotype[which(is.na(Dat$Map$Genotype))] <- "AgarPlant"
Dat$Map$Genotype <- Dat$Map$Genotype %>% factor(levels = c("AgarPlant","Col-0","phf1","phr1/phl1"))

Dat$Map$Fraction <- Dat$Map$Fraction %>% gsub(pattern = "Agar",replacement = "AgarPlant") %>% 
  factor
#Compute the Depth per Sample
Dat$Map$Depth <- colSums(Dat$Tab)

#Compute the Log Depth per Sample
Dat$Map$LDepth <-log(colSums(Dat$Tab))

#Remove samples with less than 1000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])

#Rarefaction
set.seed(130817)
#Remove the OTUs
#Dat_temp <- remove_taxons(Dat = Dat,taxons = as.character(Dat$Tax$ID[ grep(pattern = "OTU",x = Dat$Tax$ID) ]))
Dat_temp <- Dat
Dat_temp$Map$Counts_Rar<-colSums(Dat_temp$Tab)
Dat_temp <- subset(Dat_temp, Counts_Rar > 1000, drop = TRUE, clean = TRUE)
Dat_rar <- rarefaction(x = Dat_temp,sample = 1000)
Dat_rar <- clean(Dat_rar)


###Relative abundance 
Dat_ra <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat_ra$Tab <- scale(x = Dat_ra$Tab,center = F,scale = colSums(Dat_ra$Tab))


##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_amplicon<-list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_ra)
filename <- "../cleandata/dat_hallepi_syncom_piDO_useq97.RDS"
saveRDS(object = Dat_amplicon,file = filename)
rm(list=ls())
