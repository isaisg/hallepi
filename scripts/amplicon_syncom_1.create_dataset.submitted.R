library(ohchibi)
library(gridExtra)


#Read the data
Tab_usearch <- read.table(file = "../rawdata/matrix_id98_phosphate_syncom.txt",
                          header = T,sep = "\t",row.names = 1,
                          quote = "",check.names = F,comment.char = "")

Dat_otus <-read.am(file = "../rawdata/otu_table_phosphate_syncom.txt",
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


#Remove the ones that hit to contaminants
Tab_usearch <- Tab_usearch[,-grep(pattern = "contaminant",colnames(Tab_usearch))]

#Remove the unmapped
Tab_usearch <- Tab_usearch[,-which(colnames(Tab_usearch)=="Unmapped")]


#Read the Metadata
Map <- read.table(file = "../rawdata/metadata_phosphate_syncom.txt",
                  header = T,sep = "\t",quote = "",comment.char = "")



#Use the Map to put both tables in the same order
Tab_usearch <- match(Map$ID_Matrix,rownames(Tab_usearch)) %>% 
  Tab_usearch[.,]

Tab_otus <- match(Map$OTU_Header,colnames(Tab_otus)) %>% 
  Tab_otus[,.] %>% t
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


#Compute the Depth per Sample
Dat$Map$Depth <- colSums(Dat$Tab)

#Compute the Log Depth per Sample
Dat$Map$LDepth <-log(colSums(Dat$Tab))
Dat$Map$Pi <- Dat$Map$condition %>% gsub(pattern = "P",replacement = "") %>%
  gsub(replacement = NA,pattern = "Inoculum") %>% as.numeric

#Readjust the typebyTissue 
Dat$Map$Rep <- factor(gsub(pattern = "R",replacement = "Rep",x = as.character(Dat$Map$Rep)))
Dat$Map$condition <- Dat$Map$condition %>% gsub(pattern = "P",replacement = " micromolar Pi") %>% 
  factor(levels = c("Inoculum","0 micromolar Pi","10 micromolar Pi","30 micromolar Pi",
                    "50 micromolar Pi","100 micromolar Pi","1000 micromolar Pi"))


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
filename <- "../cleandata/dat_hallepi_syncom_bacteria_useq97.RDS"
saveRDS(object = Dat_amplicon,file = filename)

