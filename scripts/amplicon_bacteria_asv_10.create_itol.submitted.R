library(ohchibi)
library(data.tree)
library(paletteer)
library(palettesPM)


df_res <- read.table(file = "../cleandata/df_dds_res_amplicon_bacteria_asv_fraction_asvlevel.tsv",
                     header = T,sep = "\t",quote = "",comment.char = "",stringsAsFactors = F)

#Create phylogenetic tree
df_res$Phylum[is.na(df_res$Phylum)] <- "Unknown"
df_res$Class[is.na(df_res$Class)] <- "Unknown"
df_res$Order[is.na(df_res$Order)] <- "Unknown"
df_res$Family[is.na(df_res$Family)] <- "Unknown"
df_res$Genus[is.na(df_res$Genus)] <- "Unknown"
#Clean the shigella so itol does not break
df_res$Genus <- df_res$Genus %>%
  gsub(pattern = "\\/",replacement = "") 


df_res$Phylum <- paste0("p_",df_res$Phylum)
df_res$Class <- paste0(df_res$Phylum,"|c_",df_res$Class)
df_res$Order <- paste0(df_res$Class,"|o_",df_res$Order)
df_res$Family <- paste0(df_res$Order,"|f_",df_res$Family)
df_res$Genus <- paste0(df_res$Family,"|g_",df_res$Genus)
df_res$ASV_Id_Long <- paste0(df_res$Genus,"|",df_res$ASV_Id)
df_res$ASV_Id_Long <- df_res$ASV_Id_Long %>% 
  gsub(pattern = "\\(|\\)",replacement = "")
df_res$pathString <- paste0("Root/",df_res$Phylum,"/",df_res$Class,"/",
                            df_res$Order,"/",df_res$Family,"/",df_res$Genus,"/",
                            df_res$ASV_Id_Long)
df_res$pathString <- df_res$pathString %>% 
  gsub(pattern = "\\(|\\)",replacement = "")

# Create phylogenetic tree
df_tax <- df_res %>% subset(Contrast == "Root-Soil") %>%
  droplevels

#Clean for weird characters so the tree costruction works
population <- as.Node(df_tax)
tree <- data.tree::as.phylo.Node(population)
outdir <- "../cleandata/itol"
dir.create(outdir)
write.tree(phy = tree,file = "../cleandata/itol/3874_asvs.newick")



#Colored ranges
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_colored_ranges.txt"
mline <- c("TREE_COLORS","Separator TAB",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")

df_temp <- df_res %>% subset(Contrast == "Root-Soil") %>%
  droplevels
ids <- df_temp$Phylum %>% unique %>% na.omit %>% as.character %>%
  gsub(pattern = "p_",replacement = "")
ids_c <- ids
ids_c[which(!(ids_c %in% pm.names.phyla()))] <- "Other"

phyla_colors <- pm.colors.phyla()
paleta <- match(ids_c,names(pm.colors.phyla())) %>% phyla_colors[.] %>%
  as.vector
df <- data.frame(Id = paste0("p_",ids), type = rep("clade",length(ids)),color = paleta,
                 ltype = rep("normal",length(ids)), width = rep(2,length(ids)))


write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Dataset Style Colors
#Colored ranges
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_labelbranch_colors.txt"
mline <- c("DATASET_STYLE","Separator TAB","DATASET_LABEL\tPhylum Labels","COLOR\tred",
           "LEGEND_TITLE\tPhylum","LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1",
           "LEGEND_COLORS\t#41F0AC\t#0000C0\t#FFFF00\t#008000\t#BCFFEE\t#B856D7\t#FF0D00\t#8B4513\t#FF8000\t#8D8DFF\t#000000",
           "LEGEND_LABELS\tAcidobacteria\tActinobacteria\tBacteroidetes\tChloroflexi\tCyanobacteria\tFirmicutes\tGemmatimonadetes\tPatescibacteria\tProteobacteria\tVerrucomicrobia\tOther",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")

df_temp <- df_res %>% subset(Contrast == "Root-Soil") %>%
  droplevels
df_temp <- data.frame(Id = df_temp$ASV_Id_Long,
                      Phylum = df_temp$Phylum %>% gsub(pattern = "p_",replacement = ""),stringsAsFactors = F)
df_temp$Phylum[which(!(df_temp$Phylum %in% pm.names.phyla()))] <- "Other"
phyla_colors <- pm.colors.phyla()
paleta <- match(df_temp$Phylum,names(pm.colors.phyla())) %>% phyla_colors[.] %>%
  as.vector
df_temp$Colors <- paleta
df <- data.frame(Id = df_temp$Id,Type = rep("branch",nrow(df_temp)),What = rep("node",nrow(df_temp)),
                 Color = df_temp$Colors,Size = rep(2,nrow(df_temp),Line = rep("normal",nrow(df_temp))))
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Define labels size as 0
df <- data.frame(Id = df_temp$Id,Type = rep("label",nrow(df_temp)),What = rep("node",nrow(df_temp)),
                 COLOR = rep("white",nrow(df_temp)),Size = rep(0,nrow(df_temp)),Style = rep("normal",nrow(df_temp)))
write.table(x = df,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)



#Heatmap for Root soil 
df_temp <- df_res
df_temp$log2FoldChange[which(df_temp$padj > 0.1)] <- 0
df_temp$log2FoldChange[which(is.na(df_temp$padj))] <- 0
df_temp$log2FoldChange[which(df_temp$log2FoldChange>= 10)] <- 10
df_temp$log2FoldChange[which(df_temp$log2FoldChange<=-10)] <- -10
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_heatmap_rootsoil.txt"
mline <- c("DATASET_HEATMAP","Separator TAB","DATASET_LABEL\tRootvsSoil log2FoldChange","COLOR\t#B27612","FIELD_LABELS\tRoot vs Soil",
           "STRIP_WIDTH\t50","MARGIN\t10","SIZE_FACTOR\t30","AUTO_LEGEND\t0",
           "LEGEND_TITLE\tlog2FoldChange",
           "COLOR_MIN\t#123AB2","COLOR_MAX\t#FF3503","USE_MID_COLOR\t1","COLOR_MID\t#D9D9D9","USER_MID_VALUE\t0",
           "USER_MIN_VALUE\t-10","USER_MAX_VALUE\t10",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df_lf <- df_temp %>% subset(Contrast == "Root-Soil") %>%
  droplevels 
df_lf <- df_lf[,c(15,3)]
write.table(x = df_lf,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

#Heatmap for Shoot soil 
df_temp <- df_res 
df_temp$log2FoldChange[which(df_temp$padj > 0.1)] <- 0
df_temp$log2FoldChange[which(is.na(df_temp$padj))] <- 0
df_temp$log2FoldChange[which(df_temp$log2FoldChange>= 10)] <- 10
df_temp$log2FoldChange[which(df_temp$log2FoldChange<=-10)] <- -10
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_heatmap_shootsoil.txt"
mline <- c("DATASET_HEATMAP","Separator TAB","DATASET_LABEL\tShootvsSoil log2FoldChange","COLOR\t#29AB5C","FIELD_LABELS\tShoot vs Soil",
           "STRIP_WIDTH\t50","MARGIN\t10","SIZE_FACTOR\t30","AUTO_LEGEND\t0",
           "LEGEND_TITLE\tlog2FoldChange",
           "COLOR_MIN\t#123AB2","COLOR_MAX\t#FF3503","USE_MID_COLOR\t1","COLOR_MID\t#D9D9D9","USER_MID_VALUE\t0",
           "USER_MIN_VALUE\t-10","USER_MAX_VALUE\t10",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df_lf <- df_temp %>% subset(Contrast == "Shoot-Soil") %>%
  droplevels 
df_lf <- df_lf[,c(15,3)]
write.table(x = df_lf,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)
# 
#Heatmap for Shoot Root
df_temp <- df_res
df_temp$log2FoldChange[which(df_temp$padj > 0.1)] <- 0
df_temp$log2FoldChange[which(is.na(df_temp$padj))] <- 0
df_temp$log2FoldChange[which(df_temp$log2FoldChange>= 10)] <- 10
df_temp$log2FoldChange[which(df_temp$log2FoldChange<=-10)] <- -10
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_heatmap_shootroot.txt"
mline <- c("DATASET_HEATMAP","Separator TAB","DATASET_LABEL\tRootvsShoot log2FoldChange","COLOR\tpurple","FIELD_LABELS\tRoot vs Shoot",
           "STRIP_WIDTH\t50","MARGIN\t10","SIZE_FACTOR\t30","AUTO_LEGEND\t0",
           "LEGEND_TITLE\tlog2FoldChange",
           "COLOR_MIN\t#123AB2","COLOR_MAX\t#FF3503","USE_MID_COLOR\t1","COLOR_MID\t#D9D9D9","USER_MID_VALUE\t0",
           "USER_MIN_VALUE\t-10","USER_MAX_VALUE\t10",
           "DATA")
write(x = mline,file = outfile,append = F,sep = "")
df_lf <- df_temp %>% subset(Contrast == "Root-Shoot") %>%
  droplevels
df_lf <- df_lf[,c(15,3)]
write.table(x = df_lf,file = outfile,append = T,quote = F,sep = "\t",row.names = F,col.names = F)

### Popup information ###
df_tax <- df_res 

#Define function to get popup information 
getpopupinfo <- function(df = df_tax, contrasts = NULL,thres_padj = 0.1){
  ids <- df$pathString %>% unique
  Res_contrasts  <- NULL
  Res_num <- NULL
  for(id in ids){
    nt <- df %>% subset(pathString == id) %$% ASV_Id %>% unique %>% length
    inds <- id %>% strsplit(split = "\\/") %>% unlist
    ide_raw <- inds %>% length %>% inds[.]
    ide_clean <- ide_raw %>% strsplit(split = "\\|") %>% unlist
    ide_clean <- ide_clean %>% length %>% ide_clean[.]
    ide_clean <- ide_clean %>% gsub(pattern = "^[p|c|o|f|g]_",replacement = "")
    Res_num <- data.frame(Id = ide_raw,Label =  ide_clean, Total = nt) %>%
      rbind(Res_num,.)
    for(contrast in contrasts){
      df_temp <-  df %>% subset(pathString == id) %>%
        droplevels
      df_temp <- df_temp %>% subset(Contrast == contrast) %>% droplevels
      ntotal <- df_temp$ASV_Id %>% unique %>% length
      enriched <- df_temp %>% subset(padj < thres_padj & log2FoldChange >0) %>%
        nrow
      depleted <- df_temp %>% subset(padj < thres_padj & log2FoldChange <0) %>%
        nrow
      ns <- ntotal -enriched -depleted
      Res_contrasts <- data.frame(Id = ide_raw,Label =  ide_clean, Contrast = contrast,
                                  Enriched = enriched, Depleted = depleted, NS =  ns, Total = ntotal) %>%
        rbind(Res_contrasts,.)
      
    }
  }
  lRes <- list(Res_num = Res_num, Res_contrasts = Res_contrasts)
  return(lRes)
}
contrasts <- c("Root-Soil","Shoot-Soil","Root-Shoot")

#ASV Level
pop_asv  <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Genus Level
nst_vec <- NULL
for(st in df_tax$pathString){
  temp <- st %>% strsplit(split = "\\/") %>% unlist
  nst <- temp[-length(temp)] %>% paste(collapse = "/")
  nst_vec <- c(nst_vec,nst)
}
df_tax$pathString <- nst_vec
pop_genus <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Family Level
nst_vec <- NULL
for(st in df_tax$pathString){
  temp <- st %>% strsplit(split = "\\/") %>% unlist
  nst <- temp[-length(temp)] %>% paste(collapse = "/")
  nst_vec <- c(nst_vec,nst)
}
df_tax$pathString <- nst_vec
pop_family <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Order Level
nst_vec <- NULL
for(st in df_tax$pathString){
  temp <- st %>% strsplit(split = "\\/") %>% unlist
  nst <- temp[-length(temp)] %>% paste(collapse = "/")
  nst_vec <- c(nst_vec,nst)
}
df_tax$pathString <- nst_vec
pop_order <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Class Level
nst_vec <- NULL
for(st in df_tax$pathString){
  temp <- st %>% strsplit(split = "\\/") %>% unlist
  nst <- temp[-length(temp)] %>% paste(collapse = "/")
  nst_vec <- c(nst_vec,nst)
}
df_tax$pathString <- nst_vec
pop_class <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Phylum Level
nst_vec <- NULL
for(st in df_tax$pathString){
  temp <- st %>% strsplit(split = "\\/") %>% unlist
  nst <- temp[-length(temp)] %>% paste(collapse = "/")
  nst_vec <- c(nst_vec,nst)
}
df_tax$pathString <- nst_vec
pop_phylum <- getpopupinfo(df = df_tax,contrasts = contrasts)

#Merge all the tables into structure to create the itol file
Res_num <- rbind(pop_asv$Res_num,pop_genus$Res_num,pop_family$Res_num,
                 pop_order$Res_num,pop_class$Res_num,pop_phylum$Res_num)


Res_contrasts <- rbind(pop_asv$Res_contrasts,pop_genus$Res_contrasts,pop_family$Res_contrasts,
                       pop_order$Res_contrasts,pop_class$Res_contrasts,pop_phylum$Res_contrasts)

#Create files for popup info for Number of elements
outfile <- "../cleandata/itol/itol_amplicon_bacteria_asv_asvlevel_popup.txt"
mline <- c("POPUP_INFO","Separator TAB","DATA")
write(x = mline,file = outfile,append = F,sep = "")
for(i in 1:nrow(Res_num)){
  line_num <- Res_num[i,] 
  lines_contrasts <- which(Res_contrasts$Id %in% line_num$Id) %>% Res_contrasts[.,]
  #Print general information about Name,Taxonomy and num of ASVs
  cat(as.character(line_num$Id),"\tMetadata\t","<h1>",as.character(line_num$Label),"</h1>",
      "<h3>Taxonomy: ",as.character(line_num$Id),"</h3>",
      "<h3>Number of ASVs: ",line_num$Total,"</h3>",
      sep = "",file = outfile,append = TRUE)
  #Print information about the enrichments 
  for(j in 1:nrow(lines_contrasts)){
    line_c <- lines_contrasts[j,]
    cat("<h4>Contrast: <b>",as.character(line_c$Contrast),"</b> | ",
        gsub(pattern = "-.*",replacement = "",as.character(line_c$Contrast)),"(Number ASV Enriched): ",line_c$Enriched," | ",
        gsub(pattern = ".*-",replacement = "",as.character(line_c$Contrast)),"(Number ASV Enriched): ",line_c$Depleted," | ",
        "NotSignificant(Number ASV): ",line_c$NS,"</h4><br>",
        sep = "",file = outfile,append = TRUE)
  }
  
  cat(file = outfile,sep = "\n",append = TRUE)
  
}

