#libraries
library(dplyr)
library(ggplot2)
library(vroom)
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(plyr)
library(dslabs)
library(ggbreak)
library(circlize)
library(patchwork)
library(factoextra)
library(rlist)
library(ComplexUpset)
library(UpSetR)
library(grid)
library(gridExtra)
library(magick)
library(forcats)
library(ggbeeswarm)
library(readxl)
library(ggmap)
library(ggtree)
library(treeio)
library(plotly)
library(cowplot)
library(ggtreeExtra)
library(GenomicRanges)
library(ggnewscale)
library(scatterpie)
library(gggenes)
library(rnaturalearth)
library(ape)
library(phytools)

#FUNCTIONS
#Necessary for merging SVs close to each other and similar in size (listed already in the mergl list) in the PCA data frame
update_sv_pca <- function(sv_pca, mergl) {
  for (col in names(mergl)) {
    if (col %in% colnames(sv_pca)) {
      linked_cols <- mergl[[col]]
      # Filter linked_cols to include only those that exist in sv_pca
      existing_linked_cols <- linked_cols[linked_cols %in% colnames(sv_pca)]
      if (length(existing_linked_cols) > 0) {
        sv_pca <- sv_pca %>%
          mutate(!!sym(col) := ifelse(rowSums(select(., all_of(existing_linked_cols))) > 0, 1, !!sym(col))) %>%
          select(-all_of(existing_linked_cols))
      }
    }
  }
  return(sv_pca)
}
#Necessary for splitting the gff file into genes we can add to a plot
parse_info <- function(info_string) {
  # Split the string by semicolon
  info_split <- strsplit(info_string, ";")[[1]]
  
  # Split each part into key-value pairs and remove whitespace
  info_pairs <- strsplit(info_split, "=")
  info_pairs <- lapply(info_pairs, function(pair) {
    if(length(pair) == 2) return(setNames(pair[2], trimws(pair[1])))
    return(NULL)
  })
  
  # Combine into a single named list
  info_list <- do.call(c, info_pairs)
  
  return(info_list)
}
#Useful for creating a dataframe based on SVtype and genomic region you want to scan
svfe_type_rbind = function(df,type,start=1,end=4410532,move=5000,window=50000) {
  for (i in seq(start,end,move)) {
    df = rbind(df,c(i,i+window,as.numeric(sum(sv_type_mcnv[sv_type_mcnv$POS>i & sv_type_mcnv$POS<i+window & sv_type_mcnv$SVTYPE==type,]$n)),type,NA))
    df$Start = as.numeric(df$Start)
    df$End = as.numeric(df$End)
    df$n = as.numeric(df$n)
    df$SVTYPE = as.character(df$SVTYPE)
  }
  df[df$SVTYPE==type,]=df[df$SVTYPE==type,] %>% mutate(SVFE=n/mean(n))
  df$SVFE = as.numeric(df$SVFE)
  return(df)
}
# Function to join SVs and create mCNV nomenclature
join_svs_create_mcnv <- function(sv_type) {
  i <- 1
  while (i < nrow(sv_type)) {
    dup <- 0
    merg <- c()
    n <- i + 1
    while (n <= nrow(sv_type)) {
      # Check for NA values and break if conditions are not met
      if (is.na(sv_type[n, 2]) || is.na(sv_type[i, 2]) || is.na(sv_type[i, 4]) || is.na(sv_type[n, 4])) {
        break
      }
      
      if (sv_type[n, 2] > sv_type[i, 2] + 25) {
        if (length(merg) > 0) {
          mergl[[sv_type[i, 1]]] <- merg
        }
        break
      }
      
      if (sv_type[n, 3] == "DUP" || sv_type[i, 3] == "DUP") {
        dup <- 1
      }
      
      if ((sv_type[i, 4] + 50 < sv_type[n, 4] && dup == 0) || (sv_type[i, 4] - 50 > sv_type[n, 4] && dup == 0)) {
        if (length(merg) > 0) {
          mergl[[sv_type[i, 1]]] <- merg
        }
        break
      }
      
      if (dup == 1) {
        sv_type[i, 3] <- "CNV"
        sv_type[i, 1] <- gsub("DEL|INS|DUP", "CNV", sv_type[i, 1])
        sv_type[i, 4] <- max(sv_type[i, 4], sv_type[n, 4])
      }
      
      sv_type[i, 5] <- sv_type[i, 5] + sv_type[n, 5]
      merg <- c(merg, sv_type[n, 1])
      sv_type <- sv_type[-n, ]
      n <- n - 1  # Adjust index after row removal
      n <- n + 1
    }
    i <- i + 1
  }
  return(list(sv_type = sv_type, mergl = mergl))
}

#this function will help us find homoplasic events across the phylogeny
detect_homoplasy <- function(tree, sv_matrix) {
  homoplasic_svs <- list()
  for (sv in colnames(sv_matrix)) {
    pattern <- sv_matrix[, sv]
    names(pattern) <- rownames(sv_matrix)
    
    # Perform ancestral state reconstruction
    recon <- ace(pattern, tree, type = "discrete", model = "ER")
    
    # Check for homoplasy
    lik_anc <- recon$lik.anc
    homoplasy_detected <- any(lik_anc[, 1] > 0.4 & lik_anc[, 1] < 0.6) || any(lik_anc[, 2] > 0.4 & lik_anc[, 2] < 0.6)
    
    if (homoplasy_detected) {
      homoplasic_svs[[sv]] <- recon
    }
  }
  return(homoplasic_svs)
}

# Function to check overlap and increment count
overlap_count <- function(gene_start, gene_end, sv_start, sv_end) {
  # Check if there is any overlap
  overlap <- (sv_start <= gene_end) & (sv_end >= gene_start)
  return(overlap)
}

esx_plot <- function(gene,genes_to_plot) {
  return(essentially %>% mutate(Plot = case_when(Name %in% genes_to_plot ~ Name,
                                          `Final Call` == "ES" ~ "ES",
                                          `Final Call` == "NE" ~ "NE"),
                         Color = case_when(Name %in% genes_to_plot ~ "Other",
                                           `Final Call` == "ES" ~ "ES",
                                           `Final Call` == "NE" ~ "NE"),
                         ESX = case_when(grepl("1",Name) ~ "ESX1",
                                         grepl("2",Name) ~ "ESX2",
                                         grepl("3",Name) ~ "ESX3",
                                         grepl("4",Name) ~ "ESX4",
                                         grepl("5",Name) ~ "ESX5",
                                         TRUE ~ "Other"),
                         sv_pop = replace(sv_pop,grepl("esxA|esxB|eccCb1|eccD1",Name), 0)) %>%
    
    filter(Name %in% genes_to_plot[grep(gene,genes_to_plot)]) %>%
    ggplot(.,aes(x=factor(Plot,levels=genes_to_plot[grep(gene,genes_to_plot)]),y=sv_pop,color = ESX)) +
    geom_point(alpha=0.7,size=3) +
    geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="ES",]$sv_pop,trim=0.1),color="#F08080",linetype="dashed") +
    geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="NE",]$sv_pop,trim=0.1),color="lightblue",linetype="dashed") +
    #geom_quasirandom() +
    scale_color_manual(values = c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD")) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(),legend.title = element_blank()) +
    #coord_cartesian(expand = FALSE) +
    ylab("Number of isolates with SVs"))
}

#We create the lineage data frame

lin = read.table("~/Documents/PHD/SV/SV_analysis/temp/final_lineages.txt",sep=" ")
illu = read.table("~/Documents/PHD/SV/SV_analysis/temp/ncbi_illu.txt",sep=" ")
colnames(lin) = c("ID","Lineage")
lin$L1 = NA
lin$L2 = NA
lin$L3 = NA
lin$L4 = NA
lin$L5 = NA


for (i in 1:nrow(lin)) {
  if (length(str_split(lin$Lineage,"[.]")[[i]]) == 1) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
  }
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 2) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 3) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 4) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
    lin[i,]$L4 = str_split(lin$Lineage,"[.]")[[i]][4]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 5) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
    lin[i,]$L4 = str_split(lin$Lineage,"[.]")[[i]][4]
    lin[i,]$L5 = str_split(lin$Lineage,"[.]")[[i]][5]
  }
}    
lin_color = c("#DDCC77","#117733","#88CCEE","#CC6677","#882255","#AA4499","#332288","#44AA99","#DDAA55")

#We load the vcf files 
files_vcf = list.files(path='~/Documents/PHD/SV/SV_analysis/vcf/longass/',pattern = ".*vcf")
files_vcf = files_vcf[!grepl(paste(illu$V1,collapse="|"),files_vcf)]
file.pathsv <<- as.vector(paste("~/Documents/PHD/SV/SV_analysis/vcf/longass", files_vcf, sep = "/"))
rep = read.table("~/Documents/PHD/SV/SV_analysis/peppe_sites.txt",sep=" ")
tablesv <- lapply(X = file.pathsv, FUN = read.table, header = FALSE,sep = "\t",comment.char='#')

#We now create a new data frame to add all the found SVs
sv_type = data.frame(ID="",POS="",SVTYPE="",SVLEN="",stringsAsFactors = FALSE)
sv_pca = data.frame(Lineage="",stringsAsFactors = FALSE)
for (n in 1:length(tablesv)) {
  for (m in 1:nrow(tablesv[[n]])) {
    if (is.na(strsplit(as.character(tablesv[[n]][m,3]),".",fixed = TRUE)[[1]][2])) {
      next
    }
    sv_type = rbind(sv_type,c(paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = "."),tablesv[[n]][m,2],strsplit(as.character(tablesv[[n]][m,3]),".",fixed = TRUE)[[1]][1],strsplit(as.character(tablesv[[n]][m,3]),".",fixed = TRUE)[[1]][2]))
    if (!(paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = ".") %in% colnames(sv_pca)) & grepl("Na",paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = "."),fixed=TRUE)==FALSE) {
      sv_pca[,paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = ".")] = 0
    }
  }
}

sv_type = sv_type[2:nrow(sv_type),]
sv_type_len = sv_type %>% add_count(SVLEN,SVTYPE) %>% distinct(SVLEN,SVTYPE,.keep_all = TRUE)
sv_type = sv_type %>% add_count(ID) %>% distinct(ID,.keep_all = TRUE)
sv_type = sv_type %>% mutate_at(c('POS','SVLEN'),as.numeric)
sv_type_len = sv_type_len %>% mutate_at(c('POS','SVLEN'),as.numeric)
sv_type = sv_type %>% mutate_at(c('POS','SVLEN','n'),as.numeric)

#We add mCNV and merge SVs that are close to each other
mergl = list()

sv_type2 = sv_type
sv_type = sv_type[order(sv_type$POS,sv_type$SVTYPE),]

# Apply the function
result <- join_svs_create_mcnv(sv_type)

# Extract results
sv_type_mcnv <- result$sv_type
mergl <- result$mergl

sv_type_mcnv = sv_type
merglcnv = mergl

#SV fold-enrichment figure
all_svs = data.frame(Start=as.numeric(),End=as.numeric(),n=as.numeric(),SVTYPE=as.character(),stringsAsFactors = FALSE)
for (i in seq(1,4410532,5000)) {
  all_svs = rbind(all_svs,c(i,i+50000,sum(sv_type_mcnv[sv_type_mcnv$POS>i & sv_type_mcnv$POS<i+50000,]$n),"All SVs"))
}
colnames(all_svs)=c("Start","End","n","SVTYPE")
all_svs$n = as.numeric(all_svs$n)
all_svs=all_svs %>% mutate(SVFE=n/mean(n))

all_svs = svfe_type_rbind(df=all_svs,type="DEL")
all_svs = svfe_type_rbind(df=all_svs,type="DUP")
all_svs = svfe_type_rbind(df=all_svs,type="INV")
all_svs = svfe_type_rbind(df=all_svs,type="INS")
all_svs = svfe_type_rbind(df=all_svs,type="CNV")

all_svs_z = data.frame(Start=as.numeric(),End=as.numeric(),n=as.numeric(),SVTYPE=as.character(),stringsAsFactors = FALSE)
for (i in seq(3900000,3975000,500)) {
  all_svs_z = rbind(all_svs_z,c(i,i+1000,sum(sv_type_mcnv[sv_type_mcnv$POS>i & sv_type_mcnv$POS<i+1000,]$n),"All SVs"))
  colnames(all_svs_z)=c("Start","End","n","SVTYPE")
  all_svs_z$n = as.numeric(all_svs_z$n)
  all_svs_z$Start = as.numeric(all_svs_z$Start)
  all_svs_z$End = as.numeric(all_svs_z$End)
  all_svs_z$SVTYPE = as.character(all_svs_z$SVTYPE)
}
all_svs_z[1,1] = 3900000
all_svs_z[1,2] = 3900500
all_svs_z[1,3] = 0
all_svs_z=all_svs_z %>% mutate(SVFE=n/mean(n))
all_svs_z = svfe_type_rbind(df=all_svs_z,type="DEL",start=3900000,end=3975000,move = 500,window = 1000)
all_svs_z = svfe_type_rbind(df=all_svs_z,type="DUP",start=3900000,end=3975000,move = 500,window = 1000)
all_svs_z = svfe_type_rbind(df=all_svs_z,type="INS",start=3900000,end=3975000,move = 500,window = 1000)
all_svs_z = svfe_type_rbind(df=all_svs_z,type="INV",start=3900000,end=3975000,move = 500,window = 1000)
all_svs_z = svfe_type_rbind(df=all_svs_z,type="CNV",start=3900000,end=3975000,move = 500,window = 1000)
#ALL SVs ACROSS THE MTB GENOME

svfe_plot = ggplot(all_svs,aes(x=Start,y=SVFE,color=SVTYPE)) +
  geom_point(size=0.5) +
  geom_line(size=3) +
  theme_pubr() +
  geom_area(aes(fill=SVTYPE),alpha=0.7,show.legend = FALSE) +
  coord_cartesian(ylim = c(0,NA),expand = FALSE) +
  geom_segment(aes(x=0,y=1,xend=4411532,yend=1),color="black") +
  xlab("Mtb Genome location") +
  ylab("SV fold-enrichment") +
  scale_color_manual(values=c("darkgray","#71E38C","#2376B2","#D43925","#D474E0","#FA931E")) +
  scale_fill_manual(values=c("darkgray","#71E38C","#2376B2","#D43925","#D474E0","#FA931E")) +
  theme(legend.title = element_blank(),strip.background = element_rect(fill = "white",color="white"),
        strip.text.x = element_blank(),axis.text = element_text(size = 30),axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),axis.ticks = element_line(linewidth = 3),legend.key.size = unit(3,"cm")) +
  facet_wrap(~SVTYPE,nrow=6,ncol=1,scales = "free")

#Now we can zoom in to a region rich in SVs and also look whether they affect any genes
genes = read.table("~/Downloads/Mycobacterium_tuberculosis_H37Rv_gff_v4.1.gff",sep=" ")
genes$V9 = as.character(genes$V8)
# Apply the function to the data frame and convert to a wide format
genes_names <- genes %>%
  dplyr::mutate(parsed_info = purrr::map(V9, parse_info)) %>%
  unnest_wider(parsed_info) %>%
  filter(Name != "ncRv3520")

gr_annotations <- GRanges(seqnames = Rle("NC_000962.3"), 
                          ranges = IRanges(start = genes_names$V4, end = genes_names$V5),
                          gene = genes_names$Name)

annotation_df <- as.data.frame(gr_annotations)
annotation_df$midpoint <- (annotation_df$start + annotation_df$end) / 2


svfe_zoom = ggplot(all_svs_z,aes(x=Start,y=SVFE,color=SVTYPE)) +
  geom_point(size=0.5) +
  geom_line(size=3) +
  theme_pubr() +
  geom_area(aes(fill=SVTYPE),alpha=0.7,show.legend = FALSE) +
  coord_cartesian(xlim=c(3920000,3960500),ylim = c(0,NA),expand = FALSE) +
  geom_segment(aes(x=0,y=1,xend=4411532,yend=1),color="black") +
  xlab("Mtb Genome location") +
  ylab("SV fold-enrichment") +
  scale_color_manual(values=c("darkgray","#71E38C","#2376B2","#D43925","#D474E0","#FA931E")) +
  scale_fill_manual(values=c("darkgray","#71E38C","#2376B2","#D43925","#D474E0","#FA931E")) +
  theme(legend.title = element_blank(),strip.background = element_rect(fill = "white",color="white"),
        strip.text.x = element_blank(),axis.title.x = element_blank(),legend.text = element_blank(),legend.key = element_blank(),
        legend.ticks = element_blank(),plot.margin = margin(0.5, 1.25, , , "cm"),axis.text.x = element_text(size = 25,angle = 5,vjust = 0.7, hjust = 0.7),axis.text.y = element_text(size = 30), 
        axis.title.y = element_text(size = 40), axis.ticks = element_line(linewidth = 3),legend.key.size = unit(3,"cm")) +
  guides(fill="none",color="none") +
  facet_wrap(~SVTYPE,nrow=6,ncol=1,scales = "free")

anno_plot = ggplot(annotation_df) +
  geom_segment(aes(x = start, xend = end, y = 0, yend = 0), 
               color = 'black', size = 2) +
  geom_text(aes(x = midpoint, y = -0.035, label = gene), 
            angle = 90, vjust = 0, hjust = 0, size=8) +
  coord_cartesian(xlim = c(3920000, 3960500), ylim = c(-0.05, 0), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )
zoom_anno = svfe_zoom / anno_plot + plot_layout(heights = c(6,1))
Fig2D = plot_grid(svfe_plot,zoom_anno,ncol=2,nrow=1,rel_widths = c(1.15,1))





#In this loop we won't create the mCNV as it is not informative at the isolate level (inside an mCNV having more or less copies can affect the phenotype)
sv_type = sv_type2
sv_type = sv_type[order(sv_type$POS,sv_type$SVTYPE),]
# First we have to cluster the bp distances
for (i in 2:nrow(sv_type)) {
  if (sv_type[i,2]-25 < sv_type[i-1,2]) {
    sv_type[i,2] = sv_type[i-1,2]
  }
}
#Now we can cluster based on SV size
sv_type = sv_type[order(sv_type$POS,sv_type$SVLEN),]
for (i in 1:nrow(sv_type)) {
  if (i == nrow(sv_type)) {
    break
  }
  dup=0
  ins=0
  del=0
  merg=c()
  for (n in i+1:nrow(sv_type)) {
    #print(sv_type[n,2])
    #print(sv_type[i,2]+25)
    if (sv_type[i+1,2] > sv_type[i,2]+25) {
      if (length(merg)>0) {
        mergl[[sv_type[i,1]]] =merg
      }
      break
    }
    if (sv_type[i,4]+50 < sv_type[i+1,4] || sv_type[i,4]-50 > sv_type[i+1,4]) {
      if (length(merg)>0) {
        mergl[[sv_type[i,1]]] =merg
      }
      break
    }
    if (sv_type[i,3]=="DEL" && sv_type[i+1,3]=="INS" || sv_type[i,3]=="DEL" && sv_type[i+1,3]=="INS" || sv_type[i+1,3]=="DEL" && sv_type[i,3]=="INS" || sv_type[i+1,3]=="DEL" && sv_type[i,3]=="DUP") {
      if (length(merg)>0) {
        mergl[[sv_type[i,1]]] =merg
      }
      break
    }
    #x=x+1
    sv_type[i,5] = sv_type[i,5] + sv_type[i+1,5]
    merg=append(merg,sv_type[i+1,1])
    n=i+1
    sv_type= sv_type[-n,]
  }
}

###
# COUNT AND SIZE PLOT
###
sv_type_len_mcnv = aggregate(sv_type_mcnv$n, by = list(sv_type_mcnv$SVTYPE, sv_type_mcnv$SVLEN), FUN = sum)
#sv_type_len_mcnv = sv_type_mcnv %>% add_count(SVLEN,SVTYPE) %>% distinct(SVLEN,SVTYPE,.keep_all = TRUE) 
colnames(sv_type_len_mcnv) = c("SVTYPE","SVLEN","nn")
sv_type_len_mcnv = sv_type_len_mcnv[order(sv_type_len_mcnv$SVLEN,sv_type_len_mcnv$SVTYPE),]
ctsz = data.frame(Start=as.numeric(),End=as.numeric(),n=as.numeric(),SVTYPE=as.character(),stringsAsFactors = FALSE)
for (i in seq(50,70000,50)) {
  ctsz = rbind(ctsz,c(i,i+50,sum(sv_type_len_mcnv[sv_type_len_mcnv$SVLEN>i & sv_type_len_mcnv$SVLEN<i+50 & sv_type_len_mcnv$SVTYPE=="CNV",]$nn),"CNV"))
}
for (i in seq(50,70000,50)) {
  ctsz = rbind(ctsz,c(i,i+50,sum(sv_type_len_mcnv[sv_type_len_mcnv$SVLEN>i & sv_type_len_mcnv$SVLEN<i+50 & sv_type_len_mcnv$SVTYPE=="DEL",]$nn),"DEL"))
}
for (i in seq(50,70000,50)) {
  ctsz = rbind(ctsz,c(i,i+50,sum(sv_type_len_mcnv[sv_type_len_mcnv$SVLEN>i & sv_type_len_mcnv$SVLEN<i+50 & sv_type_len_mcnv$SVTYPE=="DUP",]$nn),"DUP"))
}
for (i in seq(50,70000,50)) {
  ctsz = rbind(ctsz,c(i,i+50,sum(sv_type_len_mcnv[sv_type_len_mcnv$SVLEN>i & sv_type_len_mcnv$SVLEN<i+50 & sv_type_len_mcnv$SVTYPE=="INV",]$nn),"INV"))
}
for (i in seq(50,70000,50)) {
  ctsz = rbind(ctsz,c(i,i+50,sum(sv_type_len_mcnv[sv_type_len_mcnv$SVLEN>i & sv_type_len_mcnv$SVLEN<i+50 & sv_type_len_mcnv$SVTYPE=="INS",]$nn),"INS"))
}
colnames(ctsz)=c("Start","End","n","SVTYPE")
ctsz$n = as.numeric(ctsz$n)
ctsz$Start = as.numeric(ctsz$Start)
ctsz$End = as.numeric(ctsz$End)

Fig2A = ggplot(ctsz,aes(x=Start,y=n,color=SVTYPE,group=SVTYPE)) +
  theme_pubr() +
  geom_point() + 
  geom_line(size=4) +
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  #scale_x_continuous(breaks = seq(50,43450,by=1000)) +
  #scale_x_break(breaks = c(28000,37000)) +
  scale_color_manual(values=c("#71E38C","#2376B2","#D43925","#D474E0","#FA931E")) +
  scale_x_continuous(trans='log10',breaks = c(50,100,200,400,800,1500,3250,7000,15000,30000,70000)) +
  theme(legend.title = element_blank(),axis.text = element_text(size = 40),axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),axis.ticks = element_line(linewidth = 3),legend.key.size = unit(3,"cm")) +
  ylab("Number of SVs")+
  xlab("SV length")

#PCA
#With this loop we create a dataframe with all the SVs, if an isolate has the SV, it will be marked as 1, if not, 0
for (n in 1:length(tablesv)) {
  sv_pca[nrow(sv_pca)+1,] = 0
  rownames(sv_pca)[nrow(sv_pca)] = str_replace(files_vcf[n],"_ref.vcf","")
  sv_pca[nrow(sv_pca),1] = lin[lin$ID == str_replace(files_vcf[n],"_ref.vcf",""),]$L1
  for (m in 1:nrow(tablesv[[n]])) {
    if (paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = ".") %in% colnames(sv_pca)) {
      sv_pca[nrow(sv_pca),paste(tablesv[[n]][m,2],tablesv[[n]][m,3],sep = ".")] = 1
    }
  }
}

sv_pca = sv_pca[2:nrow(sv_pca),]


#we cluster SVs that are close to each other similar in size and the same type of SV (data already extracted from the previous loop)
sv_pca2 = sv_pca
x=0
sv_pca = update_sv_pca(sv_pca,mergl)

model = prcomp(sv_pca[,2:ncol(sv_pca)])
pca_output = data.frame(model$x,Lineage=sv_pca$Lineage)
pca_output$L8 = ""
for (i in 1:nrow(pca_output)) {
  if (lin[lin$ID==rownames(pca_output)[i],]$L1 == 8) {
    pca_output[i,"L8"] = 1
  }
  else {
    pca_output[i,"L8"] = 0
  }
}
fviz_pca_ind(model)
Fig2B=ggplot(pca_output,aes(x=PC1,y=PC2,color=Lineage,shape=L8,size=L8)) +
  geom_point() +
  theme_pubr() +
  scale_color_manual(values=lin_color) +
  ylab("PC2(8.8%)")+
  xlab("PC1(13%)") +
  scale_shape_manual(values=c(19,18)) +
  guides(shape=FALSE,size=FALSE,color=guide_legend(title = "Lineage")) +
  scale_size_manual(values = c(7,13)) +
  theme(axis.text = element_text(size = 40),axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),legend.title = element_text(size = 40),axis.ticks = element_line(linewidth = 3)) +
  guides(color = guide_legend(override.aes = list(size = 8)))

#Dim 2 1761789 holds TbD1 DEL/INS --> accounts for 2% of variability in Dim 2 (separates L2,3,4 from L1)
png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/pc2_contrib.png",height = 500,width = 500)
fviz_contrib(model, choice = "var", axes = 2, fill = "#66a182", sort.val = c("desc"), top = 10)
dev.off()

#SINGLETONS
singletons = data.frame("Isolate"=as.character(),"Lineage"=as.numeric(),"Type"=as.character(),stringsAsFactors = FALSE)
for (i in 1:nrow(sv_pca)) {
  for (x in 2:ncol(sv_pca)) {
    if (sv_pca[i,x]==1) {
      if (table(sv_pca[,x])[[2]]==1) {
        singletons = rbind(singletons,c(rownames(sv_pca)[i],sv_pca[i,1],"Singleton"))
        singletons[,2] = as.numeric(singletons[,2])
      }
      else if (length(table(sv_pca[sv_pca$Lineage==sv_pca[i,1],x]))>1 & length(table(sv_pca[sv_pca$Lineage!=sv_pca[i,1],x]))==1) {
        singletons = rbind(singletons,c(rownames(sv_pca)[i],sv_pca[i,1],"Lineage-specific"))
        singletons[,2] = as.numeric(singletons[,2])
      }
      else {
        singletons = rbind(singletons,c(rownames(sv_pca)[i],sv_pca[i,1],"Common"))
        singletons[,2] = as.numeric(singletons[,2])
      }
    }
  }
}
colnames(singletons) = c("Isolate","lineage","Type")
sing_p = ggplot(singletons,aes(x=fct_infreq(Isolate),fill=Type)) +
  geom_bar(color="transparent",width=2)+
  facet_grid(~lineage,space="free",scales="free") +
  theme_pubr()+
  theme(strip.background = element_rect(fill = "white",color="white"),strip.text.x = element_text(size=12),axis.text.x = element_blank(),axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.title = element_blank()) +
  scale_fill_manual(values = c("#FFDAB9","#CCCCFF","#AFEEEE")) +
  coord_cartesian(ylim = c(0,NA),expand=FALSE) +
  geom_hline(data = singletons %>%
               add_count(Isolate) %>%
               select(Isolate,lineage,n) %>%
               distinct() %>%
               group_by(lineage) %>%
               dplyr::summarise(mean_count = median(n)), aes(yintercept = mean_count,group=lineage),color="gray",linetype="dashed") +
  scale_y_continuous(breaks = seq(0,220,by=10))

sing = singletons[,c(1,2)] %>% group_by(Isolate) %>% count()


## PLOT MAP next to phylo tree
world <- ne_countries(scale = "medium", returnclass = "sf")
country=read.csv("~/Downloads/SV_metadata_country.csv",na.strings=c("","NA"))
coords=read.csv("~/Downloads/country-coord.csv",header = T)
count_lin = merge(country,sing,by.x=1,by.y=1) %>%
  merge(.,coords,by.x=2,by.y=1)
iso = read.table("~/Documents/PHD/SV/SV_analysis/temp/final_isolates_NN.txt")
count_lin=count_lin[count_lin$id%in%iso$V1,]
count_lin = count_lin[!grepl(paste(illu$V1,collapse="|"),count_lin$id),]
#Create a map object with all countries of the world
world_map <- map_data("world", exact=F)
country_lin = merge(count_lin[,-c(4)] %>% group_by(Location,lineage) %>% tally(),coords,by.x=1,by.y=1)
#country_lin[country_lin$Location=="United States" & country_lin$lineage==4,]$Longitude..average. = -115
country_lin_map = country_lin %>% pivot_wider(names_from = lineage, values_from = n)
country_lin_map[is.na(country_lin_map)] = 0
colnames(country_lin_map) = c("Location","Alpha.2.code","Alpha.3.code","Numeric.code","Latitude..average.","Longitude..average.","C","A","B","D","E","G","F","H","I")
country_lin_map$n = rowSums(country_lin_map[,7:15])/4
#Plot
map = ggplot() + 
  coord_fixed() + 
  xlab("") + 
  ylab("") + 
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group),fill="grey", size=0) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),legend.title = element_text(size = 20),
        panel.border=element_blank(),legend.text = element_text(size = 20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()) +
  labs(color="Lineage")+
  geom_scatterpie(aes(x=Longitude..average., y=Latitude..average., group=Location), data=country_lin_map,cols=LETTERS[1:9]) +
  #geom_point(data = country_lin[order(country_lin$n,decreasing=T),],aes(x=Longitude..average.,y=Latitude..average.,fill=as.factor(lineage),size=n),pch=21,position=position_jitter(width=2,height=2)) +
  scale_fill_manual(values=lin_color,labels=c("1","2","3","4","5","6","7","8","9")) +
  guides(fill=guide_legend(title="Lineage"))
#guides(fill=guide_legend(title="Lineage"))
#scale_size_continuous(range = c(8, 20))

##  
#PHYLO TREE
##

tree=read.newick("~/Documents/PHD/SV/SV_analysis/tree/RAxML_bestTree.mtb_graph_tree_sv_filtx2")
empty_label_nodes <- which(tree$tip.label == "")
if (length(empty_label_nodes) > 0) {
  tree <- drop.tip(tree, empty_label_nodes)
}
ancestral_sv = read.table("~/Documents/PHD/SV/SV_analysis/temp/ancestral_svs.tsv",sep="\t",header = FALSE)
anc_merg = merge(ancestral_sv,lin,by.x=1,by.y=1)
anc_merg = anc_merg[,1:5]
colnames(anc_merg)= c("Isolate","freq","lin","lineage","lineage2")
anc_merg=anc_merg[!grepl(paste(illu$V1,collapse="|"),anc_merg$Isolate),]
good_tree = tree
tips_remove= c("GCA_000735915.1","SRR12395420","SRR12395448","SRR12395076","SRR12395142","02MTB0534_S66_L001","02MTB0313_S32_L001",
               "02MTB0308_S91_L001","02MTB0558_S33_L001","02MTB0048","02MTB0593_S24_L001","02MTB0154_S148_L002","02MTB1214","02MTB0702_S190_L002",
               "02MTB0184_S115_L002","02MTB0585_S47_L001","02MTB1026","SRR12395378","SRR12395406","SRR12395432","SRR12395408","SRR12395365","SRR12395091",
               "SRR12395191","SRR6117388","GCA_000389905.1","GCA_932530315.1","02MTB0037_37_S40_L001","02MTB0409_S114_L002","SRR12395434","SRR12395403",
               "SRR12395073","SRR12395353","SRR12395456","SRR12395094","SRR12395454","SRR12395451","SRR12395312","SRR12395155","SRR12395226","ERR12328527",
               "SRR12395098","SRR12395175","02MTB0197_S164_L002","SRR12395151","SRR12395404","SRR12395425","ERR12328525","SRR12395436","SRR12395160",
               "SRR12395174","SRR12395133","02MTB1607","SRR12395419","SRR12395230","SRR12395407","SRR12395309","SRR12395439","SRR12395407","SRR12395377",
               "SRR12395430","SRR12395363","AUSMDU00099754","GCA_000277105.1","GCA_000422125.1","SRR12395054","SRR12395142","SRR12395222","SRR12395276",
               "SRR12395280","AUSMDU00099757","GCA_001997625.1","GCA_001998175.1","GCA_030285725.1","SRR12395153","SRR12395215","SRR12395382","SRR12395364")
pruned_tree = drop.tip(good_tree,tips_remove) %>% drop.tip(.,illu$V1)
#write.tree(pruned_tree,"~/Documents/PHD/SV/SV_analysis/tree/751_pruned_final_tree.nwk")
Fig2C = ggtree(pruned_tree,layout = "circular") %<+% lin + geom_tippoint(aes(colour=L1), size = 1) + scale_color_manual(values = lin_color)+
  new_scale_fill() +
  geom_fruit(data=anc_merg,geom=geom_col,mapping = aes(x=freq,y=Isolate,fill=as.factor(lineage)),
             axis.params=list(axis = "x",text.size= 1,nbreak= 3),
             grid.params=list()) +
  scale_fill_manual(values = lin_color) +
  guides(fill="none",color="none") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"))
#
ggsave(filename = "~/Downloads/try.png", plot = Fig2C, width=10, height=10)

x = image_read("~/Downloads/try.png")
y=image_trim(x)





#PLOT ISOLATES THAT HAVE HOMOPLASIC SV
sv_pca_order = sv_pca[pruned_tree$tip.label,-1]
zero_cols <- which(colSums(sv_pca_order == 0) == nrow(sv_pca_order))
sv_pca_order_nozero = sv_pca_order[,-zero_cols]
homoplasic_svs = detect_homoplasy(pruned_tree,sv_pca_order_nozero[rownames(sv_pca_order_nozero)%in%pruned_tree$tip.label,])
homo_df = data.frame(names(homoplasic_svs))
colnames(homo_df) = "ID"
homo_df = homo_df%>% rowwise() %>% mutate(Position=as.numeric(strsplit(ID, "\\.")[[1]][1]),
                             Length=as.numeric(str_extract(ID, "(?<=\\.)\\d+$")),
                             Type=str_extract(ID, "(INS|DEL|DUP|INV)")) %>%
  ungroup()
homo_dens = homo_df %>% ggplot(.,aes(x=Position)) + geom_bar(fill="lightblue") +theme_pubr()
homo_df$END <- ifelse(homo_df$Type %in% c("INS", "DUP"),
                      homo_df$Position,
                      homo_df$Position + homo_df$Length)
homo_regions = GRanges(seqnames = "NC_000962.3", ranges = IRanges(start = homo_df$Position, end = homo_df$END))
overlaps=findOverlaps(homo_regions,gr_annotations)
overlapping_genes <- mcols(gr_annotations)[subjectHits(overlaps), "gene"]
table(overlapping_genes)[order(table(overlapping_genes))]

tree_data <- pruned_tree %>% 
  as_tibble() %>% 
  mutate(color = ifelse(label %in% rownames(sv_pca[sv_pca$"4053018.DUP.111"==1,]), "SV", "No-SV"))
ggtree(pruned_tree,layout = "circular") + 
  geom_tippoint(aes(color=tree_data$color)) +
  scale_color_manual(values = c("gray","#2171B5"))


tree_data <- pruned_tree %>% 
  as_tibble() %>% 
  mutate(color = ifelse(label %in% rownames(sv_pca[sv_pca$"2025842.DEL.3125"==1,]), "SV", "No-SV"))
ggtree(pruned_tree,layout = "circular") + 
  geom_tippoint(aes(color=tree_data$color)) +
  scale_color_manual(values = c("gray","#2171B5"),name='ppe25-ppe27 DEL') +
  geom_strip('GCA_030866965.1', 'ERR12743230', barsize=2, color='#CC6677', 
             label="L4.4", offset.text=.004,fontsize = 5)

#espI SVs
tree_data <- pruned_tree %>% 
  as_tibble() %>% 
  mutate(color = case_when(label %in% rownames(sv_pca[sv_pca$"4353406.DUP.51"==1,]) ~ "4353406.DUP.51",
                           label %in% rownames(sv_pca[sv_pca$"4353305.DEL.102"==1,]) ~ "4353305.DEL.102",
                           label %in% rownames(sv_pca[sv_pca$"4353314.DEL.51"==1,]) ~ "4353314.DEL.51",
                           label %in% rownames(sv_pca[sv_pca$"4353197.DEL.199"==1,]) ~ "4353197.DEL.199",
                           label %in% rownames(sv_pca[sv_pca$"4353197.DEL.199"==1,]) ~ "4353197.DEL.199",
                           label %in% rownames(sv_pca[sv_pca$"4353197.DEL.97"==1,]) ~ "4353197.DEL.97",
                           TRUE ~ "No-SV"))
ggtree(pruned_tree,layout = "circular") + 
  geom_tippoint(aes(color=tree_data$color)) +
  scale_color_manual(values = c("#2171B5","#D43925","lightblue2","green3","blue3","gray"))

#L1.2.1-only SVs
tree_data <- pruned_tree %>% 
  drop.tip(.,pruned_tree$tip.label[pruned_tree$tip.label %in% lin[lin$Lineage!="1.2.1.2" & lin$Lineage!="1.2.1.2.1" & lin$Lineage!="1.2.1",]$ID]) %>%
  #drop.tip(.,c("SRR12395199","SRR12395427","SRR12395416","SRR12801738")) %>%
  as_tibble() %>% 
  mutate(color = case_when(label %in% rownames(sv_pca[sv_pca$"1077957.DEL.2145"==1,]) ~ "1077957.DEL.2145",
                           label %in% rownames(sv_pca[sv_pca$"1078517.DEL.297"==1,]) ~ "1078517.DEL.297",
                           label %in% rownames(sv_pca[sv_pca$"1078521.DEL.1020"==1,]) ~ "1078521.DEL.1020",
                           label %in% rownames(sv_pca[sv_pca$"1078521.DEL.1189"==1,]) ~ "1078521.DEL.1189",
                           label %in% rownames(sv_pca[sv_pca$"1078521.DEL.1456"==1,]) ~ "1078521.DEL.1456",
                           label %in% rownames(sv_pca[sv_pca$"1078516.DEL.1586"==1,]) ~ "1078516.DEL.1586",
                           TRUE ~ "No-SV"))
ggtree(pruned_tree%>% drop.tip(.,pruned_tree$tip.label[pruned_tree$tip.label %in% lin[lin$Lineage!="1.2.1.2" & lin$Lineage!="1.2.1.2.1" & lin$Lineage!="1.2.1",]$ID]),
        layout = 'circular') + 
  geom_tippoint(aes(color=tree_data$color)) +
  scale_color_manual(values = c("#87CEFA", "#1E90FF", "#0ce7f2", "#4682B4", "#7B68EE", "#483D8B","gray"),name="ctpV DELs") +
  geom_strip('GCA_014899465.1', 'SRR12395061', barsize=2, color='#DDCC77', 
             label="L1.2.1", offset.text=.004,fontsize = 5)


#Look at number of SVs inside specific genes + essential genes and compare
essential = read.table("~/Downloads/METADATA.GENE.ESSENTIALITY.COMAS.2010.txt",sep=" ")
essential = essential %>% filter(V2 == "essential")
essentially = read_excel("~/Downloads/mbo002173137st3.xlsx",skip = 1)
essentially = essentially %>% select(`ORF ID`,Name,`Final Call`) %>% 
  filter(`Final Call` == "ES")
essentially = merge(essentially,essential,by.y=1,by.x=1) %>%
  select(-V2)
essentially$sv = 0
essentially$sv_pop = 0
essentially$lin = NA
essentially = rbind(essentially,genes_names %>% filter(!(Locus %in% essentially$`ORF ID`)) %>% select(Locus,Name) %>% mutate("Final Call"="NE","sv"=0,"sv_pop"=0,"lin"=NA) %>% dplyr::rename(`ORF ID` = Locus))
essentially = merge(essentially,genes_names[,c("V4","V5","Locus")],by.x=1,by.y=3)
for (i in 2:ncol(sv_pca)) {
  print(i)
  sv_start = sv_type[sv_type$POS<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])+25 & sv_type$POS>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])-25 & sv_type$SVTYPE==str_split(colnames(sv_pca)[i],"[.]")[[1]][2] & sv_type$SVLEN>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])-50 & sv_type$SVLEN<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])+50,]$POS
  sv_end <- ifelse(sv_type[sv_type$POS<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])+25 & sv_type$POS>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])-25 & sv_type$SVTYPE==str_split(colnames(sv_pca)[i],"[.]")[[1]][2] & sv_type$SVLEN>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])-50 & sv_type$SVLEN<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])+50,]$SVTYPE == "DEL",
                   sv_type[sv_type$POS<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])+25 & sv_type$POS>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])-25 & sv_type$SVTYPE==str_split(colnames(sv_pca)[i],"[.]")[[1]][2] & sv_type$SVLEN>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])-50 & sv_type$SVLEN<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])+50,]$POS + sv_type[sv_type$POS<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])+25 & sv_type$POS>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])-25 & sv_type$SVTYPE==str_split(colnames(sv_pca)[i],"[.]")[[1]][2] & sv_type$SVLEN>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])-50 & sv_type$SVLEN<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])+50,]$SVLEN,
                   sv_type[sv_type$POS<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])+25 & sv_type$POS>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][1])-25 & sv_type$SVTYPE==str_split(colnames(sv_pca)[i],"[.]")[[1]][2] & sv_type$SVLEN>as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])-50 & sv_type$SVLEN<as.numeric(str_split(colnames(sv_pca)[i],"[.]")[[1]][3])+50,]$POS)
  for (j in 1:nrow(essentially)) {
    if (j == "Rvnr01") {
      next
    }
    #print(sv_start)
    #print(sv_end)
    #print(i)
    #print(j)
    if (sv_start <= essentially[j,]$V5 & sv_end >= essentially[j,]$V4) {
      essentially[j,]$sv = essentially[j,]$sv + 1
      essentially[j,]$sv_pop = essentially[j,]$sv_pop + nrow(sv_pca[sv_pca[,i]==1,])
      if (is.na(essentially[j,]$lin)==TRUE) {
        essentially[j,]$lin = paste(unique(sv_pca[sv_pca[,i]==1,]$Lineage),collapse = ",")
      }
      else {
        essentially[j,]$lin = paste(unique(unlist(c(unique(sv_pca[sv_pca[,i]==1,]$Lineage),str_split(essentially[j,]$lin,",")))),collapse = ",")
      }
    }
    else if (sv_start < essentially[j,]$V4 & sv_end >= essentially[j,]$V4 & sv_end <= essentially[j,]$V5) {
      essentially[j,]$sv = essential[essentially$`ORF ID`==j,]$sv + 1
      essentially[j,]$sv_pop = essential[j,]$sv_pop + nrow(sv_pca[sv_pca[,i]==1,])
      if (is.na(essentially[j,]$lin)==TRUE) {
        essentially[j,]$lin = paste(unique(sv_pca[sv_pca[,i]==1,]$Lineage),collapse = ",")
      }
      else {
        essentially[j,]$lin = paste(unique(unlist(c(unique(sv_pca[sv_pca[,i]==1,]$Lineage),str_split(essentially[j,]$lin,",")))),collapse = ",")
      }
    }
  }  
}
#write.table(essentially,"~/Documents/PHD/SV/SV_analysis/temp/essentially_noillu.csv",sep=",")
#essentially2 = essentially
essentially = read_csv("~/Documents/PHD/SV/SV_analysis/temp/essentially.csv")
genes_to_plot = c("ES","NE","rpoB",   # Rifampicin resistance
  "katG",   # Isoniazid resistance
  "inhA",   # Isoniazid resistance (including promoter region mutations)
  "gyrA",   # Fluoroquinolone resistance (e.g., moxifloxacin, ofloxacin)
  "gyrB",   # Fluoroquinolone resistance
  "embB",   # Ethambutol resistance
  "pncA",   # Pyrazinamide resistance
  "ethA",   # Ethionamide resistance
  "rrs",    # Aminoglycoside resistance (e.g., streptomycin, kanamycin, amikacin)
  "tlyA",   # Capreomycin resistance
  "rpsL",   # Streptomycin resistance
  "gidB",   # Streptomycin resistance
  "eis",    # Kanamycin resistance
  "ndh",    # Isoniazid resistance
  "ahpC",   # Isoniazid resistance (compensatory mutations)
  "fabG1", 
  "espE","espF","espG1","espH","eccA1","eccB1","eccCa1","eccCb1","PE35","PPE68","esxA","esxB","espI","eccD1","espJ","espK","espL","espB","eccE1","mycP1", #ESX-1
  "eccA2","eccE2","mycP2","eccD2","Rv3888c","espG2","esxC","esxD","PPE69","PE36","eccC2","eccB2", #ESX-2
  "eccA3","eccB3","eccC3","Pe5","PPE4","esxG","esxH","espG3","eccD3","mycP3","eccE3", #ESX-3
  "esxT","esxU","Rv3446c","eccC4","eccD4","mycP4","eccB4", #ESX-4
  "eccB5", "eccC5", "cyp143", "Rv1786", "PPE25", "PE18", "PPE26", "PPE27",
  "PE19", "esxM", "esxN", "Rv1794", "eccD5", "mycP5", "eccE5", "eccA5",
  "Rv1034c", "Rv1035c", "Rv1036c", "esxI", "esxJ", "PPE15", "PE8", "Rv1041c",
  "Rv1042c", "Rv1194c", "PE13", "PPE18", "esxK", "esxL", "Rv1199c",
  "Rv3618", "esxV", "esxW", "PPE65", "PE32", "lpqG")
  
ecc_plot = esx_plot("ecc|mycP",genes_to_plot)

esp_plot = essentially %>% mutate(Plot = case_when(Name %in% genes_to_plot ~ Name,
                                        `Final Call` == "ES" ~ "ES",
                                        `Final Call` == "NE" ~ "NE"),
                       Color = case_when(Name %in% genes_to_plot ~ "Other",
                                         `Final Call` == "ES" ~ "ES",
                                         `Final Call` == "NE" ~ "NE"),
                       ESX = case_when(grepl("esp(E|F|G1|H|I|J|K|L|B)",Name) ~ "ESX1",
                                       grepl("espG2",Name) ~ "ESX2",
                                       grepl("espG3",Name) ~ "ESX3",
                                       TRUE ~ "Other")) %>%
  
  filter(Name %in% genes_to_plot[grep("esp",genes_to_plot)]) %>%
  ggplot(.,aes(x=factor(Plot,levels=genes_to_plot[grep("esp",genes_to_plot)]),y=sv_pop,color = ESX)) +
  geom_point(alpha=0.7,size=3) +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="ES",]$sv_pop,trim=0.1),color="#F08080",linetype="dashed") +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="NE",]$sv_pop,trim=0.1),color="lightblue",linetype="dashed") +
  #geom_quasirandom() +
  scale_color_manual(values = c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(),legend.title = element_blank()) +
  #coord_cartesian(expand = FALSE) +
  ylab("Number of isolates with SVs") +
  scale_y_continuous(breaks = seq(0,200,by=10))

esx1_plot = essentially %>% mutate(Plot = case_when(Name %in% genes_to_plot ~ Name,
                                                   `Final Call` == "ES" ~ "ES",
                                                   `Final Call` == "NE" ~ "NE"),
                                  Color = case_when(Name %in% genes_to_plot ~ "Other",
                                                    `Final Call` == "ES" ~ "ES",
                                                    `Final Call` == "NE" ~ "NE"),
                                  ESX = case_when(grepl("esx(A|B)",Name) ~ "ESX1",
                                                  grepl("esx(C|D)",Name) ~ "ESX2",
                                                  grepl("esx(G|H)",Name) ~ "ESX3",
                                                  grepl("esx(T|U)",Name) ~ "ESX4",
                                                  grepl("esx(M|N)",Name) ~ "ESX5",
                                                  grepl("esx(I|J)",Name) ~ "ESX5a",
                                                  grepl("esx(K|L)",Name) ~ "ESX5b",
                                                  grepl("esx(V|W)",Name) ~ "ESX5c",
                                                  TRUE ~ "Other"),
                                  sv_pop = replace(sv_pop,grepl("esxA|esxB|eccCb1|eccD1",Name), 0)) %>%
  filter(Name %in% genes_to_plot[grep("esx",genes_to_plot)]) %>%
  ggplot(.,aes(x=factor(Plot,levels=genes_to_plot[grep("esx",genes_to_plot)]),y=sv_pop,color = ESX)) +
  geom_point(alpha=0.7,size=3) +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="ES",]$sv_pop,trim=0.1),color="#F08080",linetype="dashed") +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="NE",]$sv_pop,trim=0.1),color="lightblue",linetype="dashed") +
  #geom_quasirandom() +
  scale_color_manual(values = c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#6A0DAD","#7B68EE","#A55194")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(),legend.title = element_blank()) +
  #coord_cartesian(expand = FALSE) +
  ylab("Number of isolates with SVs")

pe_plot = essentially %>% mutate(Plot = case_when(Name %in% genes_to_plot ~ Name,
                                                   `Final Call` == "ES" ~ "ES",
                                                   `Final Call` == "NE" ~ "NE"),
                                  Color = case_when(Name %in% genes_to_plot ~ "Other",
                                                    `Final Call` == "ES" ~ "ES",
                                                    `Final Call` == "NE" ~ "NE"),
                                  ESX = case_when(grepl("PE35|PPE68",Name) ~ "ESX1",
                                                  grepl("PE36|PPE69",Name) ~ "ESX2",
                                                  grepl("PPE4",Name) ~ "ESX3",
                                                  grepl("PPE18|PE13",Name) ~ "ESX5b",
                                                  grepl("PPE25|PPE26|PPE27|PE18|PE19",Name) ~ "ESX5",
                                                  grepl("PPE15|PE8",Name) ~ "ESX5a",
                                                  grepl("PPE65|PE32",Name) ~ "ESX5c",
                                                  TRUE ~ "Other")) %>%
  filter(Name %in% genes_to_plot[grep("PE|PPE",genes_to_plot)]) %>%
  ggplot(.,aes(x=factor(Plot,levels=genes_to_plot[grep("PE|PPE",genes_to_plot)]),y=sv_pop,color = ESX)) +
  geom_point(alpha=0.7,size=3) +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="ES",]$sv_pop,trim=0.1),color="#F08080",linetype="dashed") +
  geom_hline(yintercept = mean(essentially[essentially$`Final Call`=="NE",]$sv_pop,trim=0.1),color="lightblue",linetype="dashed") +
  #geom_quasirandom() +
  scale_color_manual(values = c("#1F77B4","#FF7F0E","#2CA02C","#9467BD","#6A0DAD","#7B68EE","#A55194")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(),legend.title = element_blank()) +
  #coord_cartesian(expand = FALSE) +
  ylab("Number of isolates with SVs")
  
pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/esx_conservation.pdf",width = 10, height = 8)
plot_grid(ecc_plot,esp_plot,esx1_plot,pe_plot,nrow=2,ncol=2,labels = c("a","b","c","d"))
dev.off()


#Homoplasic SVs in genes, not lineage-specific and not due to unmappable region
marin = read.table("~/Documents/PHD/SV/SV_analysis/DR/RLC_Regions.Plus.LowPmapK50E4.H37Rv.txt",sep="\t",col.names = c("start","end","length"))
marin = marin[,c(1,2)]
essentially_pos = essentially
for (i in 1:nrow(marin)) {
  essentially_pos = essentially_pos %>% filter(V4-3000>marin[i,2] &V5+3000>marin[i,2] | V4+3000<marin[i,1] &V5-3000<marin[i,1])
}
essentially_pos = essentially_pos[essentially_pos$sv_pop>9 & grepl(",",essentially_pos$lin),]
bin_counts=table(cut(essentially_pos$V4,
          breaks = seq(0, max(essentially_pos$V4) + 50000, by = 50000),
          labels = seq(0, max(essentially_pos$V4), by = 50000),
          include.lowest = TRUE))
bin_df <- as.data.frame(bin_counts)
colnames(bin_df) <- c("Bin", "Frequency")

png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/homo_fig.png",width = 700,height = 450)
ggplot(bin_df, aes(x = Bin, y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Chromosome Position (Bin Start)",
       y = "Frequency of genes with homoplasic SVs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#HOMOPLASIC SVs mapped to TREE
tree_data <- pruned_tree %>% 
  as_tibble() %>% 
  mutate(Rv3177 = case_when(label %in% rownames(sv_pca[sv_pca$"3544187.DEL.1329"==1,])~ "SV in Rv3177", 
                           label %in% rownames(sv_pca[sv_pca$"3544276.DEL.1330"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3544651.DEL.1074"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3545556.DEL.1327"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3545945.DEL.4582"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.INS.51"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.63"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.126"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.209"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.INS.627"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.731"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.1327"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.2027"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.2809"==1,])~ "SV in Rv3177",
                           label %in% rownames(sv_pca[sv_pca$"3546266.DEL.3185"==1,])~ "SV in Rv3177",
                           TRUE~" "),
         ctpG = case_when(label %in% rownames(sv_pca[sv_pca$"2234963.DEL.2571"==1,])~ "SV in ctpG", 
                            label %in% rownames(sv_pca[sv_pca$"2235169.DEL.63"==1,])~ "SV in ctpG",
                            label %in% rownames(sv_pca[sv_pca$"2235169.DEL.2048"==1,])~ "SV in ctpG",
                            label %in% rownames(sv_pca[sv_pca$"2235810.DEL.3313"==1,])~ "SV in ctpG",
                            label %in% rownames(sv_pca[sv_pca$"2236115.DEL.67"==1,])~ "SV in ctpG",
                            label %in% rownames(sv_pca[sv_pca$"2237052.DEL.3655"==1,])~ "SV in ctpG",
                            TRUE~" "),
         mmpS1 = case_when(label %in% rownames(sv_pca[sv_pca$"483248.INS.1257"==1,])~ "SV in mmpS1", 
                          label %in% rownames(sv_pca[sv_pca$"483248.INS.1359"==1,])~ "SV in mmpS1",
                          label %in% rownames(sv_pca[sv_pca$"483533.INS.1358"==1,])~ "SV in mmpS1",
                          TRUE~" "),
         glnA3 = case_when(label %in% rownames(sv_pca[sv_pca$"2127066.DEL.975"==1,])~ "SV in glnA3", 
                           label %in% rownames(sv_pca[sv_pca$"2127982.DEL.993"==1,])~ "SV in glnA3",
                           label %in% rownames(sv_pca[sv_pca$"2128379.DEL.1204"==1,])~ "SV in glnA3",
                           label %in% rownames(sv_pca[sv_pca$"2125462.DEL.3172"==1,])~ "SV in glnA3",
                           TRUE~" "))

homotree_fig = ggtree(pruned_tree, layout = "circular") %<+% lin + geom_tippoint(aes(colour=L1))+ scale_color_manual(values = lin_color)+ 
  geom_fruit(data=tree_data[!is.na(tree_data$label),],aes(y = label, fill=Rv3177),height=0.5,width = 0.007,geom='geom_tile') +
  geom_fruit(data=tree_data[!is.na(tree_data$label),],aes(y = label, fill=ctpG),height=0.5,width = 0.007,geom='geom_tile') +
  geom_fruit(data=tree_data[!is.na(tree_data$label),],aes(y = label, fill=mmpS1),height=0.5,width = 0.007,geom='geom_tile') +
  geom_fruit(data=tree_data[!is.na(tree_data$label),],aes(y = label, fill=glnA3),height=0.5,width = 0.007,geom='geom_tile') +
scale_fill_manual(values = c("SV in mmpS1"="#2CA02C","SV in ctpG"="#FF7F00","SV in Rv3177"="#2171B5"," "="white","SV in glnA3" = "#E31A1C")) + theme(legend.title = element_blank())

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/homotree_fig.pdf",height = 10,width = 10)
homotree_fig
dev.off()


#Plot Manhattan plot!!
pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/rnaseqa.pdf",width = 8,height = 5)
ggplot(essentially %>% 
         mutate(Sel = case_when(lin == "1" | lin=="2"|lin=="3"|lin=="4"|lin=="5"|lin=="6"|lin=="7"|lin=="8"|lin=="9" & sv > 1 ~ "Lineage-specific positive selection", 
                                TRUE ~ "No"))%>% filter(Sel=="Lineage-specific positive selection"), aes(x = V4, y = sv, color = Sel)) +
  geom_point(aes(fill = Sel),size=1.75,shape=21,color="black") +
  geom_text(data = . %>% filter(Name %in% c("ctpV")), 
            aes(label = Name), vjust = -1, hjust = 0.5,color="black", size=6) +
  geom_text(data = . %>% filter(Name %in% c("ethA")), 
            aes(label = Name), vjust = -1, hjust = 0.5,color="black", size=6) +
  geom_text(data = . %>% filter(Name %in% c("Rv2765")), 
            aes(label = Name), vjust = -1, hjust = 0.5,color="black", size=6) +
  geom_text(data = . %>% filter(Name %in% c("Rv2337c")), 
            aes(label = Name), vjust = -1, hjust = 0.5,color="black", size=6) +
  geom_text(data = . %>% filter(Name %in% c("thyA")), 
            aes(label = Name), vjust = -1, hjust = 0.5,color="black", size=6) +
  scale_color_manual(values = c("Lineage-specific positive selection" = "#DDCC77")) +
  theme_pubr() +
  ylim(2,7) +
  labs(x="Genome position",y="Number of unique SVs per gene") +
  theme(legend.title = element_blank(),legend.text = element_text(size = 14))
dev.off()

#PLOT FIGURE 2
panel1=plot_grid(Fig2A,Fig2B,labels = c('a','b'),nrow = 1,ncol = 2,label_size = 40)

panel2 = plot_grid(panel1,plot_list(y),labels = c('','c'),nrow = 2,ncol = 1,rel_heights = c(1,2),label_size = 40)

Fig2 = plot_grid(panel2,Fig2D,labels = c('','d'),nrow = 2,ncol = 1,rel_heights = c(2,1),label_size = 40)

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/Figure1.pdf",height = 60,width = 38)
Fig2
dev.off()


###
#GGGENES FIGURE - Visualization using the PRG
###
# ESX5
esx5 <- data.frame(
  molecule = c(rep("PPE25;PE18;PPE26;PPE27 (H37Rv)",4), "PPE27", rep("PPE25;(PE18;PPE26;PPE27)x2",7),rep("PPE25;PE18-INS;PPE26;PPE27",6),rep("DEL-PPE26;PPE27",2),rep("DEL-PPE25;PE18;PPE26;DEL-PPE27",4)),
  gene = c("PPE25", "PE18", "PPE26", "PPE27", "PPE27","PPE25", "PE18", "PPE26", "PPE27", "PE18", "PPE26", "PPE27","PPE25","PE18","INS","PE18", "PPE26", "PPE27","PPE26","PPE27","PPE25","PE18", "PPE26", "PPE27"),
  start = c(1, 1177, 1490, 3125, 1, 1, 1177, 1490, 3125,4491, 4804, 6439,1,1177,1363,2719,3139,4774,1,874,1,422, 735, 2370),
  end = c(1098, 1476, 2671, 4177, 1052, 1098, 1476, 2671, 4177,4790, 5985, 7491,1098,1362,2718,2826,4320,5826,420,1926,343,721, 1916, 2878)
)
cb_palette <- c("gray","#D55E00", "#E69F00", "#56B4E9", "#009E73")
esx5$molecule = factor(esx5$molecule, levels=c("DEL-PPE26;PPE27","DEL-PPE25;PE18;PPE26;DEL-PPE27","DEL-PPE25;PE18;PPE26;PPE27","PPE25;(PE18;PPE26;PPE27)x2","PPE25;PE18-INS;PPE26;PPE27","PPE27","PPE25;PE18;PPE26;PPE27 (H37Rv)"))
# Create the gene plot
esx5_graph <- ggplot(esx5, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height =grid::unit(4, "mm"),  arrowhead_width = grid::unit(1, "mm")) +
  theme_genes() +
  scale_fill_manual(values = cb_palette) +
  theme(legend.title = element_blank(), axis.title.y = element_blank(),axis.line.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_blank(), legend.text = element_text(size=16)) +
  scale_y_discrete(expand = c(1.4,0))

#
# ESX5C
esx5c <- data.frame(
  molecule = c(rep("Rv3618;esxV;esxW;PPE65;PE32;lpqG (H37Rv)",6),rep("DEL-lpqG",1), rep("Rv3618;esxV;esxW;DEL-PPE65;DEL-PE32;lpqG",6),rep("Rv3618;esxV;esxW;PPE65;lpqG",5)),
  gene = c("Rv3618", "esxV", "esxW", "PPE65", "PE32","lpqG","lpqG", "Rv3618", "esxV", "esxW", "PPE65", "PE32","lpqG","Rv3618","esxV","esxW","PPE65","lpqG"),
  start = c(1, 1572, 1896, 3195, 3495, 3824,1,1, 1572, 1896, 2497, 2623, 2952,1, 1572, 1896, 3195, 3236),
  end = c(1188, 1287, 1599, 1953, 3195, 4547,516, 1188, 1287, 1599, 1953, 2497, 3675,1188, 1287, 1599, 1953, 3959)
)
cbc_palette <- c("#007943","#00441B","#D7191C","#66C2A5","#1F78B4","#CC79A7")
esx5c$molecule = factor(esx5c$molecule,levels=c("Rv3618;esxV;esxW;PPE65;lpqG","Rv3618;esxV;esxW;DEL-PPE65;DEL-PE32;lpqG","DEL-lpqG","Rv3618;esxV;esxW;PPE65;PE32;lpqG (H37Rv)"))
# Create the gene plot
esx5c_graph <- ggplot(esx5c, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height =grid::unit(4, "mm"),  arrowhead_width = grid::unit(1, "mm")) +
  theme_genes() +
  scale_fill_manual(values = cbc_palette) +
  theme(legend.title = element_blank(), axis.title.y = element_blank(),axis.line.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_blank(), legend.text = element_text(size=16)) +
  scale_y_discrete(expand = c(1.5,0))

# ESX5B
esx5b <- data.frame(
  molecule = c(rep("Rv1194c;PE13;PPE18;esxK;esxL;Rv1199c (H37Rv)",6), rep("Rv1194c;PE13;PPE18-DEL;Rv1199c",4),rep("Rv1194c;Rv1199c",2)),
  gene = c("Rv1194c","PE13","PPE18","esxK","esxL","Rv1199c","Rv1194c","PE13","PPE18","Rv1199c","Rv1194c","Rv1199c"),
  start = c(1267, 1757, 2104, 3415, 3763, 5364,1267, 1757, 2104, 4281,1236,2484),
  end = c(1, 2057, 3280, 3712, 4048, 4116,1, 2057, 3023, 3033,1,1236)
)
cbb_palette <- c("#1ABC9C","#FF6F61","#B39DDB","#9ACD32","#00CED1","#FF19AF")
esx5b$molecule = factor(esx5b$molecule,levels=c("Rv1194c;Rv1199c","Rv1194c;PE13;PPE18-DEL;Rv1199c","Rv1194c;PE13;PPE18;esxK;esxL;Rv1199c (H37Rv)"))
# Create the gene plot
esx5b_graph <- ggplot(esx5b, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height =grid::unit(4, "mm"),  arrowhead_width = grid::unit(1, "mm")) +
  theme_genes() +
  scale_fill_manual(values = cbb_palette) +
  theme(legend.title = element_blank(), axis.title.y = element_blank(),axis.line.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_blank(), legend.text = element_text(size=16)) +
  scale_y_discrete(expand = c(1.5,0))
#

# ESX5A
esx5a <- data.frame(
  molecule = c(rep("esxI;esxJ;PPE15;PE8 (H37Rv)",4), rep("esxI;esxJ;INS;PPE15;PE8",5),rep("esxI;INS;esxJ;PPE15;PE8",5),rep("esxI;esxJ-INS;PPE15;PE8",6)),
  gene = c("esxI","esxJ","PPE15","PE8","esxI","esxJ","INS","PPE15","PE8","esxI","INS","esxJ","PPE15","PE8","esxI","esxJ","INS","esxJ","PPE15","PE8"),
  start = c(285, 609, 1931, 2836, 285, 609,2128, 3315, 4212, 285, 1640, 1963, 3285, 4190, 285, 560,1919,1969, 3290, 4195),
  end = c(1, 312, 755, 2008, 1, 312, 752,2131, 3384,1, 286, 1666, 2109, 3362, 1, 312,561,1920, 2114, 3367)
)
cba_palette <- c("#8E44AD","#FF4500","#D3D3D3","#C0392B","#48C9B0","#F39C12","#3498DB")
esx5a$molecule = factor(esx5a$molecule,levels=c("esxI;esxJ-INS;PPE15;PE8","esxI;INS;esxJ;PPE15;PE8","esxI;esxJ;INS;PPE15;PE8","esxI;esxJ;PPE15;PE8 (H37Rv)"))
# Create the gene plot
esx5a_graph <- ggplot(esx5a, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrow_body_height =grid::unit(4, "mm"),  arrowhead_width = grid::unit(1, "mm")) +
  theme_genes() +
  scale_fill_manual(values = cba_palette) +
  theme(legend.title = element_blank(), axis.title.y = element_blank(),axis.line.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.y = element_blank(), legend.text = element_text(size=16)) +
  scale_y_discrete(expand = c(1.5,0))
#

### GC COUNT
bkp_gc = read.table("~/Documents/PHD/SV/SV_analysis/temp/gc_bkp.txt")
random_gc = read.table("~/Documents/PHD/SV/SV_analysis/temp/gc_ran.txt")
bkp_gc$TYPE = "SV"
random_gc$TYPE = "Non-SV"
gc_plot = ggplot(rbind(bkp_gc,random_gc),aes(x=V1,fill=TYPE)) +
  geom_histogram(alpha=0.4,position = "identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_pubr() +
  xlab("GC%") +
  theme(legend.title = element_blank(),axis.text = element_text(size = 40),axis.title = element_text(size = 40),
        legend.text = element_text(size = 40),axis.ticks = element_line(linewidth = 3),legend.key.size = unit(2,"cm"))
png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/gc_plot.png",height = 1000,width = 1000)
gc_plot
dev.off()
shapiro.test(bkp_gc$V1)
shapiro.test(random_gc$V1)
wilcox.test(random_gc$V1,bkp_gc$V1)

#random gene overlap
# Function to generate random regions
generate_random_regions <- function(n, region_length, genome_size) {
  starts <- sample(1:(genome_size - region_length), size = n, replace = FALSE)
  GRanges(seqnames = "NC_000962.3", ranges = IRanges(start = starts, end = starts + region_length - 1))
}

# Permutation test
set.seed(123)  # For reproducibility
n_permutations <- 100
random_overlaps <- numeric(n_permutations)

for (i in 1:n_permutations) {
  random_regions <- generate_random_regions(100, 1, 4411532)
  random_overlaps[i] <- sum(countOverlaps(random_regions, gr_annotations))
}

# SV regions
sv_type2 = sv_type
sv_type2$END=sv_type$POS+100
sv_regions <- GRanges(seqnames = "NC_000962.3", ranges = IRanges(start = sv_type2$POS, end = sv_type2$POS))

# Observed overlaps
observed_overlap <- sum(countOverlaps(sv_regions, gr_annotations))

ggplot(data.frame(ro=random_overlaps),aes(x=ro)) + geom_histogram(bins = 15) +theme_pubr() +
  geom_vline(aes(xintercept = observed_overlap/34.05),color="red") + xlab("% of breakpoints that fall in genes")

z = (observed_overlap/34.05 -mean(random_overlaps))/sd(random_overlaps)
p_value = 2 * (1 - pnorm(abs(z)))
