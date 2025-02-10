#### Import packages ####
library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(vegan)
library(moments)
library(ggplot2)
library(rtracklayer)
library(mclust)
library(ggpubr)
library(sva)
library(ggstream)
library(GenomicRanges)
library(reshape2)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(msa)
library(pheatmap)


#### Functions ####
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
calculating_TMM_RPKM_allSamples <- function(counts_path, # path of htseq count files (.count)
                                            #sample_study, # data frame. Column "Run" is SRA accession (refer to htseq-count file name) and column "BioProject" is bio-projects
                                            gene_length, # data frame. Column "gene" is locus tag and column "length" is gene length
                                            filter_lib_size = 1000000, # a numeric that sample with library size <= 'filter_lib_size' will be excluded
                                            filter_geneLength = 150, # genes shorter than 150 bp will be removed
                                            filter_zeroExprGenes = TRUE, # whether remove genes with zero counts in all samples
                                            filter_gene = NA # a vector contains genes that are excluded
) {
  # set path
  setwd(counts_path)
  
  # read count files
  counts_data <- readDGE(list.files(pattern =".count"), columns=c(1,2))
  counts_data$counts <- counts_data$counts[1:(nrow(counts_data$counts)-5),]
  raw.counts <- as.data.frame(counts_data$counts)
  
  # filter library size
  #lib_size <- counts_data$samples
  #count_sample <- rownames(lib_size)[lib_size$lib.size > filter_lib_size] 
  #count_sample <- count_sample[count_sample %in% sample_study$Run]
  #raw.counts <- raw.counts[,colnames(raw.counts) %in% count_sample]
  
  # filter genes
  if (filter_zeroExprGenes) {raw.counts <- raw.counts[apply(raw.counts, 1, function(x) {sum(x) > 0}),]}
  gene_geneLength <- gene_length$gene[gene_length$length >= filter_geneLength]
  count_gene <- rownames(raw.counts)
  count_gene <- count_gene[count_gene %in% gene_geneLength & (!count_gene %in% filter_gene)]
  raw.counts <- raw.counts[rownames(raw.counts) %in% count_gene,]
  
  # re-shape data
  count_gene <- rownames(raw.counts)
  
  gene_length <- gene_length[gene_length$gene %in% count_gene,]
  gene_length <- as.data.frame(gene_length[match(gene_length$gene, count_gene),])
  gene_length <- data.frame(row.names = gene_length$gene, length = gene_length$length)
  
  #sample_study <- sample_study[sample_study$Run %in% count_sample,]
  #sample_project <- sample_study$BioProject[!duplicated(sample_study$BioProject)]
  
  # calculate TMM normalized RPKM 
  temp.count <- raw.counts
  
  temp.counts.data <- DGEList(counts = temp.count)
  temp.counts.data <- calcNormFactors(temp.counts.data)
  temp.counts.data$genes <- gene_length
  
  RPKM <- as.data.frame(rpkm(temp.counts.data))
  return(RPKM)
}


#### Input data ####

# gene length
gtf_Mtb <- read.table("~/Downloads/Mycobacterium_tuberculosis_H37Rv_gff_v4.1.gff",sep=" ",header = F)

gtf_Mtb$V9 = as.character(gtf_Mtb$V8)
# Apply the function to the data frame and convert to a wide format
genes_names <- gtf_Mtb %>%
  dplyr::mutate(parsed_info = purrr::map(V9, parse_info)) %>%
  unnest_wider(parsed_info) %>%
  filter(Name != "ncRv3520")

gr_annotations <- GRanges(seqnames = Rle("NC_000962.3"), 
                          ranges = IRanges(start = genes_names$V4, end = genes_names$V5),
                          gene = genes_names$Locus)

annotation_df <- as.data.frame(gr_annotations)
geneLength_Mtb <- annotation_df[,c(6,4)] %>% group_by(gene) %>% filter(n()==1) %>% ungroup()
colnames(geneLength_Mtb) = c("gene","length")


#### Calculate TMM normalized RPKM ####
# Mtb
filter_gene_mtb <- c("MTB000019","MTB000020","MTB000021") # remove ribosome RNA
filter_gene_mtb <- unique(genes_names[grepl("nc",genes_names$Name),]$Locus) # remove noncoding RNA

rpkm_Mtb <- calculating_TMM_RPKM_allSamples(counts_path = "~/Documents/PHD/SV/SV_analysis/Manila/rna_counts_unstrand", gene_length = geneLength_Mtb,
                                            filter_gene = filter_gene_mtb)
#setwd("../../../")
#write.csv(rpkm_Mtb, "./data/unfiltered_tmm_rpkm_Mtb.csv")
rpkm_Mtb <- rpkm_Mtb[rowSums(rpkm_Mtb == 0) < 2, ]
#Now lets take a look at the batch effects
manila = c(1,1,2,1,1,2,1,1,1,1,1,1,1,2,2,2,2)
names(manila) = colnames(rpkm_Mtb)
plotMDS(rpkm_Mtb,col=manila)
batch=c(1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,3,4)
names(batch)=colnames(rpkm_Mtb)
plotMDS(rpkm_Mtb,col=batch)

#We can see that there is a clear batch effect among samples from different studies. We will only keep those from ERR34*
plotMDS(rpkm_Mtb[,-c(7,8,9,10,11,12,13,14,15,16,17)],col=manila[c(1:6)])
rpkm_Mtb_nobatch = rpkm_Mtb[,-c(7,8,9,10,11,12,13,14,15,16,17)]

#Since RPKM values can vary widely, it’s important to log-transform the data to make the distribution more appropriate for analysis
rpkm_Mtb_nobatch_log = log2(rpkm_Mtb_nobatch +1)
sample_info <- data.frame(
  Sample = colnames(rpkm_Mtb_nobatch_log),
  Condition = factor(c("NonL1.2.1","NonL1.2.1","L1.2.1","NonL1.2.1","NonL1.2.1","L1.2.1"))  # Modify based on your data
)
rownames(sample_info) <- sample_info$Sample
# Create design matrix based on the conditions
design <- model.matrix(~0 + sample_info$Condition)
colnames(design) <- levels(sample_info$Condition)

# Fit the linear model
fit <- lmFit(rpkm_Mtb_nobatch_log, design)

# Create contrast matrix to specify which comparisons you want to make
contrast_matrix <- makeContrasts(L1.2.1_vs_NonL1.2.1 = L1.2.1 - NonL1.2.1, levels=design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes smoothing to the standard errors
fit2 <- eBayes(fit2)

#We extract significant differentially expressed genes
significant_genes <- topTable(fit2, adjust.method = "fdr", number = Inf)
results <- topTable(fit2, adjust="fdr", number=Inf)
significant_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
significant_genes=significant_genes[!rownames(significant_genes) %in% c("Rv2645","Rv2646c","Rv2647c","Rv2648c","Rv2649c","Rv2650c","Rv2651c","Rv2652c","Rv2653c","Rv2654c","Rv2655c",
                    "Rv2656c","Rv2657c","Rv2658c","Rv2659c","Rv1967","Rv1968","Rv1971"),]
results=results[!rownames(results) %in% c("Rv2645","Rv2646","Rv2647","Rv2648c","Rv2649c","Rv2650c","Rv2651c","Rv2652c","Rv2653c","Rv2654c","Rv2655c",
                                                                        "Rv2656c","Rv2657c","Rv2658c","Rv2659c","Rv1966","Rv1965","Rv1969",
                                          "Rv1973","Rv1974","Rv1975","Rv1976c"),]
logFC_cutoff <- 1
fdr_cutoff <- 0.05
pvalue_cutoff <- 0.05

# Add a column to classify genes based on significance and directionality
results$Category <- with(results, 
                                   ifelse(adj.P.Val < fdr_cutoff & logFC > logFC_cutoff, "Significantly Enriched",
                                          ifelse(adj.P.Val < fdr_cutoff & logFC < -logFC_cutoff, "Significantly Depleted",
                                                 #ifelse(P.Value < pvalue_cutoff & logFC > logFC_cutoff, "Moderately Enriched",
                                                        #ifelse(P.Value < pvalue_cutoff & logFC < -logFC_cutoff, "Moderately Depleted",
                                                               "Not Significant")))

# Set up custom color palette
color_palette <- c("Significantly Enriched" = "#D43925", 
                   "Significantly Depleted" = "#2376B2", 
                   #"Moderately Enriched" = "lightcoral", 
                   #"Moderately Depleted" = "lightblue", 
                   "Not Significant" = "grey")
results=merge(results,genes_names,by.x="row.names",by.y=10)

# Basic Volcano Plot
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  
  # Plot the points, colored by category
  geom_point(aes(color = Category), size = 3) +
  
  # Customize the color scale
  scale_color_manual(values = color_palette) +
  
  # Add labels to significant genes
  geom_text_repel(data = subset(results, Category %in% c("Significantly Enriched", "Significantly Depleted")), 
                  aes(label = Name), 
                  size = 4, max.overlaps = 10) +
  
  # Add custom lines for p-value and adjusted p-value cutoffs
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "black", size = 1.2) +
  #geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dotted", color = "black", size = 1, label="P-value = 0.05") +
  
  # Add vertical lines for fold change cutoffs
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", size = 1) +
  
  # Customize axis labels and title
  labs(x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  
  theme_minimal(base_size = 25) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(colour = "black", fill = NA, size = 1),legend.title = element_blank(),
        plot.margin =  unit(c(0, 1, 0, 1), "cm")) +
  
  # Set square aspect ratio
  coord_fixed()

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/rnaseqd.pdf",width = 12,height = 6)
volcano_plot
dev.off()

functions = read.csv("~/Documents/PHD/SV/SV_analysis/Manila/Mtb_gene_function.csv")
results = merge(results,functions,by.x=18,by.y=1)

results_significant = results[results$adj.P.Val<0.05 & abs(results$logFC)>1,]
results_significant = results_significant%>%mutate(sign=ifelse(logFC<0,"-","+"))
resfun = results_significant%>%group_by(Functions,sign)%>%tally()

# Summarize counts by Category and Sign
df_summary <- aggregate(n ~ Functions + sign, data = resfun, sum)

# Reshape the data to wide format
df_wide <- reshape(df_summary, idvar = "Functions", timevar = "sign", direction = "wide")

# Replace NA with 0
df_wide[is.na(df_wide)] <- 0

# Rename columns
colnames(df_wide) <- c("Category", "Negative", "Positive")

# Melting data for plotting
df_melt <- reshape2::melt(df_wide, id.vars = "Category", variable.name = "Sign", value.name = "Count")

# Adjust signs for plotting
df_melt$Count <- ifelse(df_melt$Sign == "Negative", -df_melt$Count, df_melt$Count)

# Modify the 'Category' column
df_melt$Category <- gsub("virulence detoxification adaptation", 
                    "Virulence, detoxification and adaptation", 
                    df_melt$Category)

df_melt$Category <- gsub("\\b(\\w)", "\\U\\1", df_melt$Category, perl = TRUE) # Capitalize first letter of each word

# Plotting
DEG_FUNC = ggplot(df_melt, aes(x = Category, y = Count, fill = Sign)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("Positive" = "#D43925", "Negative" = "#2376B2")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none") +
  labs(title = "Differentially Expressed Genes by Category",
       x = "", y = "Gene Count") +
  coord_flip()  # Flips the axis for better readability

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/rnaseqe.pdf",width = 8,height = 3)
DEG_FUNC
dev.off()



genes_of_interest <- c("Rv0968", "ctpV","ctpG","ctpB")  # replace with actual gene names
df_genes <- rpkm_Mtb_nobatch_log %>% rownames_to_column(var = "Gene") %>% merge(.,genes_names[,c(10,11)],by.x=1,by.y=1) %>%
  filter(Name %in% genes_of_interest)

manila_samples <- c("ERR1140772", "ERR1140791")  # Example Manila samples

# Reshape the dataframe (long format)
df_long <- df_genes %>%
  pivot_longer(cols = -c(Gene,Name), names_to = "Sample", values_to = "Expression")

# Add lineage information
df_long <- df_long %>%
  mutate(Lineage = ifelse(Sample %in% manila_samples, "L1.2.1", "Non-L1.2.1"),
         Colour = case_when(Sample %in% manila_samples & Gene == "ctpB" ~ "Dec",
                            Sample %in% manila_samples ~ "Inc",
                            TRUE ~ "No"))

# Plot the violin plot
boxplots = ggplot(df_long, aes(x = Lineage, y = Expression)) +
  #geom_boxplot(color="black",outlier.shape = NA) +
  #geom_jitter(width = 0.2, size = 1.5, aes(color = Colour), alpha = 0.6) +
  geom_point(aes(fill=Colour),shape=21,color="black",size=5) +
  facet_wrap(~ Name, scales = "free_y", nrow = 1, strip.position = "top") +  # Facet by gene
  theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.spacing = unit(0.2, "lines"), 
               panel.background=element_rect(fill="white"),
               panel.border = element_blank(),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.ticks = element_line(),
               axis.ticks.length = unit(.25, "cm"),
               axis.line = element_line(size = 0.5),
               axis.text = element_text(size = rel(1), color = 'black'),
               legend.key = element_blank(),
               legend.position = "none",
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="white",fill="white"),
               strip.text = element_text(face="bold"),
               axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  #theme_pubr() +
  labs(title = "",
       x = "", y = "log2(RPKM)") +
  scale_fill_manual(values = c("Dec" = "#2376B2","Inc" = "#D43925", "No" = "black"))

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/rnaseqf.pdf",width = 12,height = 4.2)
boxplots
dev.off()
