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
library(patchwork)
library(factoextra)
library(rlist)
library(grid)
library(purrr)
library(gridExtra)
library(magick)
library(forcats)
library(treeio)
library(plotly)
library(scales)
library(ggnewscale)
library(scatterpie)
library(VariantAnnotation)
library(parallel)
library(caret)
library(randomForest)
library(xgboost)
library(cowplot)
library(arm)
library(kableExtra)
library(formattable)
library(webshot2)
library(htmltools)
library(htmlwidgets)
library(vroom)
library(irlba)
library(logistf)
#
##FUNCTIONS
#
# Function to check overlap considering position, type, and size (for DELs)
check_overlap <- function(merged_sv, sample_sv) {
  pos_overlap <- (sample_sv$POS >= (merged_sv$POS - 50)) & (sample_sv$POS <= (merged_sv$POS + 50))
  
  type_match <- (grepl(strsplit(as.character(sample_sv$ID),".",fixed = TRUE)[[1]][2], strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][2])  | grepl("INS", strsplit(as.character(sample_sv$ID),".",fixed = TRUE)[[1]][2]) & grepl("DUP", strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][2]) |
    (grepl("DUP", strsplit(as.character(sample_sv$ID),".",fixed = TRUE)[[1]][2]) & grepl("INS", strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][2])))
  
  if (grepl("DEL", merged_sv$ID) & grepl("DEL", sample_sv$ID)) {
    merged_size <- as.numeric(strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][3])
    size_diff <- abs(merged_size - sample_sv$SIZE)
    print(pos_overlap)
    print(type_match)
    print(size_diff)
    return(any(pos_overlap & type_match & (size_diff <= 100)))
  } else {
    return(any(pos_overlap & type_match))
  }
}
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
# Function to process a single VCF file
process_mtb_vcf <- function(vcf) {
  # Extract relevant columns
  sv_data = data.frame(POS = "",
                       END = "",
                       ID = "",
                       SIZE = "",
                       stringsAsFactors = FALSE)
  for (m in 1:nrow(vcf)) {
    if (is.na(strsplit(as.character(vcf[m,3]),".",fixed = TRUE)[[1]][2])) {
    id_col <- NA
    }
    else {
      id_col <- paste(vcf[m,2],vcf[m,3],sep = ".")
    }
    pos_col <- vcf[m,2]
    end_col <- vcf[m,2]
    
    id_col[is.na(id_col)] <- "NaN"
    
    # Extract the size of DELs from the ID column
    sv_size <- ifelse(grepl("DEL", id_col), strsplit(as.character(vcf[m,3]),".",fixed = TRUE)[[1]][2], NA)
    
    sv_next <- data.frame(
      POS = pos_col,
      END = end_col,
      ID = id_col,
      SIZE = sv_size,
      stringsAsFactors = FALSE
    )
    sv_data = rbind(sv_data,sv_next)
  }
  return(sv_data)
}

# Function to merge SVs based on specified rules
merge_svs <- function(sv_data) {
  sv_data$POS = as.numeric(as.character(sv_data$POS))
  sv_data$END = as.numeric(as.character(sv_data$END))
  sv_data$SIZE = as.numeric(as.character(sv_data$SIZE))
  sv_data <- sv_data[order(sv_data$POS,sv_data$ID), ]  # Ensure data is sorted by position
  
  merged_svs <- list()
  current_sv <- sv_data[1, ]
  
  for (i in 2:nrow(sv_data)) {
    print(i)
    next_sv <- sv_data[i, ]
    
    same_type <- grepl(strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][2], strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][2])  # Match SV types (INS with INS, DEL with DEL, etc.)
    is_del <- grepl("DEL",strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][2]) & grepl("DEL", strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][2])
    is_ins_dup <- (grepl("INS", strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][2]) & grepl("DUP", strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][2])) |
      (grepl("DUP", strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][2]) & grepl("INS", strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][2]))
    
    # Check merging conditions
    if (is_del) {
      if (abs(current_sv$END - next_sv$POS) <= 50 &&
          abs(as.numeric(strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][3]) - as.numeric(strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][3])) <= 100) {
        current_sv$END <- max(current_sv$END, next_sv$END)
      } else {
        merged_svs <- rbind(merged_svs, current_sv)
        current_sv <- next_sv
      }
    }
    else {
      if (same_type && abs(current_sv$END - next_sv$POS) <= 50) {
        current_sv$END <- max(current_sv$END, next_sv$END)
      } else if (is_ins_dup && abs(current_sv$END - next_sv$POS) <= 50) {
        current_sv$END <- max(current_sv$END, next_sv$END)
      } else {
        merged_svs <- rbind(merged_svs, current_sv)
        current_sv <- next_sv
      }
    }
  }
  merged_svs <- rbind(merged_svs, current_sv)
  return(merged_svs)
}

process_sample <- function(i) {
  sv_data <- process_mtb_vcf(tabledr[[i]])
  sample_id <- sample_names[i]
  
  sv_data$POS <- as.numeric(as.character(sv_data$POS))
  sv_data$END <- as.numeric(as.character(sv_data$END))
  sv_data$SIZE <- as.numeric(as.character(sv_data$SIZE))
  
  for (j in seq_along(merged_svs$ID)) {
    merged_sv <- merged_svs[j, ]
    
    if (any(sv_data$ID == "NaN" & sv_data$POS <= merged_sv$END & sv_data$END >= merged_sv$POS)) {
      sv_matrix[sample_id, j] <- NA
    } else {
      overlap_result <- mapply(check_overlap, list(merged_sv), split(sv_data, seq(nrow(sv_data))))
      sv_matrix[sample_id, j] <- ifelse(any(overlap_result), 1, 0)
    }
  }
}

filter_columns_by_position <- function(df, range_df) {
  # Extract numeric positions from column names that match the pattern "X<number>.DEL"
  col_positions <- sapply(names(df), function(col_name) {
    # Extract numeric position, if the pattern matches, otherwise return NA
    match <- regmatches(col_name, regexpr("X[0-9]+\\.(INS|DUP|DEL)", col_name))
    if (length(match) > 0) {
      as.numeric(gsub("X([0-9]+)\\.(DEL|INS|DUP).*", "\\1", col_name))
    } else {
      NA  # Return NA for non-matching column names
    }
  })
  
  # Identify columns that match the pattern and columns that don't
  valid_positions <- !is.na(col_positions)
  non_matching_columns <- names(df)[!valid_positions]  # Keep non-matching columns
  
  # Check if position falls within any of the ranges in range_df for valid columns
  cols_to_keep <- sapply(col_positions[valid_positions], function(pos) {
    !any(pos >= range_df$start & pos <= range_df$end)
  })
  
  # Combine columns that either don't match the pattern or fall outside the range
  columns_to_keep <- c(non_matching_columns, names(df)[valid_positions][cols_to_keep])
  
  # Subset the dataframe to keep the desired columns
  df_filtered <- df[, columns_to_keep, drop = FALSE]
  
  return(df_filtered)
}

get_importance_matrix <- function(sv_matrix_pheno_lin, tbprofiler, drug="bedaquiline", nrounds=100) {
  
  # Step 1: Identify the column index of the selected drug in sv_matrix_pheno_lin
  drug_col_index <- which(colnames(sv_matrix_pheno_lin) == drug)
  
  # Step 2: Filter and mutate the data based on the drug column index
  sv_matrix_pheno_lin_mod <- sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin[[drug]]), c(drug_col_index, 23:ncol(sv_matrix_pheno_lin))] %>% 
    mutate_all(as.factor) # %>% dplyr::select(-all_of("X.1"))
  
  # Step 3: Remove columns where factors have fewer than 2 levels
  sv_matrix_final <- sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
  
  # Step 4: Add drug resistance prediction column
  sv_matrix_final_tb <- sv_matrix_final %>% 
    mutate(DR_prediction = ifelse(rownames(sv_matrix_final) %in% tbprofiler[tbprofiler$drugs == drug, ]$sample, 1, 0))
  
  # Step 5: Convert the data to numeric and adjust the drug column
  sv_matrix_final_tb_num <- sv_matrix_final_tb %>% 
    mutate_all(as.numeric) %>% 
    mutate(!!drug := ifelse(!!sym(drug) == 1, 1, 0))  # Adjust the drug column
  
  # Step 6: Convert to DMatrix for XGBoost
  dtrain <- xgb.DMatrix(
    data = as.matrix(sv_matrix_final_tb_num[, 2:ncol(sv_matrix_final_tb_num)]),
    label = sv_matrix_final_tb_num[[drug]]
  )
  
  # Step 7: Define XGBoost parameters
  params <- list(
    objective = "binary:logistic",  # Classification objective
    eval_metric = "logloss"         # Evaluation metric
  )
  
  # Step 8: Train the XGBoost model
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,  # Number of boosting rounds
    missing = NA        # Handle missing values
  )
  
  # Step 9: Get feature importance
  importance_matrix <- xgb.importance(
    feature_names = colnames(as.matrix(sv_matrix_final_tb_num[, 2:ncol(sv_matrix_final_tb_num)])),
    model = xgb_model
  )
  
  # Return the importance matrix
  return(importance_matrix)
}

# Function to check if an SV overlaps with a gene
find_gene_overlap <- function(SV_start, SV_end, genes_df) {
  overlapping_genes <- genes_df %>%
    filter((V4 <= SV_end) & (V5 >= SV_start)) %>%
    pull(Name) # Get the overlapping gene names
  if (length(overlapping_genes) == 0) {
    return("") # No overlap
  } else {
    return(paste(overlapping_genes, collapse = ",
")) # Combine multiple overlapping genes into one string
  }
}

plot_features <- function(bdq_feat,col,genes=genes_names,adj=1.1,wid=0.6,tit="BDQ",maxplot=(max(bdq_feat$Gain*100)+0.1)) {
bdq_feat <- bdq_feat %>%
  mutate(
    total_r = (freq_r / freq_r_tot),
    total_s = (freq_s / freq_s_tot),
    total = total_r + total_s,
    proportion_r = total_r / total,
    proportion_s = total_s / total
  ) %>%
  mutate(
    SV_start = as.numeric(str_extract(Feature, "(?<=X)\\d+")),
    SV_length = as.numeric(str_extract(Feature, "(?<=\\.)\\d+$")),
    SV_type = str_extract(Feature, "(INS|DEL|DUP)"),
    SV_end = case_when(
      SV_type == "DEL" ~ SV_start + SV_length,   # DEL means SV_end = SV_start + length
      TRUE ~ SV_start                             # INS and DUP have the same start and end
    )
  ) %>%
  rowwise() %>%
  mutate(Gene = find_gene_overlap(SV_start, SV_end, genes)) %>%
  ungroup()
if (tit == "INH") {
  bdq_feat[bdq_feat$Feature=="X2155862.DEL.7945","Gene"] = "...,
katG,
..."
}
if (tit == "PZA") {
  bdq_feat[bdq_feat$Feature=="X2275483.DEL.13977","Gene"] = "...,
pncA,
..."
  bdq_feat[bdq_feat$Feature=="X2287446.DEL.2968","Gene"] = "...,
pncA,
..."
}
# Prepare data for ggplot
plot_data <- bdq_feat %>%
  dplyr::select(Gain,Feature, proportion_r, proportion_s,Gene) %>%
  pivot_longer(cols = starts_with("proportion"), names_to = "Category", values_to = "Proportion")

bdq_loli = ggplot(bdq_feat, aes(x = Gain*100, y = reorder(str_replace(Feature,"X",""), Gain))) +
  geom_segment(aes(xend = 0, yend = str_replace(Feature,"X","")), color = col) +
  geom_point(size = 5, color = col) +
  labs(x = "SV Importance (%)",title = tit) +
  theme_pubr() +
  labs(y="") +
  theme(plot.title = element_text(size=45,hjust = 0.5),axis.text = element_text(size = 20), axis.ticks = element_line(linewidth = 3),
        legend.text = element_text(size = 10),axis.title.x = element_text(size=20)) +
  scale_x_continuous(limits = c(0, maxplot), expand = c(0, 0))

bdq_bar = ggplot(plot_data, aes(x = reorder(str_replace(Feature,"X",""), Gain), y = Proportion, fill = Category)) +
  geom_bar(stat = "identity",color="black",width = wid) +
  scale_fill_manual(values = c("proportion_r" = "#D43925", "proportion_s" = "#2376B2"),
                    labels = c("proportion_r" = "Resistant", "proportion_s" = "Susceptible")) +
  theme_pubr() +
  coord_flip() +
  labs(x="",y="Proportion") +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
        legend.title = element_blank(),axis.text = element_text(size = 10), axis.ticks = element_line(linewidth = 3),
        legend.text = element_text(size = 10)) +
  geom_text(aes(x = reorder(str_replace(Feature,"X",""), Gain), y = adj, label = Gene), size = 5, vjust = 0.5)  # Add gene annotations

#pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/bdq_feat.pdf",height = 5,width = 15)
return(bdq_loli)
}

get_pvalues <- function(drug,sv_matrix_pheno_lin,df) {
  drug_col_index <- which(colnames(sv_matrix_pheno_lin) == drug)
  
  # Step 2: Filter and mutate the data based on the drug column index
  sv_matrix_pheno_lin_mod <- sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin[[drug]]), c(drug_col_index, 23:ncol(sv_matrix_pheno_lin))] %>% 
    mutate_all(as.factor)
  
  # Step 3: Remove columns where factors have fewer than 2 levels
  sv_matrix_final <- sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
  
  # Step 4: Add drug resistance prediction column
  sv_matrix_final_tb <- sv_matrix_final %>% 
    mutate(DR_prediction = ifelse(rownames(sv_matrix_final) %in% tbprofiler[tbprofiler$drugs == drug, ]$sample, 1, 0))
  
  # Step 5: Convert the data to numeric and adjust the drug column
  sv_matrix_final_tb_num <- sv_matrix_final_tb %>% 
    mutate_all(as.numeric) %>% 
    mutate(!!drug := ifelse(!!sym(drug) == 1, 1, 0))
  
  for (i in 2:length(colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2]))) {
    #if (i==985){next}
    formula = as.formula(paste(drug," ~ ", colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i],"+ X + V2 + V3 + V4 + V5 + DR_prediction"))
    if (length(table(sv_matrix_final_tb_num[sv_matrix_final_tb_num[,colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i]]==2,drug])) == 1) {
      #if S, not R or only 1 R
      if ((table(sv_matrix_final_tb_num[sv_matrix_final_tb_num[,colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i]]==2,drug]))[1]==1) {next}
      if (names((table(sv_matrix_final_tb_num[sv_matrix_final_tb_num[,colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i]]==2,drug]))[1]) == 0) {next}
      tryCatch({
        results = logistf(formula, data = sv_matrix_final_tb_num, family = binomial)
      }, error = function(e) {
        NULL
      })
    }
    else if ((table(sv_matrix_final_tb_num[sv_matrix_final_tb_num[,colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i]]==2,drug]))[2]==1) {next}
    else {
      tryCatch({
        results = logistf(formula, data = sv_matrix_final_tb_num, family = binomial)
      }, error = function(e) {
        NULL
        })
    }
    #invisible(capture.output({
    new_row = data.frame(Drug=drug,SV=colnames(sv_matrix_final_tb_num[,2:ncol(sv_matrix_final_tb_num)-2])[i],Pval=summary(results)$prob[2])
    #}))
    colnames(df) = c("Drug","SV","Pval")
    colnames(new_row) = c("Drug","SV","Pval")
    df = rbind(df,new_row)
  }
  df$Total_SVs = ncol(sv_matrix_final_tb_num)-3
  df = df %>% mutate(Adj_pval=p.adjust(Pval,method="BH")) %>% filter(Adj_pval<0.05)
  return(df)
}
#This function is for filtering out SVs according to the following conditions:
#  - SV is at least double the frequency in the reistant population than in the susceptible population
#  - SV is not in PE/PPE genes
#  - SV is not small and in a complex region
#  - SV is not lineage-specific (we print the lineages that have the SVs to double check)
filter_sv <- function(dr_sv_pval,drug,marin,peppe,sv_matrix_pheno_lin) {
  for (i in dr_sv_pval$SV) {
  y=F
  x=F
  if (i == "DR_prediction" || i == "V2" || i == "V3" || i == "V4"|| i == "V5"|| i == drug) {
    next
  }
  #if SV in low mappable region and and a small DEL
  for (n in 1:nrow(marin)) {
    if (as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) > marin[n,]$start & as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) < marin[n,]$end & strsplit(i,".",fixed=T)[[1]][2] == "DEL" & as.numeric(strsplit(i,".",fixed=T)[[1]][3]) < 150) {
      y=T
    }
  }
  #if SV in PE/PPE region
  if (y==T) {next}
  for (n in 1:nrow(peppe)) {
    if (as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) > peppe[n,]$start & as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) < peppe[n,]$end) {x=T}
  }
  if (x==T) {next}
  #if SV is only in R
  if (length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))==1) {
    print(i)
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
    next
  }
  #if SV's frequency in R is not at least double that of the frequency in S
  resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", i]
  resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", i]
  susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
  if (resistant_freq/2 < susceptible_freq) {next}
  print(i)
  print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug])/table(sv_matrix_pheno_lin[[drug]]))
  print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
  print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
}
}
#Get only SVs in PE/PPE genes
filter_sv_peppe <- function(dr_sv_pval,drug,marin,peppe,sv_matrix_pheno_lin) {
  for (i in dr_sv_pval$SV) {
    y=F
    x=F
    if (i == "DR_prediction" || i == "V2" || i == "V3"|| i == "V4"|| i == "V5"|| i == drug) {
      next
    }
    #if SV in low mappable region and and a small DEL
    for (n in 1:nrow(marin)) {
      if (as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) > marin[n,]$start & as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) < marin[n,]$end & strsplit(i,".",fixed=T)[[1]][2] == "DEL" & as.numeric(strsplit(i,".",fixed=T)[[1]][3]) < 150 & dr_sv_pval[dr_sv_pval$SV==i,]$Pval>1e-15) {
        y=T
      }
    }
    #if SV in PE/PPE region
    if (y==T) {next}
    for (n in 1:nrow(peppe)) {
      if (as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) > peppe[n,]$start & as.numeric(str_replace(strsplit(i,".",fixed=T)[[1]][1],"X","")) < peppe[n,]$end) {x=T}
    }
    if (x==F) {next}
    #if SV is only in R
    if (length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))==1) {
      print(i)
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
      next
    }
    #if SV's frequency in R is not at least double that of the frequency in S
    resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", i]
    resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
    # For Susceptible isolates
    susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", i]
    susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
    if (resistant_freq/2 < susceptible_freq) {next}
    print(i)
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug])/table(sv_matrix_pheno_lin[[drug]]))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
  }
}
#Functions for plotting SV table
r_color_scale <- function(x) {
  scaled_value <- x / 1  # Assuming max value in the column is 1 (no full fill if less than 1)
  color_bar("#D43925", fun = function(x) scaled_value)(x)
}
s_color_scale <- function(x) {
  ifelse(x == 0, 
         sprintf("%.8f", x),  # Display the number if x is 0
         color_bar("#2376B2", fun = function(x) x / 1)(x))
}
scale_function <- function(x) (x - min(x)) / (max(x) - min(x))
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  #w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(f, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

filter_gene <- function(dr_sv_pval,drug,marin,peppe,sv_matrix_pheno_lin) {
  for (i in dr_sv_pval$SV) {
    y=F
    x=F
    if (i == "DR_prediction" || i == "V2" || i == "V3"|| i == "V4"|| i == "V5"|| i == drug) {
      next
    }
    #if gene is PE/PPE family
    if (grepl("PE|PPE|PGRS",i)) {next}
    #if gene is only in R
    if (length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))==1) {
      print(i)
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug]))
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
      print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
      next
    }
    #if gene's frequency in R is not at least double that of the frequency in S
    resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", i]
    resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
    # For Susceptible isolates
    susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", i]
    susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
    if (resistant_freq/2 < susceptible_freq) {next}
    print(i)
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1,drug])/table(sv_matrix_pheno_lin[[drug]]))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="R",]$X))
    print(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[,i]==1 & sv_matrix_pheno_lin[[drug]]=="S",]$X))
  }
}

perm <- function(drug,sv_matrix_pheno_lin,sv,pval) {
  drug_col_index <- which(colnames(sv_matrix_pheno_lin) == drug)
  sv_col_index <- which(colnames(sv_matrix_pheno_lin) == sv)
  
  # Step 2: Filter and mutate the data based on the drug column index
  sv_matrix_pheno_lin_mod <- sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin[[drug]]) & !is.na(sv_matrix_pheno_lin[[sv]]), c(drug_col_index,sv_col_index, 2493:ncol(sv_matrix_pheno_lin))] %>% 
    mutate_all(as.factor)
  
  # Step 3: Remove columns where factors have fewer than 2 levels
  sv_matrix_final <- sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
  
  # Step 4: Add drug resistance prediction column
  sv_matrix_final_tb <- sv_matrix_final %>% 
    mutate(DR_prediction = ifelse(rownames(sv_matrix_final) %in% tbprofiler[tbprofiler$drugs == drug, ]$sample, 1, 0))
  
  # Step 5: Convert the data to numeric and adjust the drug column
  sv_matrix_final_tb_num <- sv_matrix_final_tb %>% 
    mutate_all(as.numeric) %>% 
    mutate(!!drug := ifelse(!!sym(drug) == 1, 1, 0))
  
  perm_results = numeric()
  
  for (i in 1:200) {
    sv_matrix_final_tb_num[[drug]] = sample(sv_matrix_final_tb_num[[drug]])
    formula = as.formula(paste(drug," ~ ",sv,"+ X + V2 + V3 + V4 + V5 + DR_prediction"))
    if (length(table(sv_matrix_final_tb_num[sv_matrix_final_tb_num[,sv]==2,drug])) == 1) {
      results = logistf(formula, data = sv_matrix_final_tb_num, family = binomial)
    }
    else {
      results = logistf(formula, data = sv_matrix_final_tb_num, family = binomial)
    }
    perm_results[i] = summary(results)$prob[2]
  }
  #perm_adj = p.adjust(perm_results,method="BH")
  print(drug)
  print(sv)
  return((table(perm_results < pval))[[1]])
}





#We read the DR phenotype information
pheno = read.csv("~/Documents/PHD/SV/SV_analysis/DR/illumina.samplesheet.pass.csv",header = T,na.strings = "")
#List of isolates with a VCF file
final_list = read.table("/home/student.unimelb.edu.au/acanaldabalt/Documents/PHD/SV/SV_analysis/DR/vcf_list.txt",header = F)
#Unmappable regions of the genome
marin = read.table("~/Documents/PHD/SV/SV_analysis/DR/RLC_Regions.Plus.LowPmapK50E4.H37Rv.txt",sep="\t",col.names = c("start","end","length"))
marin = marin[,c(1,2)]
#List of PE/PPE genes and their positions
peppe = read.table("~/Documents/PHD/SV/SV_analysis/temp/peppe.txt",sep=" ",col.names = c("start","end"))
#We only keep the phenotypes of those isolates where we've managed to obtain a vcf
pheno = pheno[pheno$run %in% final_list$V1,]
#dlm = read.table("~/Documents/PHD/SV/SV_analysis/DR/dlm_genomes2.txt")
#bdq = read.table("~/Documents/PHD/SV/SV_analysis/DR/bdq_genomes2.txt")
#sv_matrix_pheno_lin_gene_final$delamanid <- ifelse(rownames(sv_matrix_pheno_lin_gene_final) %in% dlm$V1, "R", sv_matrix_pheno_lin_gene_final$delamanid)
#sv_matrix_pheno_lin_gene_final$bedaquiline <- ifelse(rownames(sv_matrix_pheno_lin_gene_final) %in% bdq$V1, "R", sv_matrix_pheno_lin_gene_final$bedaquiline)
#sv_matrix_pheno_lin$delamanid <- ifelse(rownames(sv_matrix_pheno_lin) %in% dlm$V1, "R", sv_matrix_pheno_lin$delamanid)
#sv_matrix_pheno_lin$bedaquiline <- ifelse(rownames(sv_matrix_pheno_lin) %in% bdq$V1, "R", sv_matrix_pheno_lin$bedaquiline)
#sv_matrix_pheno_lin_gene_final <- sv_matrix_pheno_lin_gene_final[!grepl("^17", rownames(sv_matrix_pheno_lin_gene_final)), ]
#sv_matrix_pheno_lin <- sv_matrix_pheno_lin[!grepl("^17", rownames(sv_matrix_pheno_lin)), ]
##
#Fig a: Distribution of DR across our dataset
##
pheno_long = pheno[,-c(4,5)] %>% pivot_longer(cols = -run,  # All columns except 'Isolate'
               names_to = "Drug",  # New column to hold the drug names
               values_to = "Phenotype")
df_summary <- as.data.frame(table(pheno_long$Drug, pheno_long$Phenotype, useNA = "ifany"))
colnames(df_summary) <- c("Drug", "Phenotype", "Count")

drug_acronyms <- c(
  "amikacin" = "AMK", "bedaquiline" = "BDQ", "capreomycin" = "CAP", "ciprofloxacin" = "CIP", 
  "clofazimine" = "CFZ", "cycloserine" = "CS", "delamanid" = "DLM", "ethambutol" = "EMB",
  "ethionamide" = "ETA", "gatifloxacin" = "GAT", "isoniazid" = "INH", "kanamycin" = "KAN", 
  "levofloxacin" = "LEV", "linezolid" = "LZD", "moxifloxacin" = "MFX", "ofloxacin" = "OFX", 
  "para.aminosalicylic_acid" = "PAS", "pyrazinamide" = "PZA", "rifabutin" = "RBT",  
  "rifampicin" = "RIF", "streptomycin" = "STR", "thioacetazone" = "TAC"
)

pie_dr = ggplot(df_summary %>% mutate(Phenotype = recode(Phenotype, "R" = "Resistant", "S" = "Susceptible")) %>% mutate(Drug = recode(Drug, !!!drug_acronyms)), 
       aes(x = "", y = Count, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # Black borders for slices
  coord_polar(theta = "y") +
  facet_wrap(~ Drug, ncol = 11) +  # Arrange in a grid, adjust ncol as needed
  theme_minimal(base_size = 15) +  # Minimalist theme with larger base font size
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 40, face = "bold"),
    strip.text = element_text(size = 30),  # Facet labels
    legend.position = "bottom",  # Legend at the bottom
    legend.title = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 30)
  ) +
  scale_fill_manual(values = c("#D43925","#2376B2","#d0d6d6"))
  
## Now we load the SV information
files_vcf = list.files(path='~/Documents/PHD/SV/SV_analysis/DR/vcf/',pattern = ".*vcf")
file.pathsv <<- as.vector(paste("~/Documents/PHD/SV/SV_analysis/DR/vcf", files_vcf, sep = "/"))
tbprofiler = read.table("~/Documents/PHD/SV/SV_analysis/DR/tbprofiler.txt",sep="\t",header = T)
tbprofiler = tbprofiler[,c(1,11:ncol(tbprofiler))]
tbprofiler$drugs = apply(tbprofiler[, 2:19], 1, function(row) {
  # Get column names for entries that are not "-"
  resistant_drugs <- names(tbprofiler)[which(row != "-") + 1]  # +1 adjusts for the sample column
  paste(resistant_drugs, collapse = ";")
})
tbprofiler = tbprofiler %>% separate_rows(drugs,sep = ";")
tbprofiler = tbprofiler[,c(1,20)]
tbp_cip = read.table("~/Documents/PHD/SV/SV_analysis/DR/tbp_cip.tsv",header = F,sep="\t",col.names = c("sample","drugs"))
tbprofiler = rbind(tbprofiler,tbp_cip)
tabledr <- lapply(X = file.pathsv, FUN = read.table, header = FALSE,sep = "\t",comment.char='#')

# Collect all SVs from all VCF files
all_svs <- list()
for (i in 1:length(tabledr)) {
  sv_data <- process_mtb_vcf(tabledr[[i]])
  all_svs <- rbind(all_svs, sv_data[sv_data$ID != "NaN", ])
}

all_svs2 = all_svs[!apply(all_svs == "", 1, all),]
#all_svs2 = all_svs[rowSums(is.na(all_svs)) < 2, ]

# Merge SVs across all VCF files
merged_svs <- merge_svs(all_svs2)

# Create the matrix with samples as rows and merged SVs as columns
sample_names <- str_replace(files_vcf,"_ref_short.vcf","")
sv_matrix <- matrix(0, nrow = length(files_vcf), ncol = nrow(merged_svs),
                    dimnames = list(sample_names, paste(merged_svs$ID, merged_svs$POS, merged_svs$END, sep=".")))

# Fill the matrix
for (i in 1:length(tabledr)) {
  print(i)
  sv_data <- process_mtb_vcf(tabledr[[i]])
  sample_id <- sample_names[i]
  
  sv_data$POS <- as.numeric(as.character(sv_data$POS))
  sv_data$END <- as.numeric(as.character(sv_data$END))
  sv_data$SIZE <- as.numeric(as.character(sv_data$SIZE))
  
  for (j in seq_along(merged_svs$ID)) {
    merged_sv <- merged_svs[j, ]
    
    if (any(sv_data$ID == "NaN" & sv_data$POS <= merged_sv$END & sv_data$END >= merged_sv$POS)) {
      sv_matrix[sample_id, j] <- NA
    } else {
      overlap_result <- mapply(check_overlap, list(merged_sv), split(sv_data, seq(nrow(sv_data))))
      sv_matrix[sample_id, j] <- ifelse(any(overlap_result), 1, 0)
    }
  }
}

# After running in HPC, read matrix
sv_matrix = read.table("~/Documents/PHD/SV/SV_analysis/DR/sv_matrix.txt")
sv_matrix_pheno = merge(pheno[,-c(4,5)],sv_matrix,by.y="row.names",by.x=1) %>% column_to_rownames(var="run")
sv_lin = read.table("~/Documents/PHD/SV/SV_analysis/DR/summary_lin.csv",sep=",",header = T)
sv_matrix_pheno_lin = merge(sv_matrix_pheno,sv_lin,by.x="row.names",by.y=1) %>% column_to_rownames(var="Row.names")
colnames(sv_lin) = c("ID","Lineage")
sv_lin$Lineage = sapply(str_split(sv_lin$Lineage, "[.]"), `[`, 1)
sv_lin$Lineage = sapply(str_split(sv_lin$Lineage, "[;]"), `[`, 1)

set.seed(123)

#We create the PCA to control for population structure
snps_pca = read.table("~/Documents/PHD/SV/SV_analysis/DR/loadings_hail2.tsv",sep="\t")
lin_color2 = c("#DDCC77","#117733","#88CCEE","#CC6677","#882255","#AA4499","#332288","#DDAA55")
pca_all = ggplot(merge(snps_pca,sv_lin[!grepl(";",sv_lin$Lineage),],by.x=1,by.y=1),aes(x=V2,y=V3,color=Lineage)) +
  geom_point() +
  theme_pubr() +
  scale_color_manual(values=lin_color) +
  ylab("PC2(13.6%)")+
  xlab("PC1(59.8%)") +
  guides(shape=FALSE,size=FALSE,color=guide_legend(title = "Lineage")) +
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 30),
        legend.text = element_text(size = 30),legend.title = element_text(size = 30),axis.ticks = element_line(linewidth = 3)) +
  guides(color = guide_legend(override.aes = list(size = 8)))

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/42k_pca.pdf",height = 10,width = 10)
pca_all
dev.off()

snps_pca = snps_pca[,1:5]
sv_matrix_pheno_lin = merge(sv_matrix_pheno_lin,snps_pca,by.x="row.names",by.y=1) %>% column_to_rownames(var="Row.names")

#Logistic regression for pval + XGBOOST 
genes = read.table("~/Downloads/Mycobacterium_tuberculosis_H37Rv_gff_v4.1.gff",sep=" ")
genes$V9 = as.character(genes$V8)
# Apply the function to the data frame and convert to a wide format
genes_names <- genes %>%
  dplyr::mutate(parsed_info = purrr::map(V9, parse_info)) %>%
  unnest_wider(parsed_info) %>%
  filter(Name != "ncRv3520")

#We now find SVs associated to DR using logistic regression
dr_sv_pval_amk = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_amk = get_pvalues(drug = "amikacin",sv_matrix_pheno_lin,dr_sv_pval_amk)
filter_sv(dr_sv_pval_amk,"amikacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_bdq = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_bdq = get_pvalues(drug = "bedaquiline",sv_matrix_pheno_lin,dr_sv_pval_bdq)
filter_sv(dr_sv_pval_bdq,"bedaquiline",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_cap = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_cap = get_pvalues(drug = "capreomycin",sv_matrix_pheno_lin,dr_sv_pval_cap)
filter_sv(dr_sv_pval_cap,"capreomycin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_cip = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_cip = get_pvalues(drug = "ciprofloxacin",sv_matrix_pheno_lin,dr_sv_pval_cip)
filter_sv(dr_sv_pval_cip,"ciprofloxacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_cfz = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_cfz = get_pvalues(drug = "clofazimine",sv_matrix_pheno_lin,dr_sv_pval_cfz)
filter_sv(dr_sv_pval_cfz,"clofazimine",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_cs = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_cs = get_pvalues(drug = "cycloserine",sv_matrix_pheno_lin,dr_sv_pval_cs)
filter_sv(dr_sv_pval_cs,"cycloserine",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_dlm = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_dlm = get_pvalues(drug = "delamanid",sv_matrix_pheno_lin,dr_sv_pval_dlm)
filter_sv(dr_sv_pval_dlm,"delamanid",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_emb = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_emb = get_pvalues(drug = "ethambutol",sv_matrix_pheno_lin,dr_sv_pval_emb)
filter_sv(dr_sv_pval_emb,"ethambutol",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_eta = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_eta = get_pvalues(drug = "ethionamide",sv_matrix_pheno_lin,dr_sv_pval_eta)
filter_sv(dr_sv_pval_eta,"ethionamide",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_gat = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_gat = get_pvalues(drug = "gatifloxacin",sv_matrix_pheno_lin,dr_sv_pval_gat)
filter_sv(dr_sv_pval_gat,"gatifloxacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_inh = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_inh = get_pvalues(drug = "isoniazid",sv_matrix_pheno_lin,dr_sv_pval_inh)
filter_sv(dr_sv_pval_inh,"isoniazid",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_kan = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_kan = get_pvalues(drug = "kanamycin",sv_matrix_pheno_lin,dr_sv_pval_kan)
filter_sv(dr_sv_pval_kan,"kanamycin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_lev = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_lev = get_pvalues(drug = "levofloxacin",sv_matrix_pheno_lin,dr_sv_pval_lev)
filter_sv(dr_sv_pval_lev,"levofloxacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_lzd = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_lzd = get_pvalues(drug = "linezolid",sv_matrix_pheno_lin,dr_sv_pval_lzd)
filter_sv(dr_sv_pval_lzd,"linezolid",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_pas = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_pas = get_pvalues(drug = "para.aminosalicylic_acid",sv_matrix_pheno_lin,dr_sv_pval_pas)
filter_sv(dr_sv_pval_pas,"para.aminosalicylic_acid",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_mfx = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_mfx = get_pvalues(drug = "moxifloxacin",sv_matrix_pheno_lin,dr_sv_pval_mfx)
filter_sv(dr_sv_pval_mfx,"moxifloxacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_ofx = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_ofx = get_pvalues(drug = "ofloxacin",sv_matrix_pheno_lin,dr_sv_pval_ofx)
filter_sv(dr_sv_pval_ofx,"ofloxacin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_pza = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_pza = get_pvalues(drug = "pyrazinamide",sv_matrix_pheno_lin,dr_sv_pval_pza)
filter_sv(dr_sv_pval_pza,"pyrazinamide",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_rbt = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_rbt = get_pvalues(drug = "rifabutin",sv_matrix_pheno_lin,dr_sv_pval_rbt)
filter_sv(dr_sv_pval_rbt,"rifabutin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_rif = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_rif = get_pvalues(drug = "rifampicin",sv_matrix_pheno_lin,dr_sv_pval_rif)
filter_sv(dr_sv_pval_rif,"rifampicin",marin,peppe,sv_matrix_pheno_lin)

dr_sv_pval_str = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_sv_pval_str = get_pvalues(drug = "streptomycin",sv_matrix_pheno_lin,dr_sv_pval_str)
filter_sv(dr_sv_pval_str,"streptomycin",marin,peppe,sv_matrix_pheno_lin)

#dr_sv_pval$Pval_adj = p.adjust(dr_sv_pval$Pval, method = "fdr")

#Now we also find SVs using xgboost and keep SVs that are found by both xgboost and logistic regression methods to keep only confident SVs - we select those SVs that passed the logistic regression filters
#This will help prioritize SVs differently and will give us the importance of TB-Profiler's predictions in each model
#AMK
amk_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="amikacin")
#DR_pred accounts for 69.8% importance, all SVs are < 1%
amk_feat = amk_feat[c(137),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$amikacin)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$amikacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$amikacin)) >= 1,
                      table(sv_matrix_pheno_lin$amikacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$amikacin)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$amikacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$amikacin)) >= 2,
                      table(sv_matrix_pheno_lin$amikacin)[[2]], 0)
)
amk_feat$Drug = "AMK"
amk_plot = plot_features(amk_feat,col="#e7298a",tit="AMK",wid=0.1,maxplot = 0.16)

#BDQ
bdq_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="bedaquiline")
#DR_pred accounts for 20.6% importance


#CAP
cap_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="capreomycin")
#DR_pred accounts for 57% importance all SVs are < 1%

#CIP
cip_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="ciprofloxacin")
#DR_pred accounts for 79% importance. No SVs that passed the filters

#CFZ
cfz_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="clofazimine")
#DR_pred accounts for 2.7% importance
cfz_feat = cfz_feat[c(83,149),]
new_row <- data.frame(Feature = "X2208005.DEL.12719", Gain = 0, Cover = 0, Frequency = 0)
cfz_feat = rbind(cfz_feat,new_row)
cfz_feat = cfz_feat %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$clofazimine)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$clofazimine)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$clofazimine)) >= 1,
                      table(sv_matrix_pheno_lin$clofazimine)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$clofazimine)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$clofazimine)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$clofazimine)) >= 2,
                      table(sv_matrix_pheno_lin$clofazimine)[[2]], 0)
)
cfz_feat$Drug = "CFZ"
cfz_plot = plot_features(cfz_feat,col="#9467bd",wid=0.2,tit="CFZ",maxplot = 0.42)

#CS
cs_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="cycloserine")
#DR_pred accounts for 0% importance, SVs did not pass filters

#DLM
#DR_pred accounts for 3% importance
dlm_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="delamanid")
dlm_feat = dlm_feat[c(136,200,214),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$delamanid)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$delamanid)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$delamanid)) >= 1,
                      table(sv_matrix_pheno_lin$delamanid)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$delamanid)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$delamanid)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$delamanid)) >= 2,
                      table(sv_matrix_pheno_lin$delamanid)[[2]], 0)
)
dlm_feat$Drug = "DLM"
dlm_plot = plot_features(dlm_feat,col="#bcbd22",wid=0.35,tit="DLM",maxplot = 0.19)

#EMB
#DR_pred accounts for 70% importance
emb_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="ethambutol")
emb_feat = emb_feat[c(41,283),]
new_row <- data.frame(Feature = "X4137120.DEL.1356", Gain = 0, Cover = 0, Frequency = 0)
emb_feat = rbind(emb_feat,new_row)
new_row <- data.frame(Feature = "X475155.DEL.1599", Gain = 0, Cover = 0, Frequency = 0)
emb_feat = rbind(emb_feat,new_row)
emb_feat = emb_feat %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethambutol)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethambutol)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$ethambutol)) >= 1,
                      table(sv_matrix_pheno_lin$ethambutol)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethambutol)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethambutol)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$ethambutol)) >= 2,
                      table(sv_matrix_pheno_lin$ethambutol)[[2]], 0)
)
emb_feat$Drug = "EMB"
emb_plot = plot_features(emb_feat,col="#17becf",wid=0.4,tit="EMB",maxplot = 0.22)

#ETA
#DR_pred accounts for 50.2% importance
eta_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="ethionamide")
eta_feat = eta_feat[c(17,477),]
new_row <- data.frame(Feature = "X1182320.DEL.56", Gain = 0, Cover = 0, Frequency = 0)
eta_feat = rbind(eta_feat,new_row)
eta_feat = eta_feat %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethionamide)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethionamide)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$ethionamide)) >= 1,
                      table(sv_matrix_pheno_lin$ethionamide)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethionamide)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$ethionamide)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$ethionamide)) >= 2,
                      table(sv_matrix_pheno_lin$ethionamide)[[2]], 0)
)
eta_feat$Drug = "ETA"
eta_plot = plot_features(eta_feat,col="#3288bd",wid=0.4,tit="ETA",maxplot = 0.42)

#GAT
#DR_pred accounts for 0% importance
gat_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="gatifloxacin")

#INH
#DR pred accounts for 85.8%
inh_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="isoniazid")
inh_feat = inh_feat[c(11,43,62,159,372),]
new_row <- data.frame(Feature = "X3069892.DEL.102", Gain = 0, Cover = 0, Frequency = 0)
inh_feat = rbind(inh_feat,new_row)
new_row <- data.frame(Feature = "X2245629.DUP.503", Gain = 0, Cover = 0, Frequency = 0)
inh_feat = rbind(inh_feat,new_row)
inh_feat = inh_feat %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$isoniazid)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$isoniazid)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$isoniazid)) >= 1,
                      table(sv_matrix_pheno_lin$isoniazid)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$isoniazid)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$isoniazid)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$isoniazid)) >= 2,
                      table(sv_matrix_pheno_lin$isoniazid)[[2]], 0)
)
inh_feat$Drug = "INH"
inh_plot = plot_features(inh_feat,col="#66c2a5",wid=0.4,tit="INH",maxplot = 0.12)

#KAN
kan_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="kanamycin")
#DR pred accounts for 68.1%
kan_feat = kan_feat[c(63,105),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$kanamycin)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$kanamycin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$kanamycin)) >= 1,
                      table(sv_matrix_pheno_lin$kanamycin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$kanamycin)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$kanamycin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$kanamycin)) >= 2,
                      table(sv_matrix_pheno_lin$kanamycin)[[2]], 0)
)
kan_feat$Drug = "KAN"
kan_plot = plot_features(kan_feat,col="#d95f02",wid=0.25,tit="KAN",maxplot = 0.17)

#LEV
lev_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="levofloxacin")
#DR pred accounts for 76.9%
lev_feat = lev_feat[c(334),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$levofloxacin)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$levofloxacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$levofloxacin)) >= 1,
                      table(sv_matrix_pheno_lin$levofloxacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$levofloxacin)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$levofloxacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$levofloxacin)) >= 2,
                      table(sv_matrix_pheno_lin$levofloxacin)[[2]], 0)
)
lev_feat$Drug = "LEV"
lev_plot = plot_features(lev_feat,col="#7570b3",wid=0.25,tit="LEV",maxplot = 0.07)

#LZD
lzd_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="linezolid")
#DR pred accounts for 22.3%, no SVs passed the filters

#MFX
mfx_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="moxifloxacin")
#DR pred accounts for 62.4%
new_row <- data.frame(Feature = "X1443721.DEL.1351", Gain = 0, Cover = 0, Frequency = 0)
mfx_feat = rbind(mfx_feat,new_row)
mfx_feat = mfx_feat[c(543),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$moxifloxacin)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$moxifloxacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$moxifloxacin)) >= 1,
                      table(sv_matrix_pheno_lin$moxifloxacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$moxifloxacin)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$moxifloxacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$moxifloxacin)) >= 2,
                      table(sv_matrix_pheno_lin$moxifloxacin)[[2]], 0)
)
mfx_feat$Drug = "MFX"
mfx_plot = plot_features(mfx_feat,col="#e8a",wid=0.35,tit="MFX",maxplot = 0.018)

#OFX
ofx_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="ofloxacin")
#DR pred accounts for 64%

#PAS
pas_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="para.aminosalicylic_acid")
#DR pred accounts for 14.4%

#PZA
pza_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="pyrazinamide")
#DR pred accounts for 69.2 %
pza_feat = pza_feat[c(7,21,23,45,95,145,196,358,411),]
new_row <- data.frame(Feature = "X2757647.DEL.55", Gain = 0, Cover = 0, Frequency = 0)
pza_feat = rbind(pza_feat,new_row)
new_row <- data.frame(Feature = "X378875.DEL.282", Gain = 0, Cover = 0, Frequency = 0)
pza_feat = rbind(pza_feat,new_row)
pza_feat = pza_feat %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$pyrazinamide)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$pyrazinamide)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$pyrazinamide)) >= 1,
                      table(sv_matrix_pheno_lin$pyrazinamide)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$pyrazinamide)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$pyrazinamide)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$pyrazinamide)) >= 2,
                      table(sv_matrix_pheno_lin$pyrazinamide)[[2]], 0)
)
pza_feat$Drug = "PZA"
pza_plot = plot_features(pza_feat,col="#b3b3b3",wid=0.35,tit="PZA",maxplot = 0.42)

#RBT
rbt_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="rifabutin")
#DR pred accounts for 0%

#RIF
rif_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="rifampicin")
#DR pred accounts for 88.9 %

#STR
str_feat = get_importance_matrix(sv_matrix_pheno_lin,tbprofiler,drug="streptomycin")
#DR pred accounts for 69.7 %
new_row <- data.frame(Feature = "X1284426.DUP.84", Gain = 0, Cover = 0, Frequency = 0)
str_feat = rbind(str_feat,new_row)
new_row <- data.frame(Feature = "X4407807.DEL.151", Gain = 0, Cover = 0, Frequency = 0)
str_feat = rbind(str_feat,new_row)
str_feat = str_feat[c(589,590),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$streptomycin)) >= 1,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$streptomycin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin$streptomycin)) >= 1,
                      table(sv_matrix_pheno_lin$streptomycin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$streptomycin)) >= 2,
                  table(sv_matrix_pheno_lin[sv_matrix_pheno_lin[[Feature]] == 1,]$streptomycin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin$streptomycin)) >= 2,
                      table(sv_matrix_pheno_lin$streptomycin)[[2]], 0)
)
str_feat$Drug = "STR"
str_plot = plot_features(str_feat,col="#2ca02c",wid=0.35,tit="STR",maxplot = 0.018)

#WHO v2 Mutations Importance
who_plot = data.frame(Importance=c(69.8,20.6,57,79,2.7,0,3,70,50.2,0,85.8,68.1,76.9,22.3,62.4,64,14.4,69.2,0,88.9,69.7),
           Drug=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MFX","OFX","PAS","PZA",
                  "RBT","RIF","STR")) %>%
  ggplot(.,aes(x=Importance,y=reorder(Drug,Importance))) +
  geom_segment(aes(xend = 0, yend = Drug, color = Drug)) +
  geom_point(size = 5, aes(color = Drug)) +
  labs(x = "TB-Profiler mutation list Importance (%)") +
  theme_pubr() +
  labs(y="") +
  theme(plot.title = element_text(size=35,hjust = 0.5)) +
  scale_color_manual(values = c("AMK" = "#e7298a","BDQ" = "#ff7f0e","CAP" = "#1f77b4","CIP" = "#d62728","CFZ"="#9467bd",
                                "CS" = "#8c564b","DLM" = "#bcbd22","EMB"="#17becf","ETA" = "#3288bd","GAT" = "#ffd92f",
                                "INH" = "#66c2a5","KAN" = "#d95f02","LEV" = "#7570b3","LZD" ="#a6d854","MFX"="#e8a",
                                "OFX"="#e377c2","PAS" = "#e78ac3","PZA" = "#b3b3b3","RBT"="#7f7f7f","RIF"="#fdae61","STR"="#2ca02c")) +
  theme(legend.position = "none",plot.title = element_text(size=45,hjust = 0.5),axis.text = element_text(size = 20), axis.ticks = element_line(linewidth = 3),
        legend.text = element_text(size = 10),axis.title.x = element_text(size=20)) +
  scale_x_continuous(limits = c(0, 100))


#Final SV Table with p-values (adjusted significant + non-significant but in associated genes)
pval_tab = do.call("rbind",list(amk_feat,cfz_feat,dlm_feat,emb_feat,eta_feat,inh_feat,kan_feat,lev_feat,mfx_feat,pza_feat,str_feat))
pval_tab = pval_tab[,c(9,1,5,6,7,8)]
for (i in 1:nrow(pval_tab)) {
  # For Resistant isolates
  drug = pval_tab[i,]$Drug
  drug_ac = data.frame(short=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MFX",
                               "OFX","PAS","PZA","RBT","RIF","STR"),
                       long=c("amikacin","bedaquiline","capreomycin","ciprofloxacin","clofazimine","cycloserine","delamanid",
                              "ethambutol","ethionamide","gatifloxacin","isoniazid","kanamycin","levofloxacin","linezolid","moxifloxacin",
                              "ofloxacin","para.aminosalicylic_acid","pyrazinamide","rifabutin","rifampicin","streptomycin"),stringsAsFactors = F)
  drug = drug_ac[drug_ac$short==drug,]$long
  sv_col = pval_tab[i,]$Feature
  resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", sv_col]
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", sv_col]
  pval_tab[i,]$freq_r = sum(resistant_isolates == 1, na.rm = TRUE)
  pval_tab[i,]$freq_s = sum(susceptible_isolates == 1, na.rm = TRUE)
  pval_tab[i,]$freq_r_tot = sum(!is.na(resistant_isolates))
  pval_tab[i,]$freq_s_tot = sum(!is.na(susceptible_isolates))
}
dr_sv_all = do.call("rbind",list(dr_sv_pval_amk%>% mutate(Drug="AMK"),dr_sv_pval_cfz%>% mutate(Drug="CFZ"),
                                 dr_sv_pval_dlm%>% mutate(Drug="DLM"),dr_sv_pval_emb%>% mutate(Drug="EMB"),dr_sv_pval_eta%>% mutate(Drug="ETA"),
                                 dr_sv_pval_inh%>% mutate(Drug="INH"),dr_sv_pval_kan%>% mutate(Drug="KAN"),dr_sv_pval_lev%>% mutate(Drug="LEV"),dr_sv_pval_mfx%>% mutate(Drug="MFX"),
                                 dr_sv_pval_pza%>% mutate(Drug="PZA"),dr_sv_pval_str%>% mutate(Drug="STR")))
pval_tab = merge(pval_tab,dr_sv_all,by.x=c(1,2),by.y=c(1,2),all.x=T)
#pval_tab_nosig = pval_tab %>% filter(Adj_pval>0.05)
#pval_tab_sig = pval_tab %>% filter(Adj_pval<0.05)
#pval_tab = pval_tab %>% mutate(Adj_pval = pmin(Pval * Total_SVs / seq_along(Pval), 1))
write.table(pval_tab,"~/Documents/PHD/SV/SV_analysis/DR/sv_dr_pval.csv",sep=",")
#write.table(pval_tab_nosig,"~/Documents/PHD/SV/SV_analysis/DR/sv_dr_pval_nosig.csv",sep=",")


#TABLE part of the plot
#First, we rbind all of the features together to plot them all in the same table
all_feat = do.call("rbind",list(amk_feat,cfz_feat,dlm_feat,emb_feat,eta_feat,inh_feat,kan_feat,lev_feat,mfx_feat,pza_feat,str_feat))
all_feat = all_feat[,c(9,1,5,6,7,8)]
for (i in 1:nrow(all_feat)) {
  # For Resistant isolates
  drug = all_feat[i,]$Drug
  drug_ac = data.frame(short=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MFX",
                                           "OFX","PAS","PZA","RBT","RIF","STR"),
                       long=c("amikacin","bedaquiline","capreomycin","ciprofloxacin","clofazimine","cycloserine","delamanid",
                       "ethambutol","ethionamide","gatifloxacin","isoniazid","kanamycin","levofloxacin","linezolid","moxifloxacin",
                       "ofloxacin","para.aminosalicylic_acid","pyrazinamide","rifabutin","rifampicin","streptomycin"),stringsAsFactors = F)
  drug = drug_ac[drug_ac$short==drug,]$long
  sv_col = all_feat[i,]$Feature
  resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", sv_col]
  resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", sv_col]
  susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
  all_feat$freq_r <- as.numeric(all_feat$freq_r)
  all_feat$freq_s <- as.numeric(all_feat$freq_s)
  all_feat$freq_r_tot <- as.numeric(all_feat$freq_r_tot)
  all_feat$freq_s_tot <- as.numeric(all_feat$freq_s_tot)
  all_feat[i,]$freq_r = resistant_freq
  all_feat[i,]$freq_s = susceptible_freq
  all_feat[i,]$freq_r_tot = sum(!is.na(resistant_isolates))
  all_feat[i,]$freq_s_tot = sum(!is.na(susceptible_isolates))
}
all_feat=all_feat %>%
  left_join(pval_tab, by=c("Drug", "Feature"))
all_feat = all_feat[order(all_feat$Adj_pval)[1:20],]
all_feat= all_feat %>% mutate(Start=as.numeric(str_extract(Feature, "(?<=X)\\d+")),
                    Length=as.numeric(str_extract(Feature, "(?<=\\.)\\d+$")),Type=str_extract(Feature, "(INS|DEL|DUP)"),
                    End = case_when(
                      Type == "DEL" ~ Start + Length,   # DEL means SV_end = SV_start + length
                      TRUE ~ Start                             # INS and DUP have the same start and end
                    ),
                    Gene = find_gene_overlap(Start, End, genes_names))
all_feat = all_feat[,c(1,14,15,16,3,5,4,6,18)]
colnames(all_feat) = c("Drug","Start","Length","Type","Resistant Frequency","Susceptible Frequency","Total Resistant","Total Susceptible","Gene/s")
all_feat$`Susceptible Frequency` <- sprintf("%.10f",all_feat$`Susceptible Frequency`)  # Show 10 decimal places
all_feat$`Susceptible Frequency` = as.numeric(all_feat$`Susceptible Frequency`)

drug_color_scale <- formatter("span",
                              style = x ~ style( display="block",
                                                 "border-radius" = "0px",
                                                 "background-color" = case_when(all_feat$Drug=="AMK" ~ "#e7298a",
                                                                all_feat$Drug=="CFZ" ~ "#9467bd",
                                                                all_feat$Drug=="DLM" ~ "#bcbd22",
                                                                all_feat$Drug=="EMB" ~ "#17becf",
                                                                all_feat$Drug=="ETA" ~ "#3288bd",
                                                                all_feat$Drug=="INH" ~ "#66c2a5",
                                                                all_feat$Drug=="KAN" ~ "#d95f02",
                                                                all_feat$Drug=="PZA" ~ "#b3b3b3",
                                                                all_feat$Drug=="LEV" ~ "#7570b3",
                                                                all_feat$Drug=="MFX" ~ "#e8a",
                                                                all_feat$Drug=="STR" ~ "#2ca02c"),
                                 "color" = "black",
                              "border-radius" = "4px"))
type_color_scale <- formatter("span",
                              style = x ~ style( display="block",
                                                 "border-radius" = "0px",
                                                 "background-color" = case_when(all_feat$Type=="DEL" ~ "lightblue",
                                                                                all_feat$Type=="DUP" ~ "#eb818a",
                                                                                all_feat$Type=="INS" ~ "#D474E0"),
                                                 "color" = "black",
                                                 "border-radius" = "4px"))

tab = formattable(all_feat, list(Length = color_bar("orange" , fun = scale_function),
                           "Resistant Frequency" = r_color_scale,
                           "Susceptible Frequency" = s_color_scale,
                           "Drug" = drug_color_scale,
                           "Type" = type_color_scale)) %>%
  as.htmlwidget(width="100%") %>%
  prependContent(tags$style("table,td,tr,th { border: 1px solid black !important;}","td, th { text-align: left !important; padding-left: 5px; }" ))
export_formattable(tab,"~/Downloads/tab_sv_20.png")

#Fig 5!
pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/Fig5_20_final_tbp.pdf",width = 25,height = 45)
plot_grid(plot_grid(pie_dr,plot_grid(amk_plot,cfz_plot,dlm_plot,emb_plot,eta_plot,inh_plot,kan_plot,lev_plot,mfx_plot,
          pza_plot,str_plot,who_plot,nrow = 3,ncol=4),nrow=2,rel_heights = c(1,3),labels = c("a","b"),label_size = 40),
          rasterGrob(png::readPNG("~/Downloads/tab_sv_20.png")),labels = c("","c"),nrow=2,ncol=1,label_size = 40,rel_heights = c(1.3,1))
dev.off()

#PE/PPE - only SVs in separate table. After running filter_sv_peppe function.
dr_sv_peppe = do.call("rbind",list(dr_sv_pval_amk%>% mutate(Drug="AMK")%>%filter(SV %in% c("X3054707.DEL.213","X2088763.DEL.348","X1217498.DEL.872","X1189598.DEL.222")),
                                 dr_sv_pval_bdq%>% mutate(Drug="BDQ")%>%filter(SV %in% c("X2796450.DEL.156")),
                                 dr_sv_pval_cfz%>% mutate(Drug="CFZ")%>%filter(SV %in% c("X2796818.DEL.147","X1982973.DEL.233")),
                                 dr_sv_pval_dlm%>% mutate(Drug="DLM")%>%filter(SV %in% c("X4032411.DEL.176")),
                                 dr_sv_pval_emb%>% mutate(Drug="EMB")%>%filter(SV %in% c("X2943895.DEL.1702","X2633894.DEL.1056","X1618350.DEL.183")),
                                 dr_sv_pval_eta%>% mutate(Drug="ETA")%>%filter(SV %in% c("X2423462.DEL.84","X836157.DEL.613")),
                                 dr_sv_pval_inh%>% mutate(Drug="INH")%>%filter(SV %in% c("X3998567.DEL.258","X3931373.DEL.389","X925788.DEL.175","X426461.DEL.216")),
                                 dr_sv_pval_kan%>% mutate(Drug="KAN")%>%filter(SV %in% c("X3949324.INS.2652","X2796450.DEL.156","X1489010.DEL.379")),
                                 dr_sv_pval_mfx%>% mutate(Drug="MFX")%>%filter(SV %in% c("X2796450.DEL.156")),
                                 dr_sv_pval_ofx%>% mutate(Drug="OFX")%>%filter(SV %in% c("X2168813.DEL.156","X1618350.DEL.183")),
                                 dr_sv_pval_pza%>% mutate(Drug="PZA")%>%filter(SV %in% c("X368038.DEL.2556","X3501665.INS.1337")),
                                 dr_sv_pval_rif%>% mutate(Drug="RIF")%>%filter(SV %in% c("X1618350.DEL.183")),
                                 dr_sv_pval_str%>% mutate(Drug="STR")%>%filter(SV %in% c("X3998567.DEL.258","X2088405.DEL.90"))))
dr_sv_peppe= dr_sv_peppe %>% mutate(Start=as.numeric(str_extract(SV, "(?<=X)\\d+")),
                              Length=as.numeric(str_extract(SV, "(?<=\\.)\\d+$")),Type=str_extract(SV, "(INS|DEL|DUP)"),
                              End = case_when(
                                Type == "DEL" ~ as.numeric(Start + Length),   # DEL means SV_end = SV_start + length
                                TRUE ~ as.numeric(Start)                             # INS and DUP have the same start and end
                              )) %>%
                              rowwise() %>%
                              mutate(Gene = find_gene_overlap(Start, End, genes_names)) %>%
                              ungroup()
write.table(dr_sv_peppe,"~/Documents/PHD/SV/SV_analysis/DR/sv_dr_pval_peppe.csv",sep=",")


#Example of finding to which drug the SV is associated in the case of MDR
which(colnames(sv_matrix_pheno_lin)=="X3663957.DUP.55")
sv_matrix_pheno_lin_mod = sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin$isoniazid) & !is.na(sv_matrix_pheno_lin$pyrazinamide), c(11,18,1984,2493,2494,2495,2496,2497)]
sv_matrix_final <- sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
summary(logistf(X3663957.DUP.55 ~ isoniazid + pyrazinamide + V2 + V3 + V4 + V5, data = sv_matrix_final, family = binomial))


#Command to find whether unexplained resistance in SVs that we've found associated to DR is due to SVs in known drug resistance-conferring genes, and not due to the new SVs
sv_matrix_pheno_lin_gene_final[rownames(sv_matrix_pheno_lin[sv_matrix_pheno_lin$X2197956.DEL.1520==1 & sv_matrix_pheno_lin$clofazimine=="R" & !is.na(sv_matrix_pheno_lin$clofazimine) & !is.na(sv_matrix_pheno_lin$X2197956.DEL.1520),])[!(rownames(sv_matrix_pheno_lin[sv_matrix_pheno_lin$X2197956.DEL.1520==1 & sv_matrix_pheno_lin$clofazimine=="R" & !is.na(sv_matrix_pheno_lin$clofazimine) & !is.na(sv_matrix_pheno_lin$X2197956.DEL.1520),]) %in% tbprofiler[tbprofiler$drugs=="clofazimine",]$sample)],"Rv0678"]

#VALIDATION OF ASSOCIATIONS THROUGH PERMUTATION TEST
pval_tab2 = pval_tab %>% mutate(Drug=case_when(Drug=="AMK" ~ "amikacin",
                                               Drug=="BDQ" ~ "bedaquiline",
                                               Drug=="CFZ" ~ "clofazimine",
                                               Drug=="DLM" ~ "delamanid",
                                               Drug=="CAP" ~ "capreomycin",
                                               Drug=="EMB" ~ "ethambutol",
                                               Drug=="ETA" ~ "ethionamide",
                                               Drug=="INH" ~ "isoniazid",
                                               Drug=="KAN" ~ "kanamycin",
                                               Drug=="MFX" ~ "moxifloxacin",
                                               Drug=="LEV" ~ "levofloxacin",
                                               Drug=="LZD" ~ "linezolid",
                                               Drug=="OFX" ~ "ofloxacin",
                                               Drug=="PZA" ~ "pyrazinamide",
                                               Drug=="STR" ~ "streptomycin",))
set.seed(123)
perm_df = data.frame(sv=as.character(),perm=as.character(),stringsAsFactors = F)
for (i in 1:nrow(pval_tab)){
  perm_df=rbind(perm_df,data.frame(sv=pval_tab2[i,1],perm=perm(pval_tab2[i,1],sv_matrix_pheno_lin,pval_tab2[i,2],pval_tab2[i,7])))
}







#
#GENE-BASED TEST - Repeat previous analysis but with genes
#

sv_matrix_pheno_lin_gene <- as.data.frame(matrix(0, nrow = nrow(sv_matrix_pheno_lin), ncol = nrow(genes_names), 
                      dimnames = list(rownames(sv_matrix_pheno_lin), c(genes_names$Name))))

extract_sv_info <- function(sv_string) {
  sv_start <- as.numeric(str_extract(sv_string, "(?<=X)\\d+"))
  sv_type <- str_extract(sv_string, "(INS|DEL|DUP)")
  sv_size <- as.numeric(str_extract(sv_string, "(?<=\\.)\\d+$"))
  sv_end <- ifelse(sv_type == "DEL",sv_start + sv_size,sv_start)
  return(list(start = sv_start, end = sv_end, type = sv_type))
}

# Loop through each SV (row) and gene, increment the gene count if it overlaps
for (sv in colnames(sv_matrix_pheno_lin[,24:ncol(sv_matrix_pheno_lin)-1])) {
  print(which(colnames(sv_matrix_pheno_lin[,24:ncol(sv_matrix_pheno_lin)-1])==sv))
  sv_info <- extract_sv_info(sv)
  sv_start <- sv_info$start
  sv_end <- sv_info$end
  # Check for overlap with each gene
  for (i in 1:nrow(genes_names)) {
    gene_start <- genes_names$V4[i]
    gene_end <- genes_names$V5[i]
    # If the SV overlaps the gene, update the gene_matrix
    if (sv_end < gene_start || sv_start > gene_end) {next}
    else {
      print(genes_names[i,]$Name)
      # Add values from sv_matrix to gene_matrix, propagate NAs
      for (isolate in 1:nrow(sv_matrix_pheno_lin)) {
        sv_value <- sv_matrix[isolate,sv]
        gene_value <- sv_matrix_pheno_lin_gene[isolate, genes_names[i,]$Name]
        # Handle NAs
        if (is.na(sv_value)) {
          sv_matrix_pheno_lin_gene[isolate, genes_names[i,]$Name] <- NA
        } else {
          sv_matrix_pheno_lin_gene[isolate, genes_names[i,]$Name] <- sv_value + gene_value
        }
      }
    }
  }
}
sv_matrix_pheno_lin_gene = read.csv("~/Documents/PHD/SV/SV_analysis/DR/sv_matrix_pheno_lin_gene.csv")

#There may be genes that have >1 so we will have to change them to 1
for (i in colnames(sv_matrix_pheno_lin_gene)) {
  print(which(colnames(sv_matrix_pheno_lin_gene)==i))
  if (length(print(table(sv_matrix_pheno_lin_gene[,i]>1))) == 2) {
    for (n in 1:nrow(sv_matrix_pheno_lin_gene)) {
      if (is.na(sv_matrix_pheno_lin_gene[n,i])) {next}
      if (sv_matrix_pheno_lin_gene[n,i] > 1) {
        sv_matrix_pheno_lin_gene[n,i] = 1
      }
    }
  }
}

sv_matrix_pheno_lin_gene_final = cbind(sv_matrix_pheno_lin[,c(1:22,ncol(sv_matrix_pheno_lin))],sv_matrix_pheno_lin_gene)

sv_matrix_pheno_lin_gene_final = merge(sv_matrix_pheno_lin_gene_final[,-which(names(sv_matrix_pheno_lin_gene_final) %in% c("X"))],snps_pca,by.x="row.names",by.y=1)%>% column_to_rownames(var="Row.names")
#sv_matrix_pheno_lin_gene_final = merge(sv_matrix_pheno_lin_gene_final,sv_lin,by.x=1,by.y=1) %>% column_to_rownames(var="Row.names")

#We now find genes with SVs associated to DR using logistic regression
dr_gene_pval_amk = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_amk = get_pvalues(drug = "amikacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_amk)
filter_gene(dr_gene_pval_amk,"amikacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_bdq = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_bdq = get_pvalues(drug = "bedaquiline",sv_matrix_pheno_lin_gene_final,dr_gene_pval_bdq)
filter_gene(dr_gene_pval_bdq,"bedaquiline",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_cap = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_cap = get_pvalues(drug = "capreomycin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_cap)
filter_gene(dr_gene_pval_cap,"capreomycin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_cip = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_cip = get_pvalues(drug = "ciprofloxacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_cip)
filter_gene(dr_gene_pval_cip,"ciprofloxacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_cfz = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_cfz = get_pvalues(drug = "clofazimine",sv_matrix_pheno_lin_gene_final,dr_gene_pval_cfz)
filter_gene(dr_gene_pval_cfz,"clofazimine",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_cs = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_cs = get_pvalues(drug = "cycloserine",sv_matrix_pheno_lin_gene_final,dr_gene_pval_cs)
filter_gene(dr_gene_pval_cs,"cycloserine",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_dlm = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_dlm = get_pvalues(drug = "delamanid",sv_matrix_pheno_lin_gene_final,dr_gene_pval_dlm)
filter_gene(dr_gene_pval_dlm,"delamanid",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_emb = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_emb = get_pvalues(drug = "ethambutol",sv_matrix_pheno_lin_gene_final,dr_gene_pval_emb)
filter_gene(dr_gene_pval_emb,"ethambutol",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_eta = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_eta = get_pvalues(drug = "ethionamide",sv_matrix_pheno_lin_gene_final,dr_gene_pval_eta)
filter_gene(dr_gene_pval_eta,"ethionamide",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_gat = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_gat = get_pvalues(drug = "gatifloxacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_gat)
filter_gene(dr_gene_pval_gat,"gatifloxacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_inh = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_inh = get_pvalues(drug = "isoniazid",sv_matrix_pheno_lin_gene_final,dr_gene_pval_inh)
filter_gene(dr_gene_pval_inh,"isoniazid",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_kan = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_kan = get_pvalues(drug = "kanamycin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_kan)
filter_gene(dr_gene_pval_kan,"kanamycin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_lev = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_lev = get_pvalues(drug = "levofloxacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_lev)
filter_gene(dr_gene_pval_lev,"levofloxacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_lzd = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_lzd = get_pvalues(drug = "linezolid",sv_matrix_pheno_lin_gene_final,dr_gene_pval_lzd)
filter_gene(dr_gene_pval_lzd,"linezolid",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_pas = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_pas = get_pvalues(drug = "para.aminosalicylic_acid",sv_matrix_pheno_lin_gene_final,dr_gene_pval_pas)
filter_gene(dr_gene_pval_pas,"para.aminosalicylic_acid",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_mfx = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_mfx = get_pvalues(drug = "moxifloxacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_mfx)
filter_gene(dr_gene_pval_mfx,"moxifloxacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_ofx = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_ofx = get_pvalues(drug = "ofloxacin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_ofx)
filter_gene(dr_gene_pval_ofx,"ofloxacin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_pza = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_pza = get_pvalues(drug = "pyrazinamide",sv_matrix_pheno_lin_gene_final,dr_gene_pval_pza)
filter_gene(dr_gene_pval_pza,"pyrazinamide",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_rbt = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_rbt = get_pvalues(drug = "rifabutin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_rbt)
filter_gene(dr_gene_pval_rbt,"rifabutin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_rif = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_rif = get_pvalues(drug = "rifampicin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_rif)
filter_gene(dr_gene_pval_rif,"rifampicin",marin,peppe,sv_matrix_pheno_lin_gene_final)

dr_gene_pval_str = data.frame(Drug=as.character(),SV=as.character(),Pval=as.numeric(),stringsAsFactors = F)
dr_gene_pval_str = get_pvalues(drug = "streptomycin",sv_matrix_pheno_lin_gene_final,dr_gene_pval_str)
filter_gene(dr_gene_pval_str,"streptomycin",marin,peppe,sv_matrix_pheno_lin_gene_final)

#Now we also find SVs using xgboost and keep SVs that are found by both xgboost and logistic regression methods to keep only confident SVs

#AMK
amk_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="amikacin")
#DR_pred accounts for 71.8% importance, all SVs are < 1%
amk_feat_gene = amk_feat_gene[c(118),]
amk_feat_gene = amk_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$amikacin)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$amikacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$amikacin)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$amikacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$amikacin)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$amikacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$amikacin)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$amikacin)[[2]], 0)
)
amk_feat_gene$Drug = "AMK"
amk_plot_gene = plot_features(amk_feat_gene,col="#e7298a",tit="AMK",wid=0.2,maxplot = 0.11)

#BDQ
bdq_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="bedaquiline")
#DR_pred accounts for 21.6% importance

#CAP
cap_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="capreomycin")
#DR_pred accounts for 59.6% importance all SVs are < 1%
cap_feat_gene = cap_feat_gene[c(9,52),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$capreomycin)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$capreomycin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$capreomycin)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$capreomycin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$capreomycin)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$capreomycin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$capreomycin)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$capreomycin)[[2]], 0)
)
cap_feat_gene$Drug = "CAP"
cap_plot_gene = plot_features(cap_feat_gene,col="#1f77b4",tit="CAP",wid=0.4,maxplot = 0.6)


#CIP
cip_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="ciprofloxacin")
#DR_pred accounts for 79% importance. No SVs that passed the filters
cip_feat_gene = cip_feat_gene[c(9,52),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ciprofloxacin)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ciprofloxacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ciprofloxacin)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$ciprofloxacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ciprofloxacin)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ciprofloxacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ciprofloxacin)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$ciprofloxacin)[[2]], 0)
)
cip_feat_gene$Drug = "CIP"
cip_plot_gene = plot_features(cip_feat_gene,col="#1f77b4",tit="CIP",wid=0.4,maxplot = 0.6)

#CFZ
cfz_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="clofazimine")
#DR_pred accounts for 2.7% importance
cfz_feat_gene = cfz_feat_gene[c(30,53,123),]
new_row <- data.frame(Feature = "Rv1949c", Gain = 0, Cover = 0, Frequency = 0)
cfz_feat_gene = rbind(cfz_feat_gene,new_row)
cfz_feat_gene = cfz_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$clofazimine)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$clofazimine)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$clofazimine)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$clofazimine)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$clofazimine)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$clofazimine)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$clofazimine)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$clofazimine)[[2]], 0)
)
cfz_feat_gene$Drug = "CFZ"
cfz_plot_gene = plot_features(cfz_feat_gene,col="#9467bd",wid=0.2,tit="CFZ",maxplot = 1.17)

#CS
cs_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="cycloserine")
#DR_pred accounts for 0% importance, SVs did not pass filters

#DLM
#DR_pred accounts for 3.4% importance
dlm_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="delamanid")
dlm_feat_gene = dlm_feat_gene[c(1,12,21,28,48,95),]
new_row <- data.frame(Feature = "ppsA", Gain = 0, Cover = 0, Frequency = 0)
dlm_feat_gene = rbind(dlm_feat_gene,new_row)
dlm_feat_gene = dlm_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$delamanid)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$delamanid)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$delamanid)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$delamanid)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$delamanid)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$delamanid)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$delamanid)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$delamanid)[[2]], 0)
)
dlm_feat_gene$Drug = "DLM"
dlm_plot_gene = plot_features(dlm_feat_gene,col="#bcbd22",wid=0.35,tit="DLM",maxplot = 13.1)

#EMB
#DR_pred accounts for 69.8% importance
emb_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="ethambutol")
emb_feat_gene = emb_feat_gene[c(71,376),]
emb_feat_gene = emb_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethambutol)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethambutol)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ethambutol)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$ethambutol)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethambutol)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethambutol)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ethambutol)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$ethambutol)[[2]], 0)
)
emb_feat_gene$Drug = "EMB"
emb_plot_gene = plot_features(emb_feat_gene,col="#17becf",wid=0.3,tit="EMB",maxplot = 0.14)

#ETA
#DR_pred accounts for 49.5% importance
eta_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="ethionamide")
eta_feat_gene = eta_feat_gene[c(74,236),]
eta_feat_gene = eta_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethionamide)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethionamide)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ethionamide)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$ethionamide)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethionamide)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ethionamide)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ethionamide)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$ethionamide)[[2]], 0)
)
eta_feat_gene$Drug = "ETA"
eta_plot_gene = plot_features(eta_feat_gene,col="#3288bd",wid=0.4,tit="ETA",maxplot = 0.19)

#GAT
#DR_pred accounts for 0% importance
gat_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="gatifloxacin")

#INH
#DR pred accounts for 84.9%
inh_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="isoniazid")
inh_feat_gene = inh_feat_gene[c(27,219,171,217,223,244),]
new_row <- data.frame(Feature = "hsdM", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0221", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0740", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0739", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0741", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv1776c", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv2767c", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv2766c", Gain = 0, Cover = 0, Frequency = 0)
inh_feat_gene = rbind(inh_feat_gene,new_row)
inh_feat_gene = inh_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$isoniazid)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$isoniazid)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$isoniazid)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$isoniazid)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$isoniazid)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$isoniazid)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$isoniazid)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$isoniazid)[[2]], 0)
)
inh_feat_gene$Drug = "INH"
inh_plot_gene = plot_features(inh_feat_gene,col="#66c2a5",wid=0.4,tit="INH",maxplot = 0.14)

#KAN
kan_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="kanamycin")
#DR pred accounts for 67.7%

#LEV
lev_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="levofloxacin")
#DR pred accounts for 76.3%
lev_feat_gene = lev_feat_gene[c(163),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$levofloxacin)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$levofloxacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$levofloxacin)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$levofloxacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$levofloxacin)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$levofloxacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$levofloxacin)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$levofloxacin)[[2]], 0)
)
lev_feat_gene$Drug = "LEV"
lev_plot_gene = plot_features(lev_feat_gene,col="#7570b3",wid=0.35,tit="LEV",maxplot = 0.07)

#LZD
lzd_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="linezolid")
#DR pred accounts for 24.6%
lzd_feat_gene = lzd_feat_gene[c(104),] %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$linezolid)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$linezolid)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$linezolid)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$linezolid)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$linezolid)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$linezolid)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$linezolid)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$linezolid)[[2]], 0)
)
lzd_feat_gene$Drug = "LZD"
lzd_plot_gene = plot_features(lzd_feat_gene,col="#a6d854",wid=0.35,tit="LZD",maxplot = 0.22)

#MFX
mfx_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="moxifloxacin")
#DR pred accounts for 61%

#OFX
ofx_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="ofloxacin")
#DR pred accounts for 66.2%
ofx_feat_gene = ofx_feat_gene[c(171),]
ofx_feat_gene = ofx_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ofloxacin)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ofloxacin)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ofloxacin)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$ofloxacin)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ofloxacin)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$ofloxacin)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$ofloxacin)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$ofloxacin)[[2]], 0)
)
ofx_feat_gene$Drug = "OFX"
ofx_plot_gene = plot_features(ofx_feat_gene,col="#e377c2",wid=0.35,tit="OFX",maxplot = 0.14)

#PAS
pas_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="para.aminosalicylic_acid")
#DR pred accounts for 15.2%

#PZA
pza_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="pyrazinamide")
#DR pred accounts for 69.7 %
pza_feat_gene = pza_feat_gene[c(7,119,128,245),]
new_row <- data.frame(Feature = "hupB", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv1506c", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0147", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0221", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0150c", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0149", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "mce3F", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "lprM", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "mce3D", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "cyp130", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
new_row <- data.frame(Feature = "Rv0310c", Gain = 0, Cover = 0, Frequency = 0)
pza_feat_gene = rbind(pza_feat_gene,new_row)
pza_feat_gene = pza_feat_gene %>% rowwise() %>% mutate(
  freq_r = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$pyrazinamide)) >= 1,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$pyrazinamide)[[1]], 0),
  freq_r_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$pyrazinamide)) >= 1,
                      table(sv_matrix_pheno_lin_gene_final$pyrazinamide)[[1]], 0),
  freq_s = ifelse(length(table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$pyrazinamide)) >= 2,
                  table(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[Feature]] == 1,]$pyrazinamide)[[2]], 0),
  freq_s_tot = ifelse(length(table(sv_matrix_pheno_lin_gene_final$pyrazinamide)) >= 2,
                      table(sv_matrix_pheno_lin_gene_final$pyrazinamide)[[2]], 0)
)
pza_feat_gene$Drug = "PZA"
pza_plot_gene = plot_features(pza_feat_gene,col="#b3b3b3",wid=0.35,tit="PZA",maxplot = 1.05)

#RBT
rbt_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="rifabutin")
#DR pred accounts for 0%

#RIF
rif_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="rifampicin")
#DR pred accounts for 88.1 %

#STR
str_feat_gene = get_importance_matrix(sv_matrix_pheno_lin_gene_final,tbprofiler,drug="streptomycin")
#DR pred accounts for 68.2 %

#WHO v2 Mutations Importance
who_plot_gene = data.frame(Importance=c(71.8,21.6,59.6,79.8,2.5,0,3.4,69.8,49.5,0,84.9,67.7,76.3,24.6,61,66.2,15.2,69.7,0,88.1,68.2),
                      Drug=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MFX","OFX","PAS","PZA",
                             "RBT","RIF","STR")) %>%
  ggplot(.,aes(x=Importance,y=reorder(Drug,Importance))) +
  geom_segment(aes(xend = 0, yend = Drug, color = Drug)) +
  geom_point(size = 5, aes(color = Drug)) +
  labs(x = "TB-Profiler mutation list Importance (%)") +
  theme_pubr() +
  labs(y="") +
  theme(plot.title = element_text(size=35,hjust = 0.5)) +
  scale_color_manual(values = c("AMK" = "#e7298a","BDQ" = "#ff7f0e","CAP" = "#1f77b4","CIP" = "#d62728","CFZ"="#9467bd",
                                "CS" = "#8c564b","DLM" = "#bcbd22","EMB"="#17becf","ETA" = "#3288bd","GAT" = "#ffd92f",
                                "INH" = "#66c2a5","KAN" = "#d95f02","LEV" = "#7570b3","LZD" ="#a6d854","MFX"="#e7298a",
                                "OFX"="#e377c2","PAS" = "#e78ac3","PZA" = "#b3b3b3","RBT"="#7f7f7f","RIF"="#fdae61","STR"="#2ca02c")) +
  theme(legend.position = "none",plot.title = element_text(size=45,hjust = 0.5),axis.text = element_text(size = 20), axis.ticks = element_line(linewidth = 3),
        legend.text = element_text(size = 10),axis.title.x = element_text(size=20)) +
  scale_x_continuous(limits = c(0, 100))

#Final SV Table with p-values (adjusted significant + non-significant but in associated genes)
pval_tab_gene = do.call("rbind",list(amk_feat_gene,cfz_feat_gene,cap_feat_gene,dlm_feat_gene,emb_feat_gene,eta_feat_gene,inh_feat_gene,lev_feat_gene,lzd_feat_gene,pza_feat_gene))
pval_tab_gene = pval_tab_gene[,c(9,1,5,6,7,8)]
for (i in 1:nrow(pval_tab_gene)) {
  # For Resistant isolates
  drug = pval_tab_gene[i,]$Drug
  drug_ac = data.frame(short=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MXF",
                               "OFX","PAS","PZA","RBT","RIF","STR"),
                       long=c("amikacin","bedaquiline","capreomycin","ciprofloxacin","clofazimine","cycloserine","delamanid",
                              "ethambutol","ethionamide","gatifloxacin","isoniazid","kanamycin","levofloxacin","linezolid","moxifloxacin",
                              "ofloxacin","para.aminosalicylic_acid","pyrazinamide","rifabutin","rifampicin","streptomycin"),stringsAsFactors = F)
  drug = drug_ac[drug_ac$short==drug,]$long
  sv_col = pval_tab_gene[i,]$Feature
  resistant_isolates <- sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[drug]] == "R", sv_col]
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[drug]] == "S", sv_col]
  pval_tab_gene[i,]$freq_r = sum(resistant_isolates == 1, na.rm = TRUE)
  pval_tab_gene[i,]$freq_s = sum(susceptible_isolates == 1, na.rm = TRUE)
  pval_tab_gene[i,]$freq_r_tot = sum(!is.na(resistant_isolates))
  pval_tab_gene[i,]$freq_s_tot = sum(!is.na(susceptible_isolates))
}
dr_gene_all = do.call("rbind",list(dr_gene_pval_amk%>% mutate(Drug="AMK"),dr_gene_pval_cap%>% mutate(Drug="CAP"),dr_gene_pval_cfz%>% mutate(Drug="CFZ"),
                                   dr_gene_pval_dlm%>% mutate(Drug="DLM"),dr_gene_pval_emb%>% mutate(Drug="EMB"),dr_gene_pval_eta%>% mutate(Drug="ETA"),
                                   dr_gene_pval_inh%>% mutate(Drug="INH"),dr_gene_pval_lev%>% mutate(Drug="LEV"),dr_gene_pval_lzd%>% mutate(Drug="LZD"),dr_gene_pval_ofx%>% mutate(Drug="OFX"),
                                   dr_gene_pval_pza%>% mutate(Drug="PZA")))
pval_tab_gene = merge(pval_tab_gene,dr_gene_all,by.x=c(1,2),by.y=c(1,2),all.x=T)
#pval_tab_gene_nosig = pval_tab_gene %>% filter(Adj_pval>0.05)
#pval_tab_gene_sig = pval_tab_gene %>% filter(Adj_pval<0.05)
#pval_tab_gene = pval_tab_gene %>% mutate(Adj_pval = pmin(Pval * Total_SVs / seq_along(Pval), 1))
pval_tab_gene=pval_tab_gene[!is.na(pval_tab_gene$Adj_pval),]
write.table(pval_tab_gene,"~/Documents/PHD/SV/SV_analysis/DR/gene_dr_pval.csv",sep=",")
#write.tab_genele(pval_tab_gene_nosig,"~/Documents/PHD/SV/SV_analysis/DR/sv_dr_pval_nosig.csv",sep=",")

#TABLE part of the plot
#First, we rbind all of the features together to plot them all in the same table
all_feat_gene = do.call("rbind",list(amk_feat_gene,cap_feat_gene,cfz_feat_gene,dlm_feat_gene,emb_feat_gene,eta_feat_gene,inh_feat_gene,lev_feat_gene,lzd_feat_gene,ofx_feat_gene,pza_feat_gene))
all_feat_gene = all_feat_gene[,c(9,1,5,6,7,8)]
for (i in 1:nrow(all_feat_gene)) {
  # For Resistant isolates
  drug = all_feat_gene[i,]$Drug
  drug_ac = data.frame(short=c("AMK","BDQ","CAP","CIP","CFZ","CS","DLM","EMB","ETA","GAT","INH","KAN","LEV","LZD","MXF",
                               "OFX","PAS","PZA","RBT","RIF","STR"),
                       long=c("amikacin","bedaquiline","capreomycin","ciprofloxacin","clofazimine","cycloserine","delamanid",
                              "ethambutol","ethionamide","gatifloxacin","isoniazid","kanamycin","levofloxacin","linezolid","moxifloxacin",
                              "ofloxacin","para.aminosalicylic_acid","pyrazinamide","rifabutin","rifampicin","streptomycin"),stringsAsFactors = F)
  drug = drug_ac[drug_ac$short==drug,]$long
  sv_col = all_feat_gene[i,]$Feature
  resistant_isolates <- sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[drug]] == "R", sv_col]
  resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[[drug]] == "S", sv_col]
  susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
  all_feat_gene$freq_r <- as.numeric(all_feat_gene$freq_r)
  all_feat_gene$freq_s <- as.numeric(all_feat_gene$freq_s)
  all_feat_gene$freq_r_tot <- as.numeric(all_feat_gene$freq_r_tot)
  all_feat_gene$freq_s_tot <- as.numeric(all_feat_gene$freq_s_tot)
  all_feat_gene[i,]$freq_r = resistant_freq
  all_feat_gene[i,]$freq_s = susceptible_freq
  all_feat_gene[i,]$freq_r_tot = sum(!is.na(resistant_isolates))
  all_feat_gene[i,]$freq_s_tot = sum(!is.na(susceptible_isolates))
}
all_feat_gene=all_feat_gene %>%
  left_join(pval_tab_gene, by=c("Drug", "Feature"))
all_feat_gene = all_feat_gene[order(all_feat_gene$Adj_pval)[1:20],]
all_feat_gene = all_feat_gene %>%
  left_join(genes_names %>% dplyr::select(Name, V4), by = c("Feature" = "Name"))# %>%
  #rename(gene_name = V4)
all_feat_gene = all_feat_gene[,c(1,2,14,3,5,4,6)]
colnames(all_feat_gene) = c("Drug","Gene","Start Position","Resistant Frequency","Susceptible Frequency","Total Resistant","Total Susceptible")
all_feat_gene$`Susceptible Frequency` <- sprintf("%.10f",all_feat_gene$`Susceptible Frequency`)  # Show 10 decimal places
all_feat_gene$`Susceptible Frequency` = as.numeric(all_feat_gene$`Susceptible Frequency`)

drug_color_scale_gene <- formatter("span",
                              style = x ~ style( display="block",
                                                 "border-radius" = "0px",
                                                 "background-color" = case_when(all_feat_gene$Drug=="AMK" ~ "#e7298a",
                                                                                all_feat_gene$Drug=="CAP" ~ "#1f77b4",
                                                                                all_feat_gene$Drug=="CFZ" ~ "#9467bd",
                                                                                all_feat_gene$Drug=="DLM" ~ "#bcbd22",
                                                                                all_feat_gene$Drug=="EMB" ~ "#17becf",
                                                                                all_feat_gene$Drug=="ETA" ~ "#3288bd",
                                                                                all_feat_gene$Drug=="INH" ~ "#66c2a5",
                                                                                all_feat_gene$Drug=="LEV" ~ "#7570b3",
                                                                                all_feat_gene$Drug=="LZD" ~ "#a6d854",
                                                                                all_feat_gene$Drug=="OFX" ~ "#e377c2",
                                                                                all_feat_gene$Drug=="PZA" ~ "#b3b3b3",
                                                                                all_feat_gene$Drug=="RIF" ~ "#fdae61"),
                                                 "color" = "black",
                                                 "border-radius" = "4px"))
tab_gene = formattable(all_feat_gene, list("Resistant Frequency" = r_color_scale,
                                 "Susceptible Frequency" = s_color_scale,
                                 "Drug" = drug_color_scale_gene)) %>%
  as.htmlwidget(width="100%") %>%
  prependContent(tags$style("table,td,tr,th { border: 1px solid black !important;}","td, th { text-align: left !important; padding-left: 5px; }" ))
export_formattable(tab_gene,"~/Downloads/tab_gene_20.png")

#Fig 5!
pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/Fig5_gene_20_tbp.pdf",width = 20,height = 35)
plot_grid(plot_grid(amk_plot_gene,cap_plot_gene,cfz_plot_gene,dlm_plot_gene,emb_plot_gene,eta_plot_gene,inh_plot_gene,lev_plot_gene,lzd_plot_gene,
                                     ofx_plot_gene,pza_plot_gene,who_plot_gene,nrow = 3,ncol=4,labels = c("a",""),label_size = 40),
          rasterGrob(png::readPNG("~/Downloads/tab_gene_20.png")),labels = c("","b"),nrow=2,ncol=1,label_size = 40,rel_heights = c(1.2,1))
dev.off()

#Example of finding to which drug the SV is associated in the case of MDR
which(colnames(sv_matrix_pheno_lin_gene_final)=="Rv2522c")
sv_matrix_pheno_lin_mod = sv_matrix_pheno_lin_gene_final[!is.na(sv_matrix_pheno_lin_gene_final$pyrazinamide) & !is.na(sv_matrix_pheno_lin_gene_final$delamanid), c(18,7,3769,4209,4210,4211,4212,4213)]
sv_matrix_final <- sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
summary(logistf(Rv2522c ~ pyrazinamide + delamanid + V2 + V3 + V4 + V5, data = sv_matrix_final, family = binomial))


#Command to find whether unexplained resistance in SVs that we've found associated to DR is due to SVs in known drug resistance-conferring genes, and not due to the new SVs
sv_matrix_pheno_lin_gene_final[rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$Rv1998c==1 & sv_matrix_pheno_lin_gene_final$clofazimine=="R" & !is.na(sv_matrix_pheno_lin_gene_final$clofazimine) & !is.na(sv_matrix_pheno_lin_gene_final$Rv1998c),])[!(rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$Rv1998c==1 & sv_matrix_pheno_lin_gene_final$clofazimine=="R" & !is.na(sv_matrix_pheno_lin_gene_final$clofazimine) & !is.na(sv_matrix_pheno_lin_gene_final$Rv1998c),]) %in% tbprofiler[tbprofiler$drugs=="clofazimine",]$sample)],"atpE"]

#VALIDATION OF ASSOCIATIONS THROUGH PERMUTATION TEST
pval_tab_gene2= pval_tab_gene %>% mutate(Drug=case_when(Drug=="AMK" ~ "amikacin",
                                               Drug=="BDQ" ~ "bedaquiline",
                                               Drug=="CFZ" ~ "clofazimine",
                                               Drug=="DLM" ~ "delamanid",
                                               Drug=="CAP" ~ "capreomycin",
                                               Drug=="EMB" ~ "ethambutol",
                                               Drug=="ETA" ~ "ethionamide",
                                               Drug=="INH" ~ "isoniazid",
                                               Drug=="KAN" ~ "kanamycin",
                                               Drug=="LEV" ~ "levofloxacin",
                                               Drug=="LZD" ~ "linezolid",
                                               Drug=="OFX" ~ "ofloxacin",
                                               Drug=="PZA" ~ "pyrazinamide",
                                               Drug=="STR" ~ "streptomycin",))
set.seed(123)
perm_df_gene = data.frame(sv=as.character(),perm=as.character(),stringsAsFactors = F)
for (i in 1:nrow(pval_tab)){
  perm_df_gene=rbind(perm_df_gene,data.frame(sv=pval_tab_gene2[i,1],perm=perm(pval_tab_gene2[i,1],sv_matrix_pheno_lin_gene_final,pval_tab_gene2[i,2],pval_tab_gene2[i,7])))
}

#hsdM frequency figure
sv_col <- "X3069892.DEL.102"

results_hsdm <- data.frame(Drug = character(), 
                             Resistance_Status = character(), 
                             Frequency = numeric(),
                             stringsAsFactors = FALSE)

for (drug in drug_cols) {
  # For Resistant isolates
  resistant_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "R", sv_col]
  resistant_freq <- sum(resistant_isolates == 1, na.rm = TRUE) / sum(!is.na(resistant_isolates))
  # Append result for resistant
  results_hsdm <- rbind(results_hsdm, data.frame(Drug = drug, 
                                                     Resistance_Status = "Resistant", 
                                                     Frequency = resistant_freq))
  # For Susceptible isolates
  susceptible_isolates <- sv_matrix_pheno_lin[sv_matrix_pheno_lin[[drug]] == "S", sv_col]
  susceptible_freq <- sum(susceptible_isolates == 1, na.rm = TRUE) / sum(!is.na(susceptible_isolates))
  # Append result for susceptible
  results_hsdm <- rbind(results_hsdm, data.frame(Drug = drug, 
                                                     Resistance_Status = "Susceptible", 
                                                     Frequency = susceptible_freq))
}
hsdm_freq = ggplot(results_hsdm%>%mutate(Drug=case_when(Drug=="amikacin"~"AMK",Drug=="bedaquiline"~"BDQ",Drug=="clofazimine"~"CFZ",Drug=="delamanid"~"DLM",
                                                            Drug=="ethambutol"~"EMB",Drug=="ethionamide"~"ETA",Drug=="isoniazid"~"INH",Drug=="kanamycin"~"KAN",
                                                            Drug=="para.aminosalicylic_acid"~"PAS",Drug=="pyrazinamide"~"PZA",Drug=="rifampicin"~"RIF",
                                                            Drug=="streptomycin"~"STR",Drug=="rifabutin"~"RBT",Drug=="ofloxacin"~"OFX",Drug=="moxifloxacin"~"MFX",
                                                            Drug=="linezolid"~"LZD",Drug=="levofloxacin"~"LEV",Drug=="gatifloxacin"~"GAT",Drug=="cycloserine"~"CS",
                                                            Drug=="capreomycin"~"CAP",Drug=="ciprofloxacin"~"CIP")), aes(x = Drug, y = Frequency, color = Resistance_Status)) +
  geom_point(position = position_dodge(width = 0.5), size = 3,shape=18) +  # Dodge the points
  theme_pubr() +
  labs(x = "", y = "Frequency of 3069892.DEL.102", 
       color = "Resistance Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(aes(xintercept = as.numeric(factor(Drug))), color = "gray90", linetype = "dashed") +
  scale_color_manual(values = c("#D43925","#2376B2"))

pdf("~/Documents/PHD/SV/SV_analysis/figures/r_figures/Fig6_hsdM.pdf",width = 6,height = 4)
hsdm_freq
dev.off()

#HOMOPLASY
essentially_44k = data.frame(Name="",svpop="",lin="",stringsAsFactors = F)
for (i in 23:ncol(sv_matrix_pheno_lin_gene_final[,23:ncol(sv_matrix_pheno_lin_gene_final)])) {
  print(i)
  essentially_44k = rbind(essentially_44k,data.frame(Name=colnames(sv_matrix_pheno_lin_gene_final)[i],
                                                     svpop=sum(sv_matrix_pheno_lin_gene_final[!is.na(sv_matrix_pheno_lin_gene_final[,i]),i]),
                                                     lin=paste(unique(sv_lin[sv_lin$ID %in% rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final[,i]==1,]),]$Lineage),collapse=",")))
}
essentially_44k=merge(essentially_44k,essentially,by.x=1,by.y=2)
essentially_44k=essentially_44k[,c(1,2,3,9,10)]
for (i in 1:nrow(marin)) {
  essentially_44k = essentially_44k %>% filter(V4-3000>marin[i,2] &V5+3000>marin[i,2] | V4+3000<marin[i,1] &V5-3000<marin[i,1])
}
essentially_44k$svpop=as.numeric(essentially_44k$svpop)
essentially_44k_pos = essentially_44k[essentially_44k$svpop>99 & grepl(",",essentially_44k$lin.x),]

four_tree = read.newick("~/Documents/PHD/SV/SV_analysis/DR/44k.vft.ml.tree")
four_tree = four_tree %>% drop.tip(.,four_tree$tip.label[!(four_tree$tip.label %in% rownames(sv_matrix_pheno_lin))])

tree_data44k <- four_tree %>% 
  as_tibble() %>% 
  mutate(ctpG = case_when(label %in% rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$ctpG==1,])~ "SV in ctpG",
                            TRUE~" "),
         glnA3 = case_when(label %in% rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$glnA3==1,])~ "SV in glnA3", 
                          TRUE~" "),
         ethA = case_when(label %in% rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$ethA==1,])~ "SV in ethA",
                          TRUE~" "),
         pncA = case_when(label %in% rownames(sv_matrix_pheno_lin_gene_final[sv_matrix_pheno_lin_gene_final$pncA==1,])~ "SV in pncA", 
                           TRUE~" "))

ggtree(four_tree,#%>% drop.tip(.,four_tree$tip.label[four_tree$tip.label %in% rownames(sv_matrix_pheno_lin_gene_final[is.na(sv_matrix_pheno_lin_gene_final$Rv3177),])]),
       layout = "circular") %<+% sv_lin + geom_tippoint(aes(colour=Lineage), size = 0.25)+ scale_color_manual(values = lin_color2) +#%>% drop.tip(.,four_tree$tip.label[four_tree$tip.label %in% sv_lin[!grepl("^1.2.1",sv_lin$X),]$ID]),layout = "circular") + 
  geom_fruit(data=tree_data44k[!is.na(tree_data44k$label),],aes(y = label, fill=pncA),height=5.5,width = 0.003,geom='geom_tile') +
  geom_fruit(data=tree_data44k[!is.na(tree_data44k$label),],aes(y = label, fill=ctpG),height=5.5,width = 0.003,geom='geom_tile') +
  geom_fruit(data=tree_data44k[!is.na(tree_data44k$label),],aes(y = label, fill=ethA),height=5.5,width = 0.003,geom='geom_tile') +
  geom_fruit(data=tree_data44k[!is.na(tree_data44k$label),],aes(y = label, fill=glnA3),height=5.5,width = 0.003,geom='geom_tile') +
  scale_fill_manual(values = c("SV in pncA"="#b3b3b3","SV in ctpG"="#FF7F00","SV in ethA"="#3288bd"," "="white","SV in glnA3" = "#E31A1C")) + theme(legend.title = element_blank())











#Extra code which didn't work
# Train the model using Random Forest
model <- train(
  amikacin ~ . + lineage,
  data = sv_matrix_final,  # Keep NAs as they are
  method = "rf",
  trControl = train_control,
  na.action = na.roughfix  # Allow NAs to be passed through to the model
)

# Get and print variable importance
importance <- varImp(model, scale = TRUE)
print(importance)


#LRS-DR -- in case I cant find enough SVs in SRS
long_meta = read.table("~/Downloads/SV_metadata1.csv",sep=",",header = T)
long_meta = long_meta[,-c(1,3:12,14:16)]
long_meta = long_meta[rowSums(is.na(long_meta)) < ncol(long_meta)-1, ]
sv_pca_meta = merge(long_meta,sv_pca,by.x=1,by.y="row.names")

control <- rfeControl(functions = rfFuncs,  # Use random forest functions
                      method = "cv",        # Cross-validation
                      number = 5)          # 10-fold cross-validation
sv_pca_meta_ami = sv_pca_meta[!is.na(sv_pca_meta$AMI_resistance),c(12,21:ncol(sv_pca_meta))]

results <- rfe(
  x = sv_pca_meta_ami[, -which(names(sv_pca_meta_ami) %in% c("AMI_resistance", "Lineage"))], 
  y = sv_pca_meta_ami$AMI_resistance,
  sizes = c(1:10, 15, 20),  # Specify the subset sizes to consider
  rfeControl = control,
  method = "rf",  # Use random forest for RFE
  trControl = trainControl(method = "cv", number = 5)
)

#GLM
drug="amikacin"
sv_matrix_pheno_lin_mod = sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin$amikacin),c(1,23:ncol(sv_matrix_pheno_lin))] %>% mutate_all(as.factor)
sv_matrix_final = sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
sv_matrix_final <- sv_matrix_final %>% 
  mutate(DR_prediction = ifelse(rownames(sv_matrix_final) %in% tbprofiler[tbprofiler$drugs == drug, ]$sample, 1, 0))
results <- lapply(names(sv_matrix_final[,-c(1,ncol(sv_matrix_final))]), function(sv) {
  formula <- as.formula(paste(drug, "~", sv, "+ X + DR_prediction"))
  #print(which(colnames(sv_matrix_final)==sv))
  glm(formula, data = sv_matrix_final, family = binomial)
})

# Extract p-values
pvals <- sapply(results, function(model) summary(model)$coefficients[2,4])

# Adjust for multiple testing (e.g., Bonferroni or FDR)
adjusted_pvals <- p.adjust(pvals, method = "fdr")

# View significant results
names(adjusted_pvals) = colnames(sv_matrix_final[,-c(1,ncol(sv_matrix_final))])
significant_sv <- names(adjusted_pvals)[adjusted_pvals < 0.05]
significant_sv


sv_matrix_pheno_lin_mod = sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin$bedaquiline),c(2,23:ncol(sv_matrix_pheno_lin))] %>% mutate_all(as.factor)
sv_matrix_final = sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
sv_matrix_final_tb=sv_matrix_final %>% mutate(DR_pred = ifelse(rownames(sv_matrix_final) %in% tbprofiler[tbprofiler$drugs=="bedaquiline",]$sample,1,0)) 
sv_matrix_final_tb_num=sv_matrix_final_tb %>% mutate_all(as.numeric) %>% mutate(bedaquiline=ifelse(bedaquiline==1,1,0))
sv_matrix_final_tb_num_marin =filter_columns_by_position(sv_matrix_final_tb_num, marin)

for (i in 20:22) {
  print(colnames(sv_matrix_pheno_lin)[i])
  sv_matrix_pheno_lin_mod = sv_matrix_pheno_lin[!is.na(sv_matrix_pheno_lin[,i]),c(i,23:ncol(sv_matrix_pheno_lin))] %>% mutate_all(as.factor)
  sv_matrix_final = sv_matrix_pheno_lin_mod[, sapply(sv_matrix_pheno_lin_mod, function(x) !(is.factor(x) && length(levels(x)) < 2))]
  
  # Identify the SVs that are only present in resistant isolates and absent in susceptible isolates
  exclusive_svs <- apply(sv_matrix_final[,-1], 2, function(col) {
    # Check the presence of SV in resistant isolates
    sum_resistant <- sum(as.numeric(col[as.character(sv_matrix_final[,1]) == "R"]), na.rm = TRUE)
    sum_susceptible <- sum(as.numeric(col[as.character(sv_matrix_final[,1]) == "S"]), na.rm = TRUE)
    #print(names(col[as.character(sv_matrix_final[,1]) == "R" & !is.na(col) & col==1]))
    #print(tbprofiler[tbprofiler$prediction=="R" & tbprofiler$drug==colnames(sv_matrix_pheno_lin)[i],]$run)
    if (identical(names(col[as.character(sv_matrix_final[,1]) == "R" & !is.na(col) & as.numeric(col)==1]),character(0) )) {
      return(FALSE)
    }
    # Check if the SV is present in more than one resistant isolate and absent in all susceptible
    #if (sum_susceptible <= 0.5 * sum_resistant & sum_resistant > 3) {
    if (sum_resistant/(sum_susceptible + sum_resistant) > 0.7 & sum_resistant>1) {
      if (all(
        names(col[as.character(sv_matrix_final[,1]) == "R" & !is.na(col) & as.numeric(col)==1]) %in% tbprofiler[tbprofiler$prediction=="R" & tbprofiler$drug==colnames(sv_matrix_pheno_lin)[i],]$run)) {
        return(FALSE)
      } else {return(TRUE)}
    } else {
      return(FALSE)
    }
  })
  exclusive_svs_names <- names(exclusive_svs)[exclusive_svs]
  for (n in exclusive_svs_names) {
    print(paste0(n," ",table(sv_matrix_final[sv_matrix_final[,n]==1,1])))
    print(table(droplevels(sv_matrix_final[sv_matrix_final[,n]==1 & sv_matrix_final[,1]=="R",ncol(sv_matrix_final)])))
    print(table(droplevels(sv_matrix_final[sv_matrix_final[,n]==1 & sv_matrix_final[,1]=="S",ncol(sv_matrix_final)])))
  }
}

# Get the names of the SVs that meet the criteria
exclusive_svs_names <- names(exclusive_svs)[exclusive_svs]

# Output the result
exclusive_svs_names

# % inside each drug, I prefer the other fig, though
df_summary %>%
  filter(!is.na(Phenotype)) %>% group_by(Drug) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>% mutate(Drug = recode(Drug, !!!drug_acronyms)) %>%
  ggplot(., aes(x = "", y = Percentage, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # Black borders for slices
  coord_polar(theta = "y") +
  facet_wrap(~ Drug, ncol = 11) +  # Arrange in a grid, adjust ncol as needed
  theme_minimal(base_size = 15) +  # Minimalist theme with larger base font size
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text (categories around the pie)
    axis.text.y = element_blank(),  # Remove y-axis text (numbers around the pie)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 14),  # Facet labels
    legend.position = "bottom",  # Legend at the bottom
    legend.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_fill_manual(values = c("#D43925","#2376B2","#d0d6d6"))  # Choose a good palette for categorical data
labs(title = "Phenotypic Distribution per Drug") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),  # Add percentage labels inside the pie slices
            position = position_stack(vjust = 0.5), size = 3, color = "black")