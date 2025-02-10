library(stringr)
library(parallel)
library(dplyr)
library(tidyverse)
#
##FUNCTIONS
#
# Function to check overlap considering position, type, and size (for DELs)
check_overlap <- function(merged_sv, sample_sv) {
  pos_overlap <- (sample_sv$POS >= (merged_sv$POS - 50)) & (sample_sv$POS <= (merged_sv$POS + 50))
  type_match <- ((strsplit(as.character(sample_sv$ID),".",fixed = TRUE)[[1]][2] == strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][2])) | (grepl("INS", sample_sv$ID) & grepl("DUP", merged_sv$ID)) | ((grepl("DUP", sample_sv$ID) & grepl("INS", merged_sv$ID)))  # Match SV types (INS with INS, DEL with DEL, etc.)
  if (grepl("DEL", merged_sv$ID) & grepl("DEL", sample_sv$ID)) {
    merged_size <- as.numeric(strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][3])
    size_diff <- abs(merged_size - sample_sv$SIZE)
    return(any(pos_overlap & type_match & (size_diff <= 100)))
  } else {
    return(any(pos_overlap & type_match))
  }
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
      pos_col <- vcf[m,2]
      end_col <- strsplit(strsplit(as.character(vcf[m,8]),";",fixed = TRUE)[[1]][1],"=",fixed = T)[[1]][2]
    }
    else {
      id_col <- paste(vcf[m,2],vcf[m,3],sep = ".")
      pos_col <- vcf[m,2]
      end_col <- vcf[m,2]
    }
    
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
    if (same_type && abs(current_sv$END - next_sv$POS) <= 50) {
      current_sv$END <- max(current_sv$END, next_sv$END)
    } else if (is_del && abs(current_sv$END - next_sv$POS) <= 50 &&
               abs(as.numeric(strsplit(as.character(current_sv$ID),".",fixed = TRUE)[[1]][3]) - as.numeric(strsplit(as.character(next_sv$ID),".",fixed = TRUE)[[1]][3])) <= 100) {
      current_sv$END <- max(current_sv$END, next_sv$END)
    } else if (is_ins_dup && abs(current_sv$END - next_sv$POS) <= 50) {
      current_sv$END <- max(current_sv$END, next_sv$END)
    } else {
      merged_svs <- rbind(merged_svs, current_sv)
      current_sv <- next_sv
    }
  }
  
  merged_svs <- rbind(merged_svs, current_sv)
  return(merged_svs)
}

load("/data/projects/punim1637/Aleix/SV/SV_analysis/180924_newmatrix.RData")

sample_names <- str_replace(files_vcf,"_ref_short.vcf","")
sv_matrix <- matrix(0, nrow = length(files_vcf), ncol = nrow(merged_svs),
                    dimnames = list(sample_names, paste(merged_svs$ID, merged_svs$POS, merged_svs$END, sep=".")))

#process_sample <- function(i) {
#  sv_data <- process_mtb_vcf(tabledr[[i]])
#  sample_id <- sample_names[i]
#  
#  sv_data$POS <- as.numeric(as.character(sv_data$POS))
#  sv_data$END <- as.numeric(as.character(sv_data$END))
#  sv_data$SIZE <- as.numeric(as.character(sv_data$SIZE))
#  
#  for (j in seq_along(merged_svs$ID)) {
#    merged_sv <- merged_svs[j, ]
#    
#    if (any(sv_data$ID == "NaN" & sv_data$POS <= merged_sv$END & sv_data$END >= merged_sv$POS)) {
#      sv_matrix[sample_id, j] <- NA
#    } else {
#      overlap_result <- mapply(check_overlap, list(merged_sv), split(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS), seq(nrow(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS)))))
#      sv_matrix[sample_id, j] <- ifelse(any(overlap_result), 1, 0)
#    }
#  }
#}

for (i in 1:1) {
  print(i)
  print(sv_matrix)
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
      overlap_result <- mapply(check_overlap, list(merged_sv), split(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS), seq(nrow(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS)))))
      sv_matrix[sample_id, j] <- ifelse(any(overlap_result), 1, 0)
    }
  }
}

check_overlap <- function(merged_sv, sample_sv) {
  pos_overlap <- (sample_sv$POS >= (merged_sv$POS - 50)) & (sample_sv$POS <= (merged_sv$POS + 50))
  type_match <- ((strsplit(as.character(sample_sv$ID),".",fixed = TRUE)[[1]][2] == strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][2])) | (grepl("INS", sample_sv$ID) & grepl("DUP", merged_sv$ID)) | ((grepl("DUP", sample_sv$ID) & grepl("INS", merged_sv$ID)))  # Match SV types (INS with INS, DEL with DEL, etc.)
  if (grepl("DEL", merged_sv$ID) & grepl("DEL", sample_sv$ID)) {
    merged_size <- as.numeric(strsplit(as.character(merged_sv$ID),".",fixed = TRUE)[[1]][3])
    size_diff <- abs(merged_size - sample_sv$SIZE)
    return(any(pos_overlap & type_match & (size_diff <= 100)))
  } else {
    return(any(pos_overlap & type_match))
  }
}

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
      pos_col <- vcf[m,2]
      end_col <- strsplit(strsplit(as.character(vcf[m,8]),";",fixed = TRUE)[[1]][1],"=",fixed = T)[[1]][2]
    }
    else {
      id_col <- paste(vcf[m,2],vcf[m,3],sep = ".")
      pos_col <- vcf[m,2]
      end_col <- vcf[m,2]
    }

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

process_sample <- function(i) {
  print(i)
  sv_data <- process_mtb_vcf(tabledr[[i]])
  sample_id <- sample_names[i]
  if (grepl("del",sample_id)) {return()}
  sv_data$POS <- as.numeric(as.character(sv_data$POS))
  sv_data$END <- as.numeric(as.character(sv_data$END))
  sv_data$SIZE <- as.numeric(as.character(sv_data$SIZE))
  
  sample_sv_matrix <- numeric(nrow(merged_svs))  # Create a vector for the current sample's row
  
  for (j in seq_along(merged_svs$ID)) {
    merged_sv <- merged_svs[j, ]
    
    if (any(sv_data$ID == "NaN" & sv_data$POS <= merged_sv$END & sv_data$END >= merged_sv$POS)) {
      sample_sv_matrix[j] <- NA
    } else {
      overlap_result <- mapply(check_overlap, list(merged_sv), split(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS), seq(nrow(sv_data[sv_data$ID != "NaN",]%>%drop_na(POS)))))
      sample_sv_matrix[j] <- ifelse(any(overlap_result), 1, 0)
    }
  }
  
  return(matrix(sample_sv_matrix,nrow=1,ncol=nrow(merged_svs),dimnames=list(sample_id,merged_svs$ID)))  # Return the row for this sample
}

#cl <- makeCluster(100)
#clusterExport(cl, list("tabledr", "sample_names", "merged_svs", "sv_matrix", "process_mtb_vcf", "check_overlap"))
#clusterEvalQ(cl, library(dplyr))
#clusterEvalQ(cl, library(tidyverse))
#clusterEvalQ(cl, library(parallel))
sv_matrix_list = lapply(seq_along(tabledr), process_sample)

# Combine the list of sample rows into a matrix
sv_matrix <- do.call(rbind, sv_matrix_list)

write.table(sv_matrix,"/data/projects/punim1637/Aleix/SV/SV_analysis/sv_matrix.txt")
