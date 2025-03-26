# Function to retrieve a dataset of interest
retrieve_data <- function(formatted_datasets, name){
    return(formatted_datasets[[name]])
}

# Function to add mutation site and gene mutation type columns
add_mut_labels <- function(data, NON_ORF){
    data$mut_site <- 1:dim(data)[1]
    data$mut_site[which(data[,2] %in% NON_ORF)] <- 'NON_ORF'
    data$mut_site[which(!(data[,2] %in% NON_ORF))]  <- 'ORF'
    
    ONO <- data.frame(table(data[,1], data[,3]))
    data$gene_mut_type <- 1:dim(data)[1]
    
    non_orf_mutated_gene <- ONO[ONO$Freq>0 & ONO$Var2=='NON_ORF' ,][,1]
    orf_mutated_gene <- ONO[ONO$Freq>0 & ONO$Var2=='ORF' ,][,1]
    both <- intersect(non_orf_mutated_gene, orf_mutated_gene)
    
    data$gene_mut_type[which(data[,1] %in% non_orf_mutated_gene)] <- 'NON_ORF'
    data$gene_mut_type[which(data[,1] %in% orf_mutated_gene)] <- 'ORF'
    data$gene_mut_type[which(data[,1] %in% both)] <- 'BOTH'
    
    return(data)
}

# Function to retrieve metaZ from precog
retrieve_metaZ <- function(name){
    precog <- read.csv('./data//PRECOG-metaZ.pcl', sep='\t')
    precog <- precog[, c(1, 2, which(colnames(precog) == name))]
    precog$type <- 1:dim(precog)[1]
    precog[which(precog[,3] >= 1.96), 4] <- 'oncogene' 
    precog[which(precog[,3] <= (-1.96)), 4] <- 'tumor suppressor' 
    precog[which(abs(precog[,3]) < (1.96)), 4] <- 'none' 
    return(precog)
}

# Function to add precog column to mut_data
add_precog <- function(df, precog){
    df$precog <- df[, 1] %in% precog$Gene[precog$type != 'none'] 
    return(df)
}

# Function to filter genes of interest
filter_genes <- function(mut_data){
    md <- median(table(mut_data[, 1]))
    min <- min(table(mut_data[, 1]))
    
    orf_precog <- data.frame(table(mut_data[, 1][(mut_data$gene_mut_type == 'ORF' & mut_data$precog)]))
    orf_precog <- as.vector(orf_precog[orf_precog[, 2] >= (md - min), 1])
    
    orf_precog_neg <- data.frame(table(mut_data[, 1][(mut_data$gene_mut_type == 'ORF' & !mut_data$precog)]))
    orf_precog_neg <- as.vector(orf_precog_neg[orf_precog_neg[, 2] >= (md - min), 1])
    
    if(length(unique(mut_data$gene_mut_type)) > 1){
        non_orf_precog <- data.frame(table(mut_data[, 1][(mut_data$gene_mut_type == 'NON_ORF' & mut_data$precog)]))
        non_orf_precog <- as.vector(non_orf_precog[non_orf_precog[, 2] >= (md + min), 1])
        
        non_orf_precog_neg <- data.frame(table(mut_data[, 1][(mut_data$gene_mut_type == 'NON_ORF' & !mut_data$precog)]))
        non_orf_precog_neg <- as.vector(non_orf_precog_neg[non_orf_precog_neg[, 2] >= (3 * md + min), 1])
        
        both_precog <- data.frame(table(mut_data[, 1][(mut_data$gene_mut_type == 'BOTH' & mut_data$precog)]))
        both_precog <- as.vector(both_precog[both_precog[, 2] >= (md - min), 1])
        
        both_precog_neg <- mut_data[(mut_data$gene_mut_type == 'BOTH' & !mut_data$precog), c(1, 3)]
        both_precog_neg <- as.vector(rownames(table(both_precog_neg))[(table(both_precog_neg)[, 1] >= (3 * md) | table(both_precog_neg)[, 2] >= (md - min))])
        
        gene_list <- c(orf_precog, orf_precog_neg, non_orf_precog, non_orf_precog_neg, both_precog, both_precog_neg)
    } else {
        gene_list <- c(orf_precog, orf_precog_neg)
    }
    
    only_precog <- setdiff(precog$Gene, gene_list)
    precog[(abs(precog[, 3]) < (2.58)) & (precog$Gene %in% only_precog),  4] <- 'none'  
    sign_only_precog <- precog$Gene[(abs(precog[, 3]) >= (2.58)) & (precog$Gene %in% only_precog)]
    gene_list <- c(gene_list, sign_only_precog)
    gene_list <- gene_list[gene_list %in% precog$Gene]
    
    return(gene_list)
}

# Function to obtain the final dataframe
result <- function(mut_data, precog, interactors, gene_list, isolation_score = FALSE) {
  result <- data.frame(gene_list)
  result$tot_int <- lapply(interactors[result$gene_list, ], length)
  itrs <- lapply(interactors[result$gene_list, ], function(x) x %in% unique(mut_data[, 1]))
  result$mut_int <- lapply(itrs, sum)
  result$tot_int <- as.double(result$tot_int)
  result$mut_int <- as.double(result$mut_int)
  result$mut_perc <- result$mut_int / result$tot_int
  result$mut_perc[is.na(result$mut_perc)] <- 0
  result$network_score <- result$mut_perc * result$mut_int

  if (length(unique(mut_data[, 4])) > 1) {
    both <- unique(mut_data[, 1])[(table(mut_data[, 4], mut_data[, 1]) > 0)[1, ]]
    non_orf <- unique(mut_data[, 1])[(table(mut_data[, 4], mut_data[, 1]) > 0)[2, ]]
    orf <- unique(mut_data[, 1])[(table(mut_data[, 4], mut_data[, 1]) > 0)[3, ]]
    result$mutation <- 1:dim(result)[1]
    result$mutation[which(result[, 1] %in% non_orf)] <- 'NON_ORF'
    result$mutation[which(result[, 1] %in% orf)] <- 'ORF'
    result$mutation[which(result[, 1] %in% both)] <- 'BOTH'
    result$mutation[(result$mutation != 'NON_ORF')  & (result$mutation != 'ORF') & (result$mutation != 'BOTH')] <- 'NONE'
  } else {
    orf <- unique(mut_data[, 1])[(table(mut_data[, 4], mut_data[, 1]) > 0)[1, ]]
    result$mutation <- 1:dim(result)[1]
    result$mutation[which(result[, 1] %in% orf)] <- 'ORF'
    result$mutation[(result$mutation != 'ORF')] <- 'NONE'
  }

  prc <- data.frame(precog$type)
  rownames(prc) <- precog$Gene
  prc$precog_metaZ <- precog[, 3]
  result$precog_type <- prc[result[, 1], 1]
  result$precog_metaZ <- prc[result[, 1], 2]
  result$isolation_score <- abs(result$precog_metaZ) / (result$network_score + 1)

  if (isolation_score == TRUE) {
    result <- result[order(-result$isolation_score), ]
    rownames(result) <- 1:dim(result)[1]

    return(result[, c(1, 6, 7, 8, 2, 3, 9)])
  } else {
    result <- result[order(-result$network_score), ]
    rownames(result) <- 1:dim(result)[1]

    return(result[, c(1, 6, 7, 8, 2, 3, 5)])
  }
}

# Function to filter results
filters <- function(res, columns, filters, filter_out = FALSE) {
  if (filter_out == FALSE) {
    for (i in 1:length(filters)) {
      res <- res[res[[columns[i]]] == filters[i], ]
    }
  } else {
    for (i in 1:length(filters)) {
      res <- res[res[[columns[i]]] != filters[i], ]
    }
  }

  rownames(res) <- 1:dim(res)[1]
  return(res)
}
                 
                 
#Obtain interactors for each gene               
get_gene_interactors <- function(gene_list, cancer_specific_interactors, mut_data, res_list, precog) {
  gene_interactors = list()
  inter <- data.frame(gene_list)
  
  # Store the list of interactors instead of just the count
  inter$total_interactors <- lapply(gene_list, function(x) cancer_specific_interactors$as.matrix.cancer_specific_interactors.[[x]]) 
  inter$mutated_interactors <- lapply(inter$total_interactors, function(x) x[x %in% unique(mut_data[, 1])])                                                                 
  gene_interactors[["inter"]] <- inter
  
  genes = res_list$PRECOG$gene_list                                      
  precog_inter = data.frame(genes)
  
  precog_inter$precog = lapply(res_list$PRECOG$gene_list, function(x) cancer_specific_interactors$as.matrix.cancer_specific_interactors.[[x]])    
  precog_inter$precog = lapply(precog_inter$precog, function(x) x[x %in% unique(precog[precog$type != "none", 1])])                           
  
  precog_inter$precog_mut = lapply(precog_inter$precog, function(x) x[x %in% unique(mut_data[mut_data$precog, 1])]) 
  
  precog_inter = precog_inter[order(precog_inter$genes),]                                 
  gene_interactors[["precog_inter"]] = precog_inter  
  
  return(gene_interactors)
} 
                                   
                                   
                                   
library(openxlsx)

generate_all_interactomes <- function(data, original, output_dir) {
  
  scenarios <- list(
    list(use_precog = TRUE, include_mutated = FALSE, suffix = "_precog"),
    list(use_precog = TRUE, include_mutated = TRUE, suffix = "_precog_mut"),
    list(use_precog = FALSE, include_mutated = FALSE, suffix = "_all_genes"),
    list(use_precog = FALSE, include_mutated = TRUE, suffix = "_all_genes_mut")
  )
  
  for (cancer_type in names(data)) {
    
    cancer_data <- data[[cancer_type]]
    if (is.null(cancer_data)) next  # Skip invalid cancer types
    
    for (scenario in scenarios) {
      
      use_precog <- scenario$use_precog
      include_mutated <- scenario$include_mutated
      suffix <- scenario$suffix
      
      sel_data <- if (use_precog) cancer_data[["precog_inter"]] else cancer_data[["inter"]]
      selection <- if (use_precog) {
        if (include_mutated) "precog_mut" else "precog"
      } else {
        if (include_mutated) "mutated_interactors" else "total_interactors"
      }
      
      if (!(selection %in% names(sel_data))) next  # Skip if selection column is missing
      
      original_type <- if (use_precog) "PRECOG" else "All_Genes"
      column_name <- if ("gene_list" %in% colnames(sel_data)) "gene_list" else "genes"
      
      nodes <- original[[cancer_type]][[original_type]]
      nodes <- nodes[nodes$gene_list %in% unique(sel_data[[column_name]]), ]
        
      # Extract all genes (original + interactors)
      all_genes <- unique(c(nodes$gene_list, unlist(sel_data[[selection]])))

      # Find missing interactors
      missing_genes <- setdiff(all_genes, nodes$gene_list)

      if (length(missing_genes) > 0) {
        # Create an empty dataframe with the same columns as `nodes`
        missing_nodes <- as.data.frame(matrix(NA, nrow = length(missing_genes), ncol = ncol(nodes)))
        colnames(missing_nodes) <- colnames(nodes)  # Ensure column names match
        missing_nodes$gene_list <- missing_genes    # Fill only gene_list column
  
       # Append new interactors to nodes
        nodes <- rbind(nodes, missing_nodes)
         }  
      
      edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
      
      for (i in seq_along(nodes$gene_list)) {
        gene <- nodes$gene_list[i]
        interactors_gene <- sel_data[sel_data[[column_name]] == gene, selection]
        
        if (length(interactors_gene) > 0) {
          interactors <- unlist(interactors_gene)
          new_edges <- data.frame(from = pmin(gene, interactors), to = pmax(gene, interactors), stringsAsFactors = FALSE)
          edges <- rbind(edges, new_edges)
        }
      }
      
      # Write to Excel
      file_name <- file.path(output_dir, paste0(cancer_type, suffix, ".xlsx"))
      wb <- createWorkbook()
      addWorksheet(wb, "Nodes")
      addWorksheet(wb, "Edges")
      writeData(wb, "Nodes", nodes)
      writeData(wb, "Edges", edges)
      saveWorkbook(wb, file_name, overwrite = TRUE)
      
      if (!file.exists(file_name)) {
        print(paste("Failed to save file:", file_name))
        }
    }
  }
}
                                   
                                   
                                   
