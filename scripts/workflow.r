#Make sure the directories in the script match yours
setwd("..")

# Load custom functions from the 'functions.R' script
source('./scripts/functions.R')  # Assuming 'functions.R' is in the 'scripts' working directory

# Read the 'NON_ORF' data from an RDS file
NON_ORF <- readRDS('./data/NON_ORF')

# Print the current system time
Sys.time()

# Read various datasets from RDS files
whole_datasets <- readRDS(file = './data/whole_datasets')
formatted_datasets <- readRDS(file = './data/formatted_datasets')
interactors <- readRDS(file = './data/biogrid_interactors')

# Create an empty list to store results
all_results <- list()
isolation_list= list()

# Loop over the names of formatted datasets
for (i in names(formatted_datasets)) {
  
  # Select the tumor of interest
  tumor_of_interest <- i
  
  # Retrieve data and add mutation labels
  data <- retrieve_data(formatted_datasets, name = tumor_of_interest)
  mut_data <- add_mut_labels(data, NON_ORF)
  
  # Retrieve precog data and add to the mutation data
  precog <- retrieve_metaZ(tumor_of_interest)
  mut_data <- add_precog(df = mut_data, precog = precog)
  
  # Filter genes
  gene_list <- filter_genes(mut_data)
  
  # Count precog interactors and select those with count >= 10
  precog_interactors <- lapply(interactors[!mut_data$precog, ], function(x) sum(x %in% unique(mut_data[mut_data$precog, 1])))
  precog_interactors <- names(precog_interactors)[precog_interactors >= 10]
  
  # Extract cancer-specific interactors
  cancer_specific_interactors <- lapply(interactors$as.matrix.interactors., function(x) x[x %in% precog_interactors])
  cancer_specific_interactors <- data.frame(as.matrix(cancer_specific_interactors))
  
  # Create a list to store results for the current tumor
  res_list <- list()
  
  # Compute network results
  res_list[['All_Genes']] <- result(mut_data, precog, interactors = cancer_specific_interactors, gene_list = gene_list)
  
  # Filter results based on precog_type (excluding 'none')
  res_list[['PRECOG']] <- filters(res_list[['All_Genes']], columns = c('precog_type'), filters = c('none'), filter_out = TRUE)
  res_list[['Non_PRECOG']] <- filters(res_list[['All_Genes']], columns = c('precog_type'), filters = c('none'), filter_out = FALSE)

  #Filter results by genes significant in precog and not mutated                                 
  res_list[["Only_PRECOG"]] =  filters(res_list[['All_Genes']], columns = c('mutation'), filters = c("NONE"), filter_out= FALSE) 

  #removing precog_type column                                      
  remove_precog_type <- function(df) { df[, !colnames(df) %in% "precog_type"] } 
  res_list <- lapply(res_list, remove_precog_type)    
                                        
  # Compute isolation results
  isolation_result <- result(mut_data, precog, interactors = cancer_specific_interactors, gene_list = gene_list, isolation_score = TRUE) 
  isolation_result <- remove_precog_type(isolation_result) 
  isolation_list[[tumor_of_interest]] <- isolation_result
  
  # Save results for the current tumor in the main list
  all_results[[tumor_of_interest]] <- res_list
}

# Save the final results list to an RDS file
saveRDS(all_results, './result/all_results.rds')
saveRDS(isolation_list, './result/isolation_list.rds')                                        

# Print the current system time again
Sys.time()
