setwd("..")
# Load custom functions from the 'functions.R' script
source('./scripts/functions.R')  # Assuming 'functions.R' is in the 'project_cancer' working directory

# Set the working directory to 'project_cancer'
# (Note: It's mentioned as a comment; make sure it's appropriately set before running the script)

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
  res_list[['all_result_network']] <- result(mut_data, precog, interactors = cancer_specific_interactors, gene_list = gene_list)
  
  # Filter results based on precog_type (excluding 'none')
  res_list[['precog_result_network']] <- filters(res_list[['all_result_network']], columns = c('precog_type'), filters = c('none'), filter_out = TRUE)
  res_list[['non_precog_result_network']] <- filters(res_list[['all_result_network']], columns = c('precog_type'), filters = c('none'), filter_out = FALSE)
  
  # Compute isolation results
  res_list[['all_result_isolation']] <- result(mut_data, precog, interactors = cancer_specific_interactors, gene_list = gene_list, isolation_score = TRUE)
  
  # Save results for the current tumor in the main list
  all_results[[tumor_of_interest]] <- res_list
}

# Save the final results list to an RDS file
saveRDS(all_results, './result/all_results')

# Print the current system time again
Sys.time()
