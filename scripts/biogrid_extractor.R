# Read the Biogrid data from a tab-separated file into a data frame
biogrid <- read.csv('data/BIOGRID-ORGANISM-Homo_sapiens-3.5.181.tab2.txt', sep = '\t')

# Create an empty list to store interactors
interactors <- list()

# Iterate over the union of Official.Symbol.Interactor.A and Official.Symbol.Interactor.B
for (i in union(biogrid$Official.Symbol.Interactor.A, biogrid$Official.Symbol.Interactor.B)) {
  
    # Extract interactors for each gene (i)
    interactors[[i]] <- union(
        biogrid$Official.Symbol.Interactor.B[biogrid$Official.Symbol.Interactor.A == i],
        biogrid$Official.Symbol.Interactor.A[biogrid$Official.Symbol.Interactor.B == i]
    )
}

# Convert the list of interactors to a matrix
interactors <- as.matrix(interactors)

# Save the matrix of interactors to an RDS file
saveRDS(object = interactors, file = '/mnt/home/ferrari/project_cancer/data/biogrid_interactors')
