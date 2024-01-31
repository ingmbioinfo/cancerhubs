library(openxlsx)
library(WriteXLS)

setwd('~/project_cancer/')

whole_datasets <- list()
formatted_datasets <- list()

# Multiple myeloma
Multiple_myeloma <- list()
df <- readRDS(file = './data/Multiple_myeloma_1')
whole_datasets[['Multiple_myeloma']] <- Multiple_myeloma
formatted_datasets[['Multiple_myeloma']] <- df[,c(1,7)]

# Breast cancer
Breast_cancer <- list()
df <- readRDS(file = './data/Breast_cancer_1')
Breast_cancer[['Breast_cancer_1']] <- df
df2 <- readRDS(file = './data/Breast_cancer_2')
Breast_cancer[['Breast_cancer_2']] <- df2
whole_datasets[['Breast_cancer']] <- Breast_cancer
df <- df[, c(7,15)]
df2 <- df2[,c(1,4)]
colnames(df)<-colnames(df2)
formatted_datasets[['Breast_cancer']] <- rbind(df, df2)

# CRC
Colon_cancer <- list()
df <- readRDS(file = './data/Colon_cancer_1')
Colon_cancer[['Colon_cancer_1']] <- df
whole_datasets[['Colon_cancer']] <- Colon_cancer
df <- df[, c(1,7)]
formatted_datasets[['Colon_cancer']] <- df

# Pancreatic cancer
Pancreatic_cancer <- list()
df <- readRDS(file = './data/Pancreatic_cancer_1')
Pancreatic_cancer[['Pancreatic_cancer_1']] <- df
df2 <- readRDS(file = './data/Pancreatic_cancer_2')
Pancreatic_cancer[['Pancreatic_cancer_2']] <- df2
df3 <- readRDS(file = './data/Pancreatic_cancer_3')
Pancreatic_cancer[['Pancreatic_cancer_3']] <- df3
whole_datasets[['Pancreatic_cancer']] <- Pancreatic_cancer
df <- df[, c(1,9)]
df2 <- df2[,c(2,9)]
df3 <- df3[,c(1,9)]
colnames(df)<-colnames(df2)
colnames(df3)<- colnames(df2)
formatted_datasets[['Pancreatic_cancer']] <- rbind(df, df2, df3)

# Prostate cancer
Prostate_cancer <- list()
df <- readRDS(file = './data/Prostate_cancer_1')
Prostate_cancer[['Prostate_cancer_1']] <- df
whole_datasets[['Prostate_cancer']] <- Prostate_cancer
df <- df[, c(1,48)]
formatted_datasets[['Prostate_cancer']] <- df


saveRDS(whole_datasets, './data/whole_datasets')

saveRDS(formatted_datasets, './data/formatted_datasets')


