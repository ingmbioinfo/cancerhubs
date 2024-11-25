<img src="cancerhubs_logo.png" align="right" alt="" width="200" />



# CancerHubs: Systematic Identification of Cancer-Related Protein Interaction Hubs

## Overview
CancerHubs is a novel computational framework designed to predict proteins and pathways involved in cancer by integrating mutational data, clinical outcome predictions, and interactomics. This method ranks genes based on the number of mutated interactors their corresponding proteins have, defining hubs of mutated proteins with potential relevance for cancer research and therapy.

## Data Sources and Methodology
- **Mutational Data**: Sourced from previously published datasets for selected cancers (Multiple Myeloma, Breast Cancer, Prostate Cancer, Colorectal Cancer, and Pancreatic Cancer).
- **Clinical Outcome Predictions**: Utilizes Precog Meta-Z data to correlate gene expression with overall survival.
- **Interactomics**: Interaction data are extracted from the BioGrid database, focusing on interactions between proteins encoded by genes found in the mutated gene lists.

## Pipeline Description
1. **Data Retrieval**: Mutation datasets are collected, and a list of mutated genes is generated for each cancer type.
2. **Clinical Outcome Correlation**: Z-scores quantifying the correlation between gene expression and patient survival are merged with mutational data.
3. **Gene Classification**: Genes are classified based on their mutation status and correlation with clinical outcomes into MUT (exclusively mutated), PRECOG (expression correlated with outcomes but not mutated), and MUT + PRECOG (both mutated and expression correlated with outcomes).
4. **Interactome Determination**: The global interactome of each cancer-related gene is defined.
5. **Network Score Calculation**: A network score is calculated for each gene based on the number of mutated interactors, ranking genes to define protein hubs.

## Repository Content
- **Scripts**: Contains R scripts used for data management, processing, and network score calculation.
  - `biogrid_extractor.r`: Extracts gene interaction information from BioGRID data.
  - `workflow.r`: Directs the analytical process, incorporating data handling and integration.
  - `functions.r`: Provides functions required for computing the network score.
- **Data**: Mutation datasets, Precog MetaZ-scores, and interactomic data from BioGRID used in the study.
- **Results**: Output file in RDS format containing the genes ranking based on the network scores for each tumor type analyzed.

## Usage
Refer to the `workflow.r` script to understand the sequence of operations performed from data preparation to network score calculation. The scripts require R programming environment and access to the datasets provided in the Data folder.

## Work in progress
This repository will be soon updated with data of other tumour types and a shiny app will be available to interactively explore our results.

## Data Availability
All underlying data and scripts are available on this GitHub repository. The datasets used in this study are derived from public sources, and the specific versions are documented within the repository.

## Supplementary Data
Supplementary data, including detailed tables of network scores and gene classifications, are provided for further exploration of cancer-specific and broad-cancer protein hubs.

## Funding
This research was funded by Associazione Italiana per la Ricerca sul Cancro (AIRC), under MFAG 2021 ID 26178 project to Nicola Manfrini.

## Contact
For questions or support, please contact manfrini@ingm.org and ferrari@ingm.org. 

## Citation
If you use CancerHubs in your research, please cite our paper:

Ivan Ferrari, Giancarlo Lai, Federica De Grossi, Stefania Oliveto, Stefano Biffo, Nicola Manfrini. "CancerHubs: A Systematic Data Mining and Elaboration Approach for Identifying Novel Cancer-Related Protein Interaction Hubs." [Journal, Volume, Year].
