<img src="cancerhubs_logo.png" align="right" alt="" width="200" />

# CancerHubs: Systematic Identification of Cancer-Related Protein Interaction Hubs

🧪 **Explore interactively**: [https://cancerhubs.app](https://cancerhubs.app)  
📄 **Published article**: [Briefings in Bioinformatics (2025)](https://doi.org/10.1093/bib/bbae635)  
📦 **Shiny App Source**: [cancerhubs_shiny](https://github.com/ingmbioinfo/cancerhubs_shiny)  
📘 **Reproducibility Package**: [cancerhubs_paper](https://github.com/ingmbioinfo/cancerhubs_paper)

---

## Overview
CancerHubs is a novel computational framework designed to predict proteins and pathways involved in cancer by integrating mutational data, clinical outcome predictions, and interactomics. This method ranks genes based on the number of mutated interactors their corresponding proteins have, defining hubs of mutated proteins with potential relevance for cancer research and therapy.

---

## Data Sources and Methodology

- **Mutational Data**: Sourced from previously published datasets for selected cancers (Multiple Myeloma, Breast Cancer, Prostate Cancer, Colorectal Cancer, and Pancreatic Cancer).
- **Clinical Outcome Predictions**: Utilizes Precog Meta-Z data to correlate gene expression with overall survival.
- **Interactomics**: Interaction data are extracted from the BioGrid database, focusing on interactions between proteins encoded by genes found in the mutated gene lists.

---

## Pipeline Description

1. **Data Retrieval**: Mutation datasets are collected, and a list of mutated genes is generated for each cancer type.
2. **Clinical Outcome Correlation**: Z-scores quantifying the correlation between gene expression and patient survival are merged with mutational data.
3. **Gene Classification**: Genes are classified based on their mutation status and correlation with clinical outcomes into MUT (exclusively mutated), PRECOG (expression correlated with outcomes but not mutated), and MUT + PRECOG (both mutated and expression correlated with outcomes).
4. **Interactome Determination**: The global interactome of each cancer-related gene is defined.
5. **Network Score Calculation**: A network score is calculated for each gene based on the number of mutated interactors, ranking genes to define protein hubs.

---

## Repository Structure

This repository contains the **core pipeline** used to process mutation data, integrate clinical scores, and calculate network-based gene rankings.

### Key Scripts

- `biogrid_extractor.r` — extracts gene interaction information from BioGRID.
- `workflow.r` — main analysis pipeline.
- `functions.r` — functions for computing the network score.

### Output

- `all_results.rds`: Core R object containing rankings by cancer and category (used in the app and downstream analyses).

---

## Explore Our Results

An interactive Shiny web application is now available at:  
👉 [https://cancerhubs.app](https://cancerhubs.app)

Use it to:
- Search gene rankings across cancers
- Visualise pan-cancer hub scores
- Browse interactors and 3D networks
- Download gene-specific data and networks

The app source code is available at:  
🔗 [https://github.com/ingmbioinfo/cancerhubs_shiny](https://github.com/ingmbioinfo/cancerhubs_shiny)

---

## Data Availability

All underlying data and scripts are available in this GitHub repository. The datasets used in this study are derived from public sources, and the specific versions are documented within the repository.

To **reproduce the results of the paper** and interactively explore them, visit the reproducibility repository:  
📘 [https://github.com/ingmbioinfo/cancerhubs_paper](https://github.com/ingmbioinfo/cancerhubs_paper)

---

## Funding

This research was funded by **Associazione Italiana per la Ricerca sul Cancro (AIRC)**, under **MFAG 2021 ID 26178** project to **Nicola Manfrini**.

---

## Contact

For questions or support, please contact:  
📧 manfrini@ingm.org  
📧 ferrari@ingm.org

---

## Citation

If you use CancerHubs in your research, please cite:

> Ivan Ferrari, Federica De Grossi, Giancarlo Lai, Stefania Oliveto, Giorgia Deroma, Stefano Biffo, Nicola Manfrini.  
> **CancerHubs: a systematic data mining and elaboration approach for identifying novel cancer-related protein interaction hubs**.  
> _Briefings in Bioinformatics_, Volume 26, Issue 1, January 2025.  
> [https://doi.org/10.1093/bib/bbae635](https://academic.oup.com/bib/article/26/1/bbae635/7918695?searchresult=1)
