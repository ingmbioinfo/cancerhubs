<img src="cancerhubs_logo.png" align="right" alt="" width="200" />

# CancerHubs: Systematic Identification of Cancer-Related Protein Interaction Hubs

🧪 **Explore interactively**: [https://cancerhubs.app](https://cancerhubs.app)  
📄 **Published article**: [Briefings in Bioinformatics (2025)](https://doi.org/10.1093/bib/bbae635)  
📦 **Shiny App Source**: [cancerhubs_shiny](https://github.com/ingmbioinfo/cancerhubs_shiny)  
📘 **Reproducibility Package**: [cancerhubs_paper](https://github.com/ingmbioinfo/cancerhubs_paper)

---

## Overview

CancerHubs is a computational method that identifies **cancer-relevant protein hubs** by integrating:
- Somatic mutation data
- Prognostic gene expression scores (PRECOG)
- Protein–protein interaction networks (BioGRID)

By ranking genes based on the **number of mutated interactors**, it prioritises central players in cancer-related pathways—beyond mutation frequency alone. This network-centric perspective allows the discovery of potential driver genes otherwise overlooked by traditional methods.

### Network Score

Each gene is assigned a **Network Score** defined as:

```
Network Score = (# Mutated Interactors)^2 / (# Total Interactors)
```

This formula reflects how strongly a gene is embedded in a network of cancer-mutated interactors. Higher scores indicate greater centrality and functional relevance in tumour contexts.

> **Note:** The Network Score is computed exclusively from somatic mutation data. It does **not** take into account copy-number variations or structural alterations.

---

## Data Sources and Methodology

- **Mutational Data**: Sourced from previously published datasets for selected cancers (Multiple Myeloma, Breast Cancer, Prostate Cancer, Colorectal Cancer, and Pancreatic Cancer).
- **Clinical Outcome Predictions**: Utilizes Precog Meta-Z data to correlate gene expression with overall survival.
- **Interactomics**: Interaction data are extracted from the BioGrid database, focusing on interactions between proteins encoded by genes found in the mutated gene lists.

---

## Pipeline Overview

- **1. Data Retrieval**  
  Mutation data are collected for each tumour type. Genes are filtered for coding/non-coding status.

- **2. Clinical Outcome Correlation**  
  Gene expression–survival associations are obtained from PRECOG and integrated as meta-Z scores.

- **3. Gene Classification**  
  Genes are labelled as:  
  - **MUT** – mutated only  
  - **PRECOG** – prognosis-associated only  
  - **MUT + PRECOG** – both

- **4. Network Construction**  
  BioGRID interactions define the gene interactome.

- **5. Network Scoring**  
  Genes are scored and ranked by the fraction and number of mutated interactors.

---

## Repository Structure

This repository contains the **core pipeline** used to process mutation data, integrate clinical scores, and calculate network-based gene rankings.

### Contents

- **/scripts**  
  Core R scripts for mutation integration and scoring.

- **/data**  
  Input datasets (formatted mutation and interactome data).

- **/result/all_results.rds**  
  Main output object containing ranked genes, used by the Shiny app.

- **/docs**  
  Supplementary data and mutation annotation summaries.

---

## 🌐 Web Application

An interactive Shiny web application is now available:  
👉 [https://cancerhubs.app](https://cancerhubs.app)

Use it to:
- Search gene rankings across cancers
- Visualise pan-cancer hub scores
- Browse interactors and 3D networks
- Download gene-specific data and networks

App source code:  
🔗 [https://github.com/ingmbioinfo/cancerhubs_shiny](https://github.com/ingmbioinfo/cancerhubs_shiny)

---

## 🔄 How to Reproduce

To run the full pipeline locally:

1. Install required packages (see headers in `functions.r`)
2. Clone this repository and set the working directory
3. Execute the pipeline:
   ```r
   source("workflow.r")
   ```

For full reproducibility with frozen data versions, refer to the companion repository:  
📘 [`cancerhubs_paper`](https://github.com/ingmbioinfo/cancerhubs_paper)

---

## Data Availability

All underlying data and scripts are available in this GitHub repository. The datasets used in this study are derived from public sources, and the specific versions are documented within the repository.

### Included Files

- [`all_results.rds`](https://github.com/ingmbioinfo/cancerhubs/blob/main/result/all_results.rds): Gene-wise rankings for each tumour type, including Network Scores and annotation categories.
- [`genes_interactors_list.rds`](https://github.com/ingmbioinfo/cancerhubs/blob/main/result/genes_interactors_list.rds): Curated interactome data mapping each gene to its interaction partners.
- [`formatted_datasets`](https://github.com/ingmbioinfo/cancerhubs/tree/main/data/formatted_datasets): Input gene tables formatted for pipeline ingestion.
- [`biogrid_interactors`](https://github.com/ingmbioinfo/cancerhubs/blob/main/data/biogrid_interactors): Full PPI interaction records derived from BioGRID.
- [`Mutational Data Summary`](https://github.com/ingmbioinfo/cancerhubs/blob/main/Mutational%20Data.pdf): PDF with curation criteria and source references for mutation datasets.

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
