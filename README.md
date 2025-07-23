SEA-AD Bioinformatic Dashboard
Overview
This project supports the analysis and visualization of single-nucleus RNA-seq and ATAC-seq data from the Seattle Alzheimer's Disease (SEA-AD) cohort. It includes workflows for processing, integrating, and visualizing multiomic data across donors and cell types, with the aim of identifying molecular signatures associated with disease progression.

Key components include:

Preprocessing & integration of RNA-seq and ATAC-seq data at the gene level

Multi-donor aggregation and filtering

Imputation and normalization of sparse matrices

Correlation and enrichment analysis

Dashboard creation for interactive exploratory analysis

Project Structure (Updated)
```bash
SEA_AD_bioinformatic_dashboard/
├── data/                         # Input and intermediate data files (not tracked)
├── streamlit_dashboard_forGit.py    # Main workflow notebook
├── streamlit_dashboard_forGit.py                    # Streamlit dashboard interface
├── results/                      # Output files (aggregated, filtered, imputed) (not tracked)
├── scripts/                      # Legacy scripts (some steps replaced by notebooks) (not tracked)
├── README.md
```
Updated Workflow
The current pipeline is implemented in the notebook:

streamlit_dashboard_forGit.py

Load individual donor H5AD files for RNA and ATAC

Subset cells by subtype and donor

Normalize and log-transform expression and accessibility data

Map ATAC peaks to genes and align with RNA features

Merge RNA and ATAC by common cells and genes

Perform imputation (e.g., SoftImpute, or deep learning models)

Save processed subsets for dashboard exploration

Interactive Dashboard
The streamlit_dashboard_forGit.py file implements an interactive Streamlit dashboard that allows users to:

Load and filter cell subsets by metadata (cell type, donor, etc.)

Define groups for comparative analysis

View dimensionality reduction plots (UMAP, PCA)

Perform gene-level comparisons and enrichment

Visualize progression-related molecular changes

Getting Started
This section will include:

Installation instructions

Environment setup (conda, pip, requirements.txt)

Example commands to run the notebook and dashboard

Notes
This section can list any assumptions, data access notes, or links to SEA-AD metadata.

