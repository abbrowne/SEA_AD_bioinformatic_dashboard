SEA-AD Bioinformatic Dashboard
Overview
This project supports the analysis and visualization of single-nucleus RNA-seq and ATAC-seq data from the Seattle Alzheimer's Disease (SEA-AD) cohort. It includes workflows for processing, integrating, and visualizing multiomic data across donors and cell types, with the aim of identifying molecular signatures associated with disease progression.

Key components include:
  1. Preprocessing & integration of RNA-seq and ATAC-seq data at the gene level
  2. Multi-donor aggregation and filtering
  3. Imputation and normalization of sparse matrices
  4. Correlation and enrichment analysis
  5. Dashboard creation for interactive exploratory analysis

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

  1. Load individual donor H5AD files for RNA and ATAC
  2. Subset cells by subtype and donor
  3. Normalize and log-transform expression and accessibility data
  4. Map ATAC peaks to genes and align with RNA features
  5. Merge RNA and ATAC by common cells and genes
  6. Perform imputation (e.g., SoftImpute, or deep learning models)
  7. Save processed subsets for dashboard exploration

Interactive Dashboard
The streamlit_dashboard_forGit.py file implements an interactive Streamlit dashboard that allows users to:

  1. Load and filter cell subsets by metadata (cell type, donor, etc.)
  2. Define groups for comparative analysis
  3. View dimensionality reduction plots (UMAP, PCA)
  4. Perform gene-level comparisons and enrichment
  5. Visualize progression-related molecular changes

Getting Started
This section will include:

  1. Installation instructions
  2. Environment setup (conda, pip, requirements.txt)
  3. Example commands to run the notebook and dashboard

Notes

SEA AD Study link: https://portal.brain-map.org/explore/seattle-alzheimers-disease/whatis

