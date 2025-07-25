# SEA-AD Bioinformatic Dashboard

An interactive Streamlit-based dashboard for exploring multiomic single-nucleus data from the Seattle Alzheimer's Disease Cohort (SEA-AD).  
This tool enables deep analysis of matched scRNA-seq and scATAC-seq data, including progression modeling, integrative differential analysis, and enrichment testing across cell types and disease states.

## Key Features
- **Data Integration**
  - Joint analysis of matched scRNA-seq and scATAC-seq at the gene level
  - Interactive filtering and grouping by:
    - Cell subtype  
    - Donor/sample metadata  
    - Quality metrics  
    - Disease progression scores  
  - Multiomic support for:
    - Differential expression & accessibility  
    - Gene set enrichment  
    - Progression/survival modeling  
    - Comparative visualization across modalities  
- **Analysis and Visualization**
  - **Dimensionality reduction**
    - Colored by metadata, RNA expression, or ATAC accessibility  
    - Joint UMAP comparisons (RNA vs ATAC modalities)  
  - **Feature comparison**
    - Compare gene expression/accessibility across metadata-defined groups  
    - Support for pseudobulk estimates with error bars  
  - **Survival/Progression modeling**
    - Correlation between features and continuous disease progression scores  
    - Visual overlays on embeddings  
  - **Differential analysis**
    - Run differential expression for RNA and differential accessibility for ATAC simultaneously  
    - LogFC comparison plots across modalities  
  - **Enrichment analysis**
    - Run ranked GSEA for RNA and ATAC differential results  
    - Normalized enrichment score comparative plots for both modalities  

## Project Structure (Updated)

```bash
SEA_AD_bioinformatic_dashboard/
├── data/                         # Input and intermediate data files (not tracked)
├── RNA_ATAC_workflow.ipynb       # Main workflow notebook
├── streamlit_dashboard_forGit.py # Streamlit dashboard interface
├── results/                      # Output files (aggregated, filtered, imputed) (not tracked)
├── scripts/                      # Legacy scripts (some steps replaced by notebooks) (not tracked)
├── README.md
```
Updated Workflow
The current pipeline is implemented in the notebook:
RNA_ATAC_workflow.ipynb

  1. Load individual donor H5AD files for RNA and ATAC
  2. Subset cells by subtype and donor
  3. Normalize and log-transform expression and accessibility data
  4. Map ATAC peaks to genes and align with RNA features
  5. Merge RNA and ATAC by common cells and genes
  7. Save processed subsets for dashboard exploration

Interactive Dashboard
The streamlit_dashboard_forGit.py file implements an interactive Streamlit dashboard that allows users to:

  1. Load and filter cell subsets by metadata (cell type, donor, etc.)
  2. Define groups for comparative analysis
  3. View dimensionality reduction plots (UMAP, PCA)
  4. Perform feature comparisons, differential results, and enrichment
  5. Visualize progression-related molecular changes

Getting Started
This section will include:

  1. Installation instructions
  2. Environment setup (conda, pip, requirements.txt)
  3. Example commands to run the notebook and dashboard

Notes

References for SEA-AD

Seattle Alzheimer’s disease - brain-map.org. https://portal.brain-map.org/explore/seattle-alzheimers-disease

Gabitto MI et al. Integrated multimodal cell atlas of Alzheimer’s disease. bioRxiv (2023) doi: 10.1101/2023.05.08.539485.

