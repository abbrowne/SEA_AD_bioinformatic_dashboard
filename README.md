# Cell-type-Specific Early Vulnerability Signatures in Alzheimer's Disease

## Project Overview
This project analyzes multimodal single-cell data from the SEA-AD dataset to identify early vulnerability signatures in specific cell types during Alzheimer's disease progression, with a focus on layer 2/3 IT excitatory neurons and Sst+ interneurons.

## Project Structure
```
.
├── data/                      # Raw and processed data
│   ├── raw/                  # Raw data files
│   ├── processed/            # Processed data files
│   └── metadata/             # Sample metadata and annotations
├── notebooks/                # Jupyter notebooks for analysis
├── scripts/                  # Python/R scripts
│   ├── preprocessing/        # Data preprocessing scripts
│   ├── analysis/            # Analysis scripts
│   └── visualization/       # Visualization scripts
├── results/                  # Analysis results
│   ├── figures/             # Generated figures
│   └── tables/              # Output tables
└── requirements.txt          # Python dependencies
```

## Setup Instructions

1. Create a conda environment:
```bash
conda create -n sea_ad python=3.9
conda activate sea_ad
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Multiomic Data Aggregation

For multiomic snRNAseq+snATACseq analysis, the sequencing data is split by type and donor. The `aggregate_multiomic.py` script aggregates data from multiple donors into separate RNA and ATAC AnnData files for independent analysis.

### Usage

```bash
python scripts/analysis/aggregate_multiomic.py \
    --metadata-file data/metadata/multiomic_metadata.csv \
    --rna-data-dir data/rna_files \
    --atac-data-dir data/atac_files \
    --output-dir results/aggregated_multiomic
```

### Expected Data Structure

```
data/
├── metadata/
│   └── multiomic_metadata.csv  # Metadata file with multiomic sample info
├── rna_files/
│   ├── donor_001_rna_data.h5ad
│   ├── donor_002_rna_data.h5ad
│   └── ...
└── atac_files/
    ├── donor_001_atac_data.h5ad
    ├── donor_002_atac_data.h5ad
    └── ...
```

### Metadata Requirements

The metadata CSV file should contain:
- `donor_id`: Unique identifier for each donor
- `sample_id`: Unique identifier for each sample/cell
- `cell_type`: Cell type annotation (optional)
- `cps_score`: Cognitive Performance Score (optional)
- `braak_stage`: Braak stage (optional)

**Note**: All samples included in the metadata file are assumed to have multiomic data.

### Output Files

The script generates:
- `aggregated_rna_data.h5ad`: Aggregated RNA-seq data from all donors
- `aggregated_atac_data.h5ad`: Aggregated ATAC-seq data from all donors
- `aggregation_summary.csv`: Summary statistics of the aggregation
- `donor_processing_report.csv`: Detailed report of each donor's processing status
- `rna_sample_metadata.csv`: RNA sample metadata
- `atac_sample_metadata.csv`: ATAC sample metadata

Files are saved without compression for easy loading in RStudio on Windows.

### Reporting Features

The aggregation script provides comprehensive reporting:

1. **Console Output**: Real-time status updates and final summary
2. **Donor Processing Report**: Detailed status for each donor including:
   - File discovery status (RNA/ATAC files found/missing)
   - Cell counts (loaded vs. subset)
   - Error messages for failed donors
3. **Summary Statistics**: Overall success rates and file counts
4. **Log Files**: Detailed logging with timestamps

### Separate AnnData Structure

The aggregated data is stored in separate AnnData objects:

```python
# RNA Data (aggregated_rna_data.h5ad)
rna_adata.X                    # RNA-seq data matrix
rna_adata.obs                  # Sample metadata (donor_id, cell_type, etc.)
rna_adata.var                  # RNA gene information

# ATAC Data (aggregated_atac_data.h5ad)
atac_adata.X                   # ATAC-seq data matrix
atac_adata.obs                 # Sample metadata (donor_id, cell_type, etc.)
atac_adata.var                 # ATAC peak information
```

**Note**: Both files contain the same cells (those with multiomic data) and can be easily loaded and analyzed separately or together. This approach ensures data integrity while providing flexibility for independent analysis of each data type.

### Testing

For testing the aggregation workflow with a subset of data, use the test script:

```bash
python scripts/analysis/test_multiomic_aggregation.py \
    --metadata-file data/metadata/multiomic_metadata.csv \
    --rna-data-dir data/rna_files \
    --atac-data-dir data/atac_files \
    --output-dir test_results/aggregated_multiomic \
    --max-donors 3
```

The test script processes only the first N donors (default: 3) and provides the same comprehensive reporting as the full script.

### Example

See `scripts/analysis/example_multiomic_aggregation.py` for a demonstration of the expected data structure and usage.

## Analysis Workflow

1. Data Preprocessing
   - Quality control and filtering
   - Cell type annotation
   - Integration of snRNA-seq and snATAC-seq data

2. Early Vulnerability Analysis
   - CPS-based stratification
   - Cell type abundance analysis
   - Differential expression analysis

3. Multi-modal Integration
   - Gene expression and chromatin accessibility correlation
   - Transcription factor analysis
   - Pathway enrichment

4. Spatial Analysis
   - MERFISH data integration
   - Layer-specific analysis
   - Spatial gene expression patterns

## Data Sources
- SEA-AD snRNA-seq and snATAC-seq data (84 donors)
- MERFISH spatial transcriptomics (27 donors)
- Matched metadata (CPS score, Braak stage, etc.) # SEA_AD_bioinformatic_dashboard
