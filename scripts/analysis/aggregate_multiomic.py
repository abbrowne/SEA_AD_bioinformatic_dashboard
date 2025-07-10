#!/usr/bin/env python3
"""
Multiomic data aggregation script for the SEA-AD project.
Aggregates snRNAseq+snATACseq data from multiple donors into a single AnnData object with multiple assays.
"""

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import logging
import argparse
from datetime import datetime
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'multiomic_aggregation_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class MultiomicAggregator:
    def __init__(self, metadata_file, rna_data_dir, atac_data_dir, output_dir):
        """
        Initialize the multiomic data aggregator.
        
        Parameters
        ----------
        metadata_file : str or Path
            Path to the metadata CSV file containing multiomic sample information
        rna_data_dir : str or Path
            Directory containing RNA-seq h5ad files
        atac_data_dir : str or Path
            Directory containing ATAC-seq h5ad files
        output_dir : str or Path
            Directory to save aggregated output files
        """
        self.metadata_file = Path(metadata_file)
        self.rna_data_dir = Path(rna_data_dir)
        self.atac_data_dir = Path(atac_data_dir)
        self.output_dir = Path(output_dir)
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.metadata = None
        self.multiomic_donors = None
        self.multiomic_samples = None
        self.aggregated_rna = None
        self.aggregated_atac = None
        
        # Initialize reporting containers
        self.donor_status = {}
        self.processing_summary = {
            'total_donors': 0,
            'successful_donors': 0,
            'failed_donors': 0,
            'total_samples': 0,
            'successful_samples': 0,
            'failed_samples': 0,
            'rna_files_found': 0,
            'atac_files_found': 0,
            'rna_files_missing': 0,
            'atac_files_missing': 0
        }
        
    def load_metadata(self):
        """Load and validate the metadata file."""
        logger.info(f"Loading metadata from {self.metadata_file}")
        
        if not self.metadata_file.exists():
            raise FileNotFoundError(f"Metadata file not found at {self.metadata_file}")
        
        try:
            self.metadata = pd.read_csv(self.metadata_file)
            logger.info(f"Loaded metadata with {len(self.metadata)} rows")
            
            # Validate required columns
            required_columns = ['donor_id', 'sample_id']
            missing_columns = [col for col in required_columns if col not in self.metadata.columns]
            if missing_columns:
                raise ValueError(f"Missing required columns in metadata: {missing_columns}")
            
            # Extract multiomic donors and samples
            # All samples in the metadata have multiomic data
            self.multiomic_donors = self.metadata['donor_id'].unique()
            self.multiomic_samples = self.metadata['sample_id'].unique()
            
            # Initialize processing summary
            self.processing_summary['total_donors'] = len(self.multiomic_donors)
            self.processing_summary['total_samples'] = len(self.multiomic_samples)
            
            logger.info(f"Found {len(self.multiomic_donors)} donors with multiomic data")
            logger.info(f"Found {len(self.multiomic_samples)} samples with multiomic data")
            
            # Initialize donor status tracking
            for donor_id in self.multiomic_donors:
                self.donor_status[donor_id] = {
                    'status': 'pending',
                    'rna_file_found': False,
                    'atac_file_found': False,
                    'rna_file_path': None,
                    'atac_file_path': None,
                    'rna_cells_loaded': 0,
                    'atac_cells_loaded': 0,
                    'rna_cells_subset': 0,
                    'atac_cells_subset': 0,
                    'error_message': None
                }
            
        except Exception as e:
            logger.error(f"Error loading metadata: {str(e)}")
            raise
    
    def find_data_files(self, donor_id):
        """
        Find RNA and ATAC data files for a given donor.
        
        Parameters
        ----------
        donor_id : str
            Donor ID to find files for
            
        Returns
        -------
        tuple
            (rna_file_path, atac_file_path) or (None, None) if files not found
        """
        # Look for RNA file
        rna_pattern = f"*{donor_id}*rna*.h5ad"
        rna_files = list(self.rna_data_dir.glob(rna_pattern))
        
        # Look for ATAC file
        atac_pattern = f"*{donor_id}*atac*.h5ad"
        atac_files = list(self.atac_data_dir.glob(atac_pattern))
        
        # Update donor status
        if rna_files:
            self.donor_status[donor_id]['rna_file_found'] = True
            self.donor_status[donor_id]['rna_file_path'] = str(rna_files[0])
            self.processing_summary['rna_files_found'] += 1
        else:
            self.donor_status[donor_id]['rna_file_found'] = False
            self.processing_summary['rna_files_missing'] += 1
            
        if atac_files:
            self.donor_status[donor_id]['atac_file_found'] = True
            self.donor_status[donor_id]['atac_file_path'] = str(atac_files[0])
            self.processing_summary['atac_files_found'] += 1
        else:
            self.donor_status[donor_id]['atac_file_found'] = False
            self.processing_summary['atac_files_missing'] += 1
        
        if not rna_files:
            logger.warning(f"No RNA file found for donor {donor_id}")
            return None, None
        
        if not atac_files:
            logger.warning(f"No ATAC file found for donor {donor_id}")
            return None, None
        
        # Use the first matching file for each type
        rna_file = rna_files[0]
        atac_file = atac_files[0]
        
        logger.info(f"Found files for donor {donor_id}:")
        logger.info(f"  RNA: {rna_file}")
        logger.info(f"  ATAC: {atac_file}")
        
        return rna_file, atac_file
    
    def load_donor_data(self, donor_id):
        """
        Load RNA and ATAC data for a specific donor.
        
        Parameters
        ----------
        donor_id : str
            Donor ID to load data for
            
        Returns
        -------
        tuple
            (rna_adata, atac_adata) or (None, None) if loading fails
        """
        rna_file, atac_file = self.find_data_files(donor_id)
        
        if rna_file is None or atac_file is None:
            self.donor_status[donor_id]['status'] = 'failed'
            self.donor_status[donor_id]['error_message'] = 'Missing data files'
            self.processing_summary['failed_donors'] += 1
            return None, None
        
        try:
            # Load RNA data
            logger.info(f"Loading RNA data for donor {donor_id}")
            rna_adata = sc.read_h5ad(rna_file)
            self.donor_status[donor_id]['rna_cells_loaded'] = rna_adata.n_obs
            
            # Load ATAC data
            logger.info(f"Loading ATAC data for donor {donor_id}")
            atac_adata = sc.read_h5ad(atac_file)
            self.donor_status[donor_id]['atac_cells_loaded'] = atac_adata.n_obs
            
            # Add donor information
            rna_adata.obs['donor_id'] = donor_id
            atac_adata.obs['donor_id'] = donor_id
            
            logger.info(f"Loaded {rna_adata.n_obs} RNA cells and {atac_adata.n_obs} ATAC cells for donor {donor_id}")
            
            return rna_adata, atac_adata
            
        except Exception as e:
            error_msg = f"Error loading data for donor {donor_id}: {str(e)}"
            logger.error(error_msg)
            self.donor_status[donor_id]['status'] = 'failed'
            self.donor_status[donor_id]['error_message'] = error_msg
            self.processing_summary['failed_donors'] += 1
            return None, None
    
    def subset_to_multiomic_samples(self, adata, donor_id):
        """
        Subset data to only samples with multiomic data for a specific donor.
        
        Parameters
        ----------
        adata : AnnData
            Single-cell dataset
        donor_id : str
            Donor ID to subset for
            
        Returns
        -------
        AnnData
            Subsetted dataset
        """
        # Get multiomic samples for this donor
        donor_samples = self.metadata[
            self.metadata['donor_id'] == donor_id
        ]['sample_id'].unique()
        
        # Find cells in the data that match multiomic sample IDs
        # Assuming sample IDs are in the index or a specific column
        if 'sample_id' in adata.obs.columns:
            sample_id_col = 'sample_id'
        else:
            # If no sample_id column, assume the index contains sample IDs
            sample_id_col = None
        
        if sample_id_col:
            multiomic_mask = adata.obs[sample_id_col].isin(donor_samples)
        else:
            multiomic_mask = adata.obs_names.isin(donor_samples)
        
        subset_adata = adata[multiomic_mask].copy()
        
        # Update donor status
        if 'rna' in adata.var_names[0] if len(adata.var_names) > 0 else False:
            self.donor_status[donor_id]['rna_cells_subset'] = subset_adata.n_obs
        else:
            self.donor_status[donor_id]['atac_cells_subset'] = subset_adata.n_obs
        
        logger.info(f"Subset donor {donor_id} to {subset_adata.n_obs} multiomic samples")
        
        return subset_adata
    
    def aggregate_data(self):
        """Aggregate multiomic data from all donors into separate RNA and ATAC AnnData objects."""
        logger.info("Starting data aggregation...")
        
        aggregated_rna_list = []
        aggregated_atac_list = []
        
        for i, donor_id in enumerate(self.multiomic_donors):
            logger.info(f"Processing donor {i+1}/{len(self.multiomic_donors)}: {donor_id}")
            
            # Load donor data
            rna_adata, atac_adata = self.load_donor_data(donor_id)
            
            if rna_adata is None or atac_adata is None:
                logger.warning(f"Skipping donor {donor_id} due to missing data files")
                continue
            
            # Subset to multiomic samples
            rna_subset = self.subset_to_multiomic_samples(rna_adata, donor_id)
            atac_subset = self.subset_to_multiomic_samples(atac_adata, donor_id)
            
            if rna_subset.n_obs == 0 or atac_subset.n_obs == 0:
                logger.warning(f"No multiomic samples found for donor {donor_id}")
                self.donor_status[donor_id]['status'] = 'failed'
                self.donor_status[donor_id]['error_message'] = 'No multiomic samples found after subsetting'
                self.processing_summary['failed_donors'] += 1
                continue
            
            # Ensure both datasets have the same cells
            common_cells = set(rna_subset.obs_names) & set(atac_subset.obs_names)
            if len(common_cells) == 0:
                logger.warning(f"No common cells found between RNA and ATAC for donor {donor_id}")
                self.donor_status[donor_id]['status'] = 'failed'
                self.donor_status[donor_id]['error_message'] = 'No common cells between RNA and ATAC'
                self.processing_summary['failed_donors'] += 1
                continue
            
            # Subset to common cells
            rna_subset = rna_subset[list(common_cells)].copy()
            atac_subset = atac_subset[list(common_cells)].copy()
            
            logger.info(f"Final subset for donor {donor_id}: {rna_subset.n_obs} cells with both RNA and ATAC data")
            
            # Add to aggregation lists
            aggregated_rna_list.append(rna_subset)
            aggregated_atac_list.append(atac_subset)
            
            # Update donor status
            self.donor_status[donor_id]['status'] = 'success'
            self.processing_summary['successful_donors'] += 1
            self.processing_summary['successful_samples'] += rna_subset.n_obs
        
        if not aggregated_rna_list:
            raise ValueError("No valid RNA data found for any donor")
        
        if not aggregated_atac_list:
            raise ValueError("No valid ATAC data found for any donor")
        
        # Concatenate all datasets
        logger.info("Concatenating RNA datasets...")
        self.aggregated_rna = sc.concat(aggregated_rna_list, join='outer', index_unique=None)
        
        logger.info("Concatenating ATAC datasets...")
        self.aggregated_atac = sc.concat(aggregated_atac_list, join='outer', index_unique=None)
        
        # Ensure both datasets have exactly the same cells
        common_cells = set(self.aggregated_rna.obs_names) & set(self.aggregated_atac.obs_names)
        if len(common_cells) != len(self.aggregated_rna.obs_names) or len(common_cells) != len(self.aggregated_atac.obs_names):
            logger.warning(f"Cell count mismatch after concatenation. RNA: {len(self.aggregated_rna.obs_names)}, ATAC: {len(self.aggregated_atac.obs_names)}, Common: {len(common_cells)}")
            logger.info("Subsetting to common cells only...")
            self.aggregated_rna = self.aggregated_rna[list(common_cells)].copy()
            self.aggregated_atac = self.aggregated_atac[list(common_cells)].copy()
        
        logger.info(f"Aggregation complete:")
        logger.info(f"  RNA dataset: {self.aggregated_rna.n_obs} cells, {self.aggregated_rna.n_vars} genes")
        logger.info(f"  ATAC dataset: {self.aggregated_atac.n_obs} cells, {self.aggregated_atac.n_vars} peaks")
        logger.info(f"  Donors: {len(self.aggregated_rna.obs['donor_id'].unique())}")
    
    def generate_processing_report(self):
        """Generate a comprehensive processing report."""
        logger.info("Generating processing report...")
        
        # Create detailed donor report
        donor_report = []
        for donor_id, status in self.donor_status.items():
            donor_report.append({
                'donor_id': donor_id,
                'status': status['status'],
                'rna_file_found': status['rna_file_found'],
                'atac_file_found': status['atac_file_found'],
                'rna_file_path': status['rna_file_path'],
                'atac_file_path': status['atac_file_path'],
                'rna_cells_loaded': status['rna_cells_loaded'],
                'atac_cells_loaded': status['atac_cells_loaded'],
                'rna_cells_subset': status['rna_cells_subset'],
                'atac_cells_subset': status['atac_cells_subset'],
                'error_message': status['error_message']
            })
        
        donor_report_df = pd.DataFrame(donor_report)
        
        # Create summary report
        summary_report = {
            'metric': [
                'total_donors',
                'successful_donors', 
                'failed_donors',
                'total_samples',
                'successful_samples',
                'failed_samples',
                'rna_files_found',
                'atac_files_found',
                'rna_files_missing',
                'atac_files_missing',
                'success_rate_donors',
                'success_rate_samples'
            ],
            'value': [
                self.processing_summary['total_donors'],
                self.processing_summary['successful_donors'],
                self.processing_summary['failed_donors'],
                self.processing_summary['total_samples'],
                self.processing_summary['successful_samples'],
                self.processing_summary['failed_samples'],
                self.processing_summary['rna_files_found'],
                self.processing_summary['atac_files_found'],
                self.processing_summary['rna_files_missing'],
                self.processing_summary['atac_files_missing'],
                f"{self.processing_summary['successful_donors']/self.processing_summary['total_donors']*100:.1f}%" if self.processing_summary['total_donors'] > 0 else "0%",
                f"{self.processing_summary['successful_samples']/self.processing_summary['total_samples']*100:.1f}%" if self.processing_summary['total_samples'] > 0 else "0%"
            ]
        }
        
        summary_report_df = pd.DataFrame(summary_report)
        
        return donor_report_df, summary_report_df
    
    def prepare_for_seurat(self, adata, data_type):
        """
        Prepare AnnData object for Seurat compatibility.
        
        Parameters
        ----------
        adata : AnnData
            AnnData object to prepare
        data_type : str
            Type of data ('rna' or 'atac')
        """
        logger.info(f"Preparing {data_type.upper()} data for Seurat compatibility...")
        
        # Ensure X contains raw counts (not normalized data)
        if data_type == 'rna':
            # For RNA, ensure we have raw counts
            if 'counts' in adata.layers:
                adata.X = adata.layers['counts'].copy()
                logger.info("Using raw counts from layers['counts']")
            else:
                logger.info("Using X matrix as counts (ensure this contains raw counts)")
        else:
            # For ATAC, ensure we have raw counts
            if 'counts' in adata.layers:
                adata.X = adata.layers['counts'].copy()
                logger.info("Using raw counts from layers['counts']")
            else:
                logger.info("Using X matrix as counts (ensure this contains raw counts)")
        
        # Ensure obs_names are unique and properly formatted
        adata.obs_names = adata.obs_names.astype(str)
        if not adata.obs_names.is_unique:
            logger.warning("Non-unique cell names found, making them unique...")
            adata.obs_names = [f"{name}_{i}" for i, name in enumerate(adata.obs_names)]
        
        # Ensure var_names are unique and properly formatted
        adata.var_names = adata.var_names.astype(str)
        if not adata.var_names.is_unique:
            logger.warning("Non-unique feature names found, making them unique...")
            adata.var_names = [f"{name}_{i}" for i, name in enumerate(adata.var_names)]
        
        # Clean up obs metadata - remove problematic columns
        problematic_cols = []
        for col in adata.obs.columns:
            # Check for columns with complex data types
            if adata.obs[col].dtype == 'object':
                # Check if column contains lists or other complex objects
                sample_val = adata.obs[col].iloc[0] if len(adata.obs) > 0 else None
                if isinstance(sample_val, (list, dict, tuple)):
                    problematic_cols.append(col)
                    logger.warning(f"Removing problematic column '{col}' (contains complex data)")
        
        if problematic_cols:
            adata.obs = adata.obs.drop(columns=problematic_cols)
        
        # Ensure all obs columns are simple data types
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                # Convert to string if possible
                try:
                    adata.obs[col] = adata.obs[col].astype(str)
                except:
                    logger.warning(f"Could not convert column '{col}' to string, removing...")
                    adata.obs = adata.obs.drop(columns=[col])
        
        # Add required metadata for Seurat
        adata.uns['seurat_version'] = '5.0.0'
        adata.uns['data_type'] = data_type
        
        # Ensure var has proper structure
        if 'gene_ids' not in adata.var.columns:
            adata.var['gene_ids'] = adata.var_names
        if 'feature_types' not in adata.var.columns:
            adata.var['feature_types'] = 'Gene Expression' if data_type == 'rna' else 'Peaks'
        
        logger.info(f"{data_type.upper()} data prepared for Seurat conversion")
        return adata
    
    def clean_anndata_for_seurat(self, adata):
        """
        Clean AnnData object for Seurat compatibility by converting all obs/var columns to strings if categorical or object type.
        """
        import pandas as pd
        logger.info("Cleaning AnnData object for Seurat compatibility...")
        for col in adata.obs.columns:
            if pd.api.types.is_categorical_dtype(adata.obs[col]) or adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype(str)
        for col in adata.var.columns:
            if pd.api.types.is_categorical_dtype(adata.var[col]) or adata.var[col].dtype == object:
                adata.var[col] = adata.var[col].astype(str)
        logger.info("AnnData object cleaned for Seurat compatibility.")
        return adata

    def save_aggregated_data(self):
        """Save aggregated data as separate RNA and ATAC files without compression for easy loading in R."""
        logger.info("Saving aggregated data...")
        
        # Clean AnnData objects for Seurat compatibility
        self.aggregated_rna = self.clean_anndata_for_seurat(self.aggregated_rna)
        self.aggregated_atac = self.clean_anndata_for_seurat(self.aggregated_atac)
        
        # Prepare data for Seurat compatibility
        rna_adata_seurat = self.prepare_for_seurat(self.aggregated_rna, 'rna')
        atac_adata_seurat = self.prepare_for_seurat(self.aggregated_atac, 'atac')
        
        # Save RNA data
        rna_output_file = self.output_dir / "aggregated_rna_data.h5ad"
        logger.info(f"Saving RNA data to {rna_output_file}")
        rna_adata_seurat.write(rna_output_file, compression=None)
        
        # Save ATAC data
        atac_output_file = self.output_dir / "aggregated_atac_data.h5ad"
        logger.info(f"Saving ATAC data to {atac_output_file}")
        atac_adata_seurat.write(atac_output_file, compression=None)
        
        # Generate and save reports
        donor_report_df, summary_report_df = self.generate_processing_report()
        
        # Save donor processing report
        donor_report_file = self.output_dir / "donor_processing_report.csv"
        donor_report_df.to_csv(donor_report_file, index=False)
        logger.info(f"Saved donor processing report to {donor_report_file}")
        
        # Save summary report
        summary_file = self.output_dir / "aggregation_summary.csv"
        summary_report_df.to_csv(summary_file, index=False)
        logger.info(f"Saved summary to {summary_file}")
        
        # Save sample metadata for RNA
        rna_sample_metadata_file = self.output_dir / "rna_sample_metadata.csv"
        rna_adata_seurat.obs.to_csv(rna_sample_metadata_file)
        logger.info(f"Saved RNA sample metadata to {rna_sample_metadata_file}")
        
        # Save sample metadata for ATAC
        atac_sample_metadata_file = self.output_dir / "atac_sample_metadata.csv"
        atac_adata_seurat.obs.to_csv(atac_sample_metadata_file)
        logger.info(f"Saved ATAC sample metadata to {atac_sample_metadata_file}")
        
        # Print summary to console
        self.print_processing_summary()
        
        logger.info("Data saving complete!")
    
    def print_processing_summary(self):
        """Print a summary of the processing results to console."""
        print("\n" + "="*60)
        print("MULTIOMIC AGGREGATION PROCESSING SUMMARY")
        print("="*60)
        print(f"Total donors processed: {self.processing_summary['total_donors']}")
        print(f"Successful donors: {self.processing_summary['successful_donors']}")
        print(f"Failed donors: {self.processing_summary['failed_donors']}")
        print(f"Success rate: {self.processing_summary['successful_donors']/self.processing_summary['total_donors']*100:.1f}%" if self.processing_summary['total_donors'] > 0 else "0%")
        print()
        print(f"Total samples: {self.processing_summary['total_samples']}")
        print(f"Successful samples: {self.processing_summary['successful_samples']}")
        print(f"Failed samples: {self.processing_summary['failed_samples']}")
        print()
        print(f"RNA files found: {self.processing_summary['rna_files_found']}")
        print(f"ATAC files found: {self.processing_summary['atac_files_found']}")
        print(f"RNA files missing: {self.processing_summary['rna_files_missing']}")
        print(f"ATAC files missing: {self.processing_summary['atac_files_missing']}")
        print()
        print("Final aggregated datasets:")
        print(f"  RNA: {self.aggregated_rna.n_obs} cells, {self.aggregated_rna.n_vars} genes")
        print(f"  ATAC: {self.aggregated_atac.n_obs} cells, {self.aggregated_atac.n_vars} peaks")
        print(f"  Donors: {len(self.aggregated_rna.obs['donor_id'].unique())}")
        print("="*60)
    
    def run_workflow(self):
        """Run the complete multiomic aggregation workflow."""
        try:
            self.load_metadata()
            self.aggregate_data()
            self.save_aggregated_data()
            
            logger.info("Multiomic aggregation workflow completed successfully!")
            
        except Exception as e:
            logger.error(f"Error in multiomic aggregation workflow: {str(e)}")
            raise

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Aggregate multiomic snRNAseq+snATACseq data from multiple donors into a single AnnData object'
    )
    parser.add_argument(
        '--metadata-file',
        required=True,
        help='Path to the metadata CSV file containing multiomic sample information'
    )
    parser.add_argument(
        '--rna-data-dir',
        required=True,
        help='Directory containing RNA-seq h5ad files'
    )
    parser.add_argument(
        '--atac-data-dir',
        required=True,
        help='Directory containing ATAC-seq h5ad files'
    )
    parser.add_argument(
        '--output-dir',
        default='results/aggregated_multiomic',
        help='Directory to save aggregated output files (default: results/aggregated_multiomic)'
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    # Initialize and run aggregation
    aggregator = MultiomicAggregator(
        metadata_file=args.metadata_file,
        rna_data_dir=args.rna_data_dir,
        atac_data_dir=args.atac_data_dir,
        output_dir=args.output_dir
    )
    
    aggregator.run_workflow() 