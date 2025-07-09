#!/usr/bin/env python3
"""
Test script for multiomic aggregation that processes only the first 3 donors.
This script is useful for testing the aggregation workflow with a subset of data.
"""

import pandas as pd
from pathlib import Path
import sys
import logging
from datetime import datetime
from aggregate_multiomic import MultiomicAggregator

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'test_multiomic_aggregation_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

class TestMultiomicAggregator(MultiomicAggregator):
    """
    Test version of MultiomicAggregator that processes only the first N donors.
    """
    
    def __init__(self, metadata_file, rna_data_dir, atac_data_dir, output_dir, max_donors=3):
        """
        Initialize the test multiomic aggregator.
        
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
        max_donors : int
            Maximum number of donors to process (default: 3)
        """
        super().__init__(metadata_file, rna_data_dir, atac_data_dir, output_dir)
        self.max_donors = max_donors
        
    def load_metadata(self):
        """Load and validate the metadata file, limiting to first N donors."""
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
            
            # Get all unique donors
            all_donors = self.metadata['donor_id'].unique()
            logger.info(f"Found {len(all_donors)} total donors in metadata")
            
            # Limit to first N donors
            self.multiomic_donors = all_donors[:self.max_donors]
            logger.info(f"TEST MODE: Processing only first {len(self.multiomic_donors)} donors: {list(self.multiomic_donors)}")
            
            # Filter metadata to only include the selected donors
            self.metadata = self.metadata[self.metadata['donor_id'].isin(self.multiomic_donors)].copy()
            
            # Extract multiomic samples for selected donors
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
    
    def print_processing_summary(self):
        """Print a summary of the processing results to console."""
        print("\n" + "="*60)
        print("TEST MULTIOMIC AGGREGATION PROCESSING SUMMARY")
        print("="*60)
        print(f"TEST MODE: Processed first {self.max_donors} donors only")
        print(f"Donors processed: {list(self.multiomic_donors)}")
        print()
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
        if hasattr(self, 'aggregated_rna') and self.aggregated_rna is not None and hasattr(self, 'aggregated_atac') and self.aggregated_atac is not None:
            print("Final aggregated datasets:")
            print(f"  RNA: {self.aggregated_rna.n_obs} cells, {self.aggregated_rna.n_vars} genes")
            print(f"  ATAC: {self.aggregated_atac.n_obs} cells, {self.aggregated_atac.n_vars} peaks")
            print(f"  Donors: {len(self.aggregated_rna.obs['donor_id'].unique())}")
        else:
            print("Final aggregated datasets: None (aggregation failed)")
        print("="*60)
        print("NOTE: This is a TEST run with limited donors!")
        print("="*60)

def run_test_aggregation(metadata_file, rna_data_dir, atac_data_dir, output_dir, max_donors=3):
    """
    Run a test aggregation with a limited number of donors.
    
    Parameters
    ----------
    metadata_file : str or Path
        Path to the metadata CSV file
    rna_data_dir : str or Path
        Directory containing RNA-seq h5ad files
    atac_data_dir : str or Path
        Directory containing ATAC-seq h5ad files
    output_dir : str or Path
        Directory to save output files
    max_donors : int
        Maximum number of donors to process
    """
    logger.info(f"Starting TEST aggregation with max {max_donors} donors")
    
    # Initialize test aggregator
    test_aggregator = TestMultiomicAggregator(
        metadata_file=metadata_file,
        rna_data_dir=rna_data_dir,
        atac_data_dir=atac_data_dir,
        output_dir=output_dir,
        max_donors=max_donors
    )
    
    # Run the workflow
    try:
        test_aggregator.run_workflow()
        logger.info("TEST aggregation completed successfully!")
        return True
    except Exception as e:
        logger.error(f"Error during TEST aggregation: {str(e)}")
        return False

def main():
    """Main function for running the test aggregation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Test multiomic aggregation with a limited number of donors'
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
        default='test_results/aggregated_multiomic',
        help='Directory to save aggregated output files (default: test_results/aggregated_multiomic)'
    )
    parser.add_argument(
        '--max-donors',
        type=int,
        default=3,
        help='Maximum number of donors to process (default: 3)'
    )
    
    args = parser.parse_args()
    
    # Run the test aggregation
    success = run_test_aggregation(
        metadata_file=args.metadata_file,
        rna_data_dir=args.rna_data_dir,
        atac_data_dir=args.atac_data_dir,
        output_dir=args.output_dir,
        max_donors=args.max_donors
    )
    
    if success:
        print(f"\n✅ TEST aggregation completed successfully!")
        print(f"Results saved to: {args.output_dir}")
    else:
        print(f"\n❌ TEST aggregation failed!")
        sys.exit(1)

if __name__ == "__main__":
    main() 