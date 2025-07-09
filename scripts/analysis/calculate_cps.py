#!/usr/bin/env python3
"""
Calculate Continuous Disease Pseudoprogression Score (CPS) from metadata.

The CPS score combines multiple measures of AD progression:
1. Braak stage (0-6)
2. Thal phase (0-5)
3. CERAD score (0-3)
4. Cognitive metrics (if available):
   - MMSE (Mini-Mental State Examination, 0-30)
   - CDR (Clinical Dementia Rating, 0-3)
   - MoCA (Montreal Cognitive Assessment, 0-30)
   - ADAS-Cog (Alzheimer's Disease Assessment Scale-Cognitive, 0-70)

The score is normalized to a 0-1 scale where:
0 = No pathology
1 = Severe pathology
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import argparse
from sklearn.preprocessing import MinMaxScaler

def normalize_cognitive_metrics(metadata_df):
    """
    Normalize cognitive metrics to a 0-1 scale where higher values indicate worse cognition.
    
    Parameters
    ----------
    metadata_df : pandas.DataFrame
        DataFrame containing cognitive metrics
    
    Returns
    -------
    numpy.ndarray
        Normalized cognitive score
    """
    scaler = MinMaxScaler()
    cognitive_scores = []
    
    # MMSE (0-30, lower is worse)
    if 'mmse_score' in metadata_df.columns:
        mmse_norm = 1 - scaler.fit_transform(metadata_df[['mmse_score']])[:, 0]  # Invert scale
        cognitive_scores.append(mmse_norm)
    
    # CDR (0-3, higher is worse)
    if 'cdr_score' in metadata_df.columns:
        cdr_norm = scaler.fit_transform(metadata_df[['cdr_score']])[:, 0]
        cognitive_scores.append(cdr_norm)
    
    # MoCA (0-30, lower is worse)
    if 'moca_score' in metadata_df.columns:
        moca_norm = 1 - scaler.fit_transform(metadata_df[['moca_score']])[:, 0]  # Invert scale
        cognitive_scores.append(moca_norm)
    
    # ADAS-Cog (0-70, higher is worse)
    if 'adas_cog_score' in metadata_df.columns:
        adas_norm = scaler.fit_transform(metadata_df[['adas_cog_score']])[:, 0]
        cognitive_scores.append(adas_norm)
    
    if not cognitive_scores:
        return None
    
    # Average all available cognitive scores
    return np.mean(cognitive_scores, axis=0)

def calculate_cps(metadata_df):
    """
    Calculate CPS score from metadata.
    
    Parameters
    ----------
    metadata_df : pandas.DataFrame
        DataFrame containing the following columns:
        - braak_stage: Braak stage (0-6)
        - thal_phase: Thal phase (0-5)
        - cerad_score: CERAD score (0-3)
        - Optional cognitive metrics:
          - mmse_score: MMSE score (0-30)
          - cdr_score: CDR score (0-3)
          - moca_score: MoCA score (0-30)
          - adas_cog_score: ADAS-Cog score (0-70)
    
    Returns
    -------
    pandas.DataFrame
        Original metadata with added cps_score column
    """
    # Required columns
    required_cols = ['braak_stage', 'thal_phase', 'cerad_score']
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")
    
    # Normalize each component to 0-1 scale
    scaler = MinMaxScaler()
    
    # Braak stage (0-6)
    braak_norm = scaler.fit_transform(metadata_df[['braak_stage']])[:, 0]
    
    # Thal phase (0-5)
    thal_norm = scaler.fit_transform(metadata_df[['thal_phase']])[:, 0]
    
    # CERAD score (0-3)
    cerad_norm = scaler.fit_transform(metadata_df[['cerad_score']])[:, 0]
    
    # Calculate composite score
    # Weights can be adjusted based on your specific needs
    weights = {
        'braak': 0.4,
        'thal': 0.3,
        'cerad': 0.3
    }
    
    cps_score = (
        weights['braak'] * braak_norm +
        weights['thal'] * thal_norm +
        weights['cerad'] * cerad_norm
    )
    
    # Add cognitive metrics if available
    #cog_norm = normalize_cognitive_metrics(metadata_df)
    cog_norm = None
    if cog_norm is not None:
        # Adjust weights to include cognitive metrics
        weights = {
            'braak': 0.3,
            'thal': 0.2,
            'cerad': 0.2,
            'cognitive': 0.3
        }
        cps_score = (
            weights['braak'] * braak_norm +
            weights['thal'] * thal_norm +
            weights['cerad'] * cerad_norm +
            weights['cognitive'] * cog_norm
        )
    
    # Add CPS score to metadata
    metadata_df['cps_score'] = cps_score
    
    return metadata_df

def main():
    parser = argparse.ArgumentParser(
        description='Calculate Continuous Disease Pseudoprogression Score (CPS)'
    )
    parser.add_argument(
        '--input-file',
        required=True,
        help='Path to input metadata file (CSV format)'
    )
    parser.add_argument(
        '--output-file',
        required=True,
        help='Path to output metadata file with CPS scores (CSV format)'
    )
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        # Read metadata
        logging.info(f"Reading metadata from {args.input_file}")
        metadata = pd.read_csv(args.input_file)
        
        # Calculate CPS scores
        logging.info("Calculating CPS scores...")
        metadata_with_cps = calculate_cps(metadata)
        
        # Save results
        logging.info(f"Saving results to {args.output_file}")
        metadata_with_cps.to_csv(args.output_file, index=False)
        
        # Print summary statistics
        logging.info("\nCPS Score Summary:")
        logging.info(f"Mean: {metadata_with_cps['cps_score'].mean():.3f}")
        logging.info(f"Std: {metadata_with_cps['cps_score'].std():.3f}")
        logging.info(f"Min: {metadata_with_cps['cps_score'].min():.3f}")
        logging.info(f"Max: {metadata_with_cps['cps_score'].max():.3f}")
        
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    main() 