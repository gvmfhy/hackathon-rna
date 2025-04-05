#!/usr/bin/env python3
"""
TCGA LIHC miRNA Data Acquisition from UCSC Xena Browser

This script downloads and processes miRNA expression data and clinical metadata
for NASH-HCC vs Normal Liver comparison from the UCSC Xena Browser.

Usage:
    python fetch_xena_data.py --output-dir data

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import requests
import gzip
import io

# Add src directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

# Xena Browser API URLs
XENA_API_BASE = "https://xena.ucsc.edu/api/"
XENA_HUB = "https://tcga.xenahubs.net"
XENA_DATASET_MIRNA = "tcga_LIHC/Xena_Matrices/tcga_LIHC_miRNA_gene_expression.txt"
XENA_DATASET_CLINICAL = "tcga_LIHC/Xena_Matrices/tcga_LIHC_clinicalMatrix"

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Download TCGA LIHC miRNA data from UCSC Xena Browser")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="data",
                      help="Output directory for downloaded data (default: data)")
    parser.add_argument("--mirna-filename", type=str, default="TCGA_LIHC_miRNA_expression.tsv",
                      help="Filename for miRNA expression data (default: TCGA_LIHC_miRNA_expression.tsv)")
    parser.add_argument("--clinical-filename", type=str, default="TCGA_LIHC_clinical.tsv",
                      help="Filename for clinical data (default: TCGA_LIHC_clinical.tsv)")
    parser.add_argument("--processed-filename", type=str, default="TCGA_LIHC_miRNA_processed.csv",
                      help="Filename for processed data (default: TCGA_LIHC_miRNA_processed.csv)")
    
    return parser.parse_args()

def download_xena_data(url, dataset):
    """
    Download data from UCSC Xena Browser.
    
    Args:
        url: Xena API base URL
        dataset: Dataset path
        
    Returns:
        DataFrame with downloaded data
    """
    print(f"Downloading data from {dataset}...")
    
    full_url = f"{url}/download/{dataset}"
    response = requests.get(full_url)
    
    if response.status_code != 200:
        print(f"Error downloading data: HTTP {response.status_code}")
        sys.exit(1)
    
    # Read the tsv data
    try:
        if dataset.endswith('.gz'):
            with gzip.open(io.BytesIO(response.content), 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            df = pd.read_csv(io.StringIO(response.text), sep='\t')
        
        print(f"Successfully downloaded data with shape: {df.shape}")
        return df
    except Exception as e:
        print(f"Error parsing downloaded data: {str(e)}")
        sys.exit(1)

def process_xena_data(mirna_data, clinical_data):
    """
    Process and merge miRNA expression and clinical data.
    
    Args:
        mirna_data: miRNA expression DataFrame
        clinical_data: Clinical data DataFrame
        
    Returns:
        Processed DataFrame ready for analysis
    """
    print("Processing and merging datasets...")
    
    # Extract sample type from clinical data
    if 'sample_type' in clinical_data.columns:
        sample_type_col = 'sample_type'
    elif 'sample_type_id' in clinical_data.columns:
        sample_type_col = 'sample_type_id'
    else:
        print("Warning: Could not find sample type column in clinical data.")
        print(f"Available columns: {', '.join(clinical_data.columns)}")
        # Try to infer from sample names - TCGA convention is that normal samples end with -11
        clinical_data['inferred_sample_type'] = clinical_data.index.map(
            lambda x: 'Solid Tissue Normal' if x.endswith('-11') else 'Primary Tumor'
        )
        sample_type_col = 'inferred_sample_type'
    
    # Create group column
    clinical_subset = clinical_data[[sample_type_col]].copy()
    clinical_subset['Group'] = clinical_subset[sample_type_col].map(
        lambda x: 'Normal' if x == 'Solid Tissue Normal' else 'NASH-HCC'
    )
    
    # Prepare miRNA data (transpose if needed)
    if 'sample' in mirna_data.columns:
        # Data is already in the right format (miRNAs as rows)
        mirna_processed = mirna_data.copy()
        id_col = 'sample'
    else:
        # Transpose data (common Xena format has samples as columns)
        mirna_processed = mirna_data.set_index(mirna_data.columns[0]).T
        mirna_processed.index.name = 'sample'
        mirna_processed = mirna_processed.reset_index()
    
    # Merge datasets
    merged_data = pd.merge(
        mirna_processed,
        clinical_subset,
        left_on='sample',
        right_index=True,
        how='inner'
    )
    
    print(f"Processed data shape: {merged_data.shape}")
    print(f"Sample groups: {merged_data['Group'].value_counts().to_dict()}")
    
    return merged_data

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Download data
    mirna_data = download_xena_data(XENA_HUB, XENA_DATASET_MIRNA)
    clinical_data = download_xena_data(XENA_HUB, XENA_DATASET_CLINICAL)
    
    # Save raw data
    mirna_output_path = output_dir / args.mirna_filename
    clinical_output_path = output_dir / args.clinical_filename
    
    mirna_data.to_csv(mirna_output_path, sep='\t', index=False)
    clinical_data.to_csv(clinical_output_path, sep='\t')
    
    print(f"Raw miRNA data saved to {mirna_output_path}")
    print(f"Raw clinical data saved to {clinical_output_path}")
    
    # Process data
    processed_data = process_xena_data(mirna_data, clinical_data)
    
    # Save processed data
    processed_output_path = output_dir / args.processed_filename
    processed_data.to_csv(processed_output_path, index=False)
    
    print(f"Processed data saved to {processed_output_path}")
    print("\nData acquisition complete. You can now run the target_prioritizer.py script with this data.")
    print(f"Example: python target_prioritizer.py --input {processed_output_path} --id-col sample")

if __name__ == "__main__":
    main() 