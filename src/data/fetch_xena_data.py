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
import time

# Add src directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

# Xena Browser API URLs
# Direct download links for TCGA LIHC data
XENA_MIRNA_URL = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FHiSeqmiRNAseq.gz"
XENA_CLINICAL_URL = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix.gz"

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

def download_data(url, output_path):
    """
    Download data from URL.
    
    Args:
        url: URL to download data from
        output_path: Path to save downloaded data
        
    Returns:
        DataFrame with downloaded data
    """
    print(f"Downloading data from {url}...")
    
    try:
        response = requests.get(url, stream=True)
        
        if response.status_code != 200:
            print(f"Error downloading data: HTTP {response.status_code}")
            sys.exit(1)
        
        # Check if gzipped
        if url.endswith('.gz'):
            with gzip.open(io.BytesIO(response.content), 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            df = pd.read_csv(io.StringIO(response.text), sep='\t')
        
        # Save raw data
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Data saved to {output_path}")
        
        return df
    except Exception as e:
        print(f"Error downloading/processing data: {str(e)}")
        sys.exit(1)

def process_data(mirna_data, clinical_data, output_path):
    """
    Process miRNA and clinical data.
    
    Args:
        mirna_data: miRNA expression DataFrame
        clinical_data: Clinical data DataFrame
        output_path: Path to save processed data
        
    Returns:
        Processed DataFrame
    """
    print("Processing data...")
    
    # Identify sample types from clinical data
    if 'sample_type' in clinical_data.columns:
        sample_type_col = 'sample_type'
    elif '_sample_type' in clinical_data.columns:
        sample_type_col = '_sample_type'
    else:
        # Try to check column names to find the right one
        print("Available clinical columns:")
        for col in clinical_data.columns:
            print(f"  - {col}")
            
        # Try to infer based on sample ID (TCGA convention)
        clinical_data['inferred_sample_type'] = clinical_data.index.map(
            lambda x: 'Solid Tissue Normal' if str(x).endswith('-11') else 'Primary Tumor'
        )
        sample_type_col = 'inferred_sample_type'
    
    # Create group column
    print(f"Using column '{sample_type_col}' to identify sample types.")
    clinical_data['Group'] = clinical_data[sample_type_col].map(
        lambda x: 'Normal' if 'normal' in str(x).lower() else 'NASH-HCC'
    )
    
    # Print sample counts
    sample_counts = clinical_data['Group'].value_counts()
    print(f"Sample group counts: {dict(sample_counts)}")
    
    # Prepare miRNA data
    # Check if transposition is needed based on data shape
    print(f"miRNA data shape: {mirna_data.shape}")
    
    if mirna_data.shape[0] > mirna_data.shape[1]:
        # More rows than columns, likely miRNAs as rows
        print("Data format: miRNAs as rows - transposing...")
        mirna_processed = mirna_data.set_index(mirna_data.columns[0]).T
    else:
        # More columns than rows, likely samples as rows
        print("Data format: Samples as rows - keeping as is...")
        mirna_processed = mirna_data
    
    # Ensure indices match
    common_samples = set(mirna_processed.index).intersection(set(clinical_data.index))
    print(f"Found {len(common_samples)} common samples between datasets.")
    
    if len(common_samples) == 0:
        print("Error: No common samples between miRNA and clinical data.")
        print("miRNA data sample IDs (first 5):", list(mirna_processed.index)[:5])
        print("Clinical data sample IDs (first 5):", list(clinical_data.index)[:5])
        sys.exit(1)
    
    # Merge datasets
    merged_data = pd.merge(
        mirna_processed.reset_index().rename(columns={'index': 'Sample_ID'}),
        clinical_data[['Group']].reset_index().rename(columns={'index': 'Sample_ID'}),
        on='Sample_ID',
        how='inner'
    )
    
    # Save processed data
    merged_data.to_csv(output_path, index=False)
    print(f"Processed data saved to {output_path}")
    
    return merged_data

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Download data
    mirna_path = output_dir / args.mirna_filename
    clinical_path = output_dir / args.clinical_filename
    processed_path = output_dir / args.processed_filename
    
    # Download and save raw data
    mirna_data = download_data(XENA_MIRNA_URL, mirna_path)
    
    # Add a small delay before downloading clinical data to avoid rate limiting
    time.sleep(1)
    
    clinical_data = download_data(XENA_CLINICAL_URL, clinical_path)
    
    # Process data
    process_data(mirna_data, clinical_data, processed_path)
    
    print("\nData acquisition complete. You can now run the target_prioritizer.py script with this data.")
    print(f"Example: python target_prioritizer.py --input {processed_path} --id-col Sample_ID")

if __name__ == "__main__":
    main() 