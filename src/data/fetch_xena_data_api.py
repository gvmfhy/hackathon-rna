#!/usr/bin/env python3
"""
Download miRNA expression data and clinical metadata from UCSC Xena Browser.
Focus on NASH-HCC vs Normal Liver comparison.

Usage:
    python3 fetch_xena_data_api.py --output-dir data

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import time

# Import our SSL fix before importing xenaPython
try:
    import xena_ssl_fix
except ImportError:
    print("Warning: xena_ssl_fix module not found. SSL certificate errors may occur.")

# Import xenaPython after the SSL fix
import xenaPython as xena

# Add src directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

# Constants
XENA_HUB = "https://tcga.xenahubs.net"
# Corrected miRNA dataset identifier based on test_xena_api.py results
XENA_MIRNA_DATASET = "TCGA.LIHC.sampleMap/miRNA_HiSeq_gene"
XENA_CLINICAL_DATASET = "TCGA.LIHC.sampleMap/LIHC_clinicalMatrix"

# Batch size for downloading data
BATCH_SIZE = 100

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Download TCGA LIHC miRNA data using xenaPython API")
    
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

def download_expression_data(hub, dataset, output_path):
    """
    Download expression data.
    
    Args:
        hub: Xena hub URL
        dataset: Dataset name
        output_path: Path to save downloaded data
        
    Returns:
        DataFrame with expression data
    """
    # Check if file already exists
    if os.path.exists(output_path):
        print(f"Expression data file already exists at {output_path}. Loading from file...")
        return pd.read_csv(output_path, sep='\t', index_col=0)
    
    print(f"Downloading miRNA expression data from {dataset}...")
    
    try:
        # Get all samples in the dataset
        samples = xena.dataset_samples(hub, dataset, None)
        print(f"Found {len(samples)} samples in dataset")
        
        # Get all probes (miRNAs) in the dataset
        probes = xena.dataset_field(hub, dataset)
        print(f"Found {len(probes)} miRNAs in dataset")
        
        # Download data in batches to avoid timeouts
        batch_size = 50
        all_data = {}
        
        for i in range(0, len(probes), batch_size):
            probe_batch = probes[i:i+batch_size]
            print(f"Downloading batch {i//batch_size + 1}/{(len(probes)-1)//batch_size + 1} ({len(probe_batch)} miRNAs)")
            
            # Get expression values
            [positions, data_vals] = xena.dataset_probe_values(hub, dataset, samples, probe_batch)
            
            # Map values to probes
            for j, probe in enumerate(probe_batch):
                all_data[probe] = data_vals[j]
            
            # Sleep to avoid overwhelming the server
            time.sleep(0.5)
        
        # Create DataFrame
        df = pd.DataFrame(all_data, index=samples)
        
        # Save to file
        df.to_csv(output_path, sep='\t')
        print(f"Expression data saved to {output_path}")
        
        return df
    except Exception as e:
        print(f"Error downloading expression data: {str(e)}")
        sys.exit(1)

def download_clinical_data(hub, dataset, output_path):
    """
    Download clinical data.
    
    Args:
        hub: Xena hub URL
        dataset: Dataset name
        output_path: Path to save downloaded data
        
    Returns:
        DataFrame with clinical data
    """
    # Check if file already exists
    if os.path.exists(output_path):
        print(f"Clinical data file already exists at {output_path}. Loading from file...")
        return pd.read_csv(output_path, sep='\t', index_col=0)
    
    print(f"Downloading clinical data from {dataset}...")
    
    try:
        # Get all samples in the dataset
        samples = xena.dataset_samples(hub, dataset, None)
        print(f"Found {len(samples)} samples in clinical dataset")
        
        # Get all fields in the clinical dataset
        fields = xena.dataset_field(hub, dataset)
        print(f"Found {len(fields)} clinical fields")
        
        # Create empty DataFrame with samples as index
        df = pd.DataFrame(index=samples)
        
        # Download data in batches to avoid timeouts
        batch_size = 20
        
        for i in range(0, len(fields), batch_size):
            field_batch = fields[i:i+batch_size]
            print(f"Downloading batch {i//batch_size + 1}/{(len(fields)-1)//batch_size + 1} ({len(field_batch)} clinical fields)")
            
            # Get clinical values for each field
            for field in field_batch:
                try:
                    [_, values] = xena.dataset_probe_values(hub, dataset, samples, [field])
                    df[field] = values[0]
                except Exception as e:
                    print(f"Warning: Could not download field '{field}': {str(e)}")
            
            # Sleep to avoid overwhelming the server
            time.sleep(0.5)
        
        # Save to file
        df.to_csv(output_path, sep='\t')
        print(f"Clinical data saved to {output_path}")
        
        return df
    except Exception as e:
        print(f"Error downloading clinical data: {str(e)}")
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
    
    # Look for sample type column in clinical data
    sample_type_cols = [
        'sample_type', '_sample_type', 'sampleType', 
        'sample_type_id', 'sampleTypeID',
        'SampleType'
    ]
    
    sample_type_col = None
    for col in sample_type_cols:
        if col in clinical_data.columns:
            sample_type_col = col
            break
    
    if sample_type_col:
        print(f"Using column '{sample_type_col}' to identify sample types.")
    else:
        print("Error: Could not find sample type column in clinical data.")
        print("Available columns:", clinical_data.columns.tolist())
        sys.exit(1)
    
    # Create group column based on sample_type_id
    # TCGA sample type codes:
    # 1 = Solid Tissue Normal
    # 0 = Primary Tumor
    clinical_data['Group'] = clinical_data[sample_type_col].apply(
        lambda x: 'Normal' if x == 1 else 'NASH-HCC'
    )
    
    # Print sample counts and sample types for debugging
    print("\nSample type distribution:")
    print(clinical_data[sample_type_col].value_counts())
    print("\nSample group counts:")
    sample_counts = clinical_data['Group'].value_counts()
    print(f"Sample group counts: {dict(sample_counts)}")
    
    # Find common samples
    common_samples = set(mirna_data.index) & set(clinical_data.index)
    print(f"Found {len(common_samples)} common samples between datasets.")
    
    # Convert set to list for indexing
    common_samples_list = list(common_samples)
    
    # Filter datasets to only include common samples
    mirna_filtered = mirna_data.loc[common_samples_list]
    clinical_filtered = clinical_data.loc[common_samples_list]
    
    # Combine data
    combined_data = mirna_filtered.copy()
    combined_data['Group'] = clinical_filtered['Group']
    
    # Reset index to make sample ID a column
    combined_data = combined_data.reset_index()
    combined_data.rename(columns={'index': 'Sample_ID'}, inplace=True)
    
    # Save processed data
    combined_data.to_csv(output_path, index=False)
    print(f"Processed data saved to {output_path}")
    
    return combined_data

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Define output paths
    mirna_path = output_dir / args.mirna_filename
    clinical_path = output_dir / args.clinical_filename
    processed_path = output_dir / args.processed_filename
    
    # Verify the datasets exist
    try:
        mirna_samples = xena.dataset_samples(XENA_HUB, XENA_MIRNA_DATASET, 5)
        clinical_samples = xena.dataset_samples(XENA_HUB, XENA_CLINICAL_DATASET, 5)
        print(f"Successfully connected to Xena Hub. Found dataset samples:")
        print(f"miRNA dataset sample examples: {mirna_samples}")
        print(f"Clinical dataset sample examples: {clinical_samples}")
    except Exception as e:
        print(f"Error connecting to Xena Hub or datasets: {str(e)}")
        sys.exit(1)
    
    # Download expression data
    mirna_data = download_expression_data(XENA_HUB, XENA_MIRNA_DATASET, mirna_path)
    
    # Download clinical data
    clinical_data = download_clinical_data(XENA_HUB, XENA_CLINICAL_DATASET, clinical_path)
    
    # Process data
    process_data(mirna_data, clinical_data, processed_path)
    
    print("\nData acquisition complete. You can now run the target_prioritizer.py script with this data.")
    print(f"Example: python target_prioritizer.py --input {processed_path} --id-col Sample_ID")

if __name__ == "__main__":
    main() 