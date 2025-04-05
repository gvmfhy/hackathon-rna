#!/usr/bin/env python3
"""
GEO miRNA Data Acquisition for NASH-HCC vs Normal Liver Comparison

This script downloads and processes miRNA expression data from the Gene Expression Omnibus (GEO)
for NASH-HCC vs Normal Liver comparison studies.

Usage:
    python fetch_geo_data.py --geo-id GSE12345 --output-dir data

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import GEOparse
import gzip
import shutil
import re
import tempfile

# Add src directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Download miRNA data from GEO for NASH-HCC vs Normal Liver comparison")
    
    # Input parameters
    parser.add_argument("--geo-id", type=str, required=True,
                      help="GEO dataset ID (e.g., GSE12345)")
    parser.add_argument("--platform-id", type=str, default=None,
                      help="GEO platform ID (e.g., GPL21263). If not provided, will use the first platform.")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="data",
                      help="Output directory for downloaded data (default: data)")
    parser.add_argument("--processed-filename", type=str, default=None,
                      help="Filename for processed data (default: {geo_id}_processed.csv)")
    
    # Processing parameters
    parser.add_argument("--normal-keywords", type=str, default="normal,control,healthy",
                      help="Comma-separated keywords to identify normal samples (default: normal,control,healthy)")
    parser.add_argument("--disease-keywords", type=str, default="NASH-HCC,HCC,nash,hepatocellular carcinoma",
                      help="Comma-separated keywords to identify disease samples (default: NASH-HCC,HCC,nash,hepatocellular carcinoma)")
    
    # Format parameters
    parser.add_argument("--id-col", type=str, default="ID_REF",
                      help="Column name for miRNA IDs (default: ID_REF)")
    
    return parser.parse_args()

def download_geo_dataset(geo_id, output_dir):
    """
    Download GEO dataset.
    
    Args:
        geo_id: GEO dataset ID (e.g., GSE12345)
        output_dir: Output directory
        
    Returns:
        GEOparse.GSE object
    """
    print(f"Downloading GEO dataset: {geo_id}...")
    
    # Create temp directory for GEO downloads
    tmp_dir = Path(tempfile.mkdtemp())
    
    try:
        # Download GEO dataset
        gse = GEOparse.get_GEO(geo=geo_id, destdir=str(tmp_dir), silent=False)
        
        # Save raw data
        raw_dir = Path(output_dir) / "raw"
        raw_dir.mkdir(exist_ok=True, parents=True)
        
        # Get soft file path from GEO
        soft_file = list(tmp_dir.glob(f"{geo_id}_*.soft.gz"))[0]
        
        # Copy to output directory
        shutil.copy(soft_file, raw_dir / soft_file.name)
        print(f"Raw data saved to {raw_dir / soft_file.name}")
        
        return gse
    except Exception as e:
        print(f"Error downloading GEO dataset: {str(e)}")
        sys.exit(1)
    finally:
        # Clean up temp files
        shutil.rmtree(tmp_dir)

def extract_sample_info(gse, normal_keywords, disease_keywords):
    """
    Extract sample information from GEO dataset.
    
    Args:
        gse: GEOparse.GSE object
        normal_keywords: List of keywords to identify normal samples
        disease_keywords: List of keywords to identify disease samples
        
    Returns:
        DataFrame with sample information
    """
    print("Extracting sample information...")
    
    # Get sample information
    samples = []
    
    for gsm_name, gsm in gse.gsms.items():
        sample_info = {
            'sample_id': gsm_name,
            'title': gsm.metadata.get('title', [''])[0],
            'source': gsm.metadata.get('source_name_ch1', [''])[0],
            'characteristics': ' '.join([' '.join(char) for char in gsm.metadata.get('characteristics_ch1', [])]),
            'description': ' '.join(gsm.metadata.get('description', [''])),
        }
        
        # Combine all text fields for keyword matching
        all_text = ' '.join([
            sample_info['title'].lower(),
            sample_info['source'].lower(),
            sample_info['characteristics'].lower(),
            sample_info['description'].lower()
        ])
        
        # Determine group based on keywords
        if any(kw.lower() in all_text for kw in normal_keywords):
            sample_info['group'] = 'Normal'
        elif any(kw.lower() in all_text for kw in disease_keywords):
            sample_info['group'] = 'NASH-HCC'
        else:
            sample_info['group'] = 'Unknown'
        
        samples.append(sample_info)
    
    # Convert to DataFrame
    samples_df = pd.DataFrame(samples)
    
    # Print group counts
    group_counts = samples_df['group'].value_counts()
    print(f"Sample groups detected: {dict(group_counts)}")
    
    if group_counts.get('Unknown', 0) > 0:
        print("Warning: Some samples could not be classified as Normal or NASH-HCC.")
        print("You may need to manually review and classify these samples.")
    
    return samples_df

def extract_expression_data(gse, gpl_id, samples_df, id_col='ID_REF'):
    """
    Extract expression data from GEO dataset.
    
    Args:
        gse: GEOparse.GSE object
        gpl_id: GEO platform ID
        samples_df: DataFrame with sample information
        id_col: Column name for miRNA IDs
        
    Returns:
        DataFrame with expression data
    """
    print(f"Extracting expression data for platform: {gpl_id}...")
    
    # Get all samples for this platform
    platform_samples = gse.pivot_samples(gpl_id)
    
    # Check if output is empty
    if platform_samples.empty:
        print(f"Error: No expression data found for platform {gpl_id}")
        sys.exit(1)
    
    # Rename index to ID_REF if needed
    if platform_samples.index.name != id_col:
        platform_samples.index.name = id_col
    
    return platform_samples

def process_geo_data(expression_data, samples_df, id_col='ID_REF'):
    """
    Process GEO expression data for differential expression analysis.
    
    Args:
        expression_data: DataFrame with expression data
        samples_df: DataFrame with sample information
        id_col: Column name for miRNA IDs
        
    Returns:
        DataFrame with processed data
    """
    print("Processing expression data...")
    
    # Filter samples to only include Normal and NASH-HCC
    valid_samples = samples_df[samples_df['group'].isin(['Normal', 'NASH-HCC'])]
    
    if len(valid_samples) == 0:
        print("Error: No valid samples found after filtering")
        sys.exit(1)
    
    # Select only relevant sample columns
    expression_subset = expression_data[valid_samples['sample_id']]
    
    # Reset index to get ID_REF as a column
    expression_subset = expression_subset.reset_index()
    
    # Melt data to long format
    melted_data = pd.melt(
        expression_subset, 
        id_vars=[id_col],
        value_vars=valid_samples['sample_id'],
        var_name='sample_id',
        value_name='expression'
    )
    
    # Merge with sample information
    merged_data = pd.merge(
        melted_data,
        valid_samples[['sample_id', 'group']],
        on='sample_id',
        how='left'
    )
    
    # Convert expression to numeric if needed
    merged_data['expression'] = pd.to_numeric(merged_data['expression'], errors='coerce')
    
    # Reshape to wide format with miRNAs as rows and samples as columns
    wide_data = merged_data.pivot_table(
        index=id_col,
        columns='sample_id',
        values='expression'
    ).reset_index()
    
    # Add group column
    wide_data['Group'] = pd.Series(dtype='object')
    
    # Convert to final format expected by target_prioritizer.py
    for i, row in samples_df.iterrows():
        if row['sample_id'] in wide_data.columns:
            wide_data.loc[0, 'Group'] = row['group']
    
    return wide_data

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Parse keywords
    normal_keywords = args.normal_keywords.split(',')
    disease_keywords = args.disease_keywords.split(',')
    
    # Download GEO dataset
    gse = download_geo_dataset(args.geo_id, output_dir)
    
    # Extract sample information
    samples_df = extract_sample_info(gse, normal_keywords, disease_keywords)
    
    # Determine platform ID if not provided
    if args.platform_id is None:
        if len(gse.gpls) == 1:
            platform_id = list(gse.gpls.keys())[0]
        else:
            print(f"Multiple platforms found: {', '.join(gse.gpls.keys())}")
            print("Please specify a platform using --platform-id")
            sys.exit(1)
    else:
        platform_id = args.platform_id
    
    # Extract expression data
    expression_data = extract_expression_data(gse, platform_id, samples_df, args.id_col)
    
    # Process data
    processed_data = process_geo_data(expression_data, samples_df, args.id_col)
    
    # Save processed data
    if args.processed_filename is None:
        processed_filename = f"{args.geo_id}_processed.csv"
    else:
        processed_filename = args.processed_filename
    
    processed_output_path = output_dir / processed_filename
    processed_data.to_csv(processed_output_path, index=False)
    
    print(f"Processed data saved to {processed_output_path}")
    print("\nData acquisition complete. You can now run the target_prioritizer.py script with this data.")
    print(f"Example: python target_prioritizer.py --input {processed_output_path} --id-col {args.id_col}")

if __name__ == "__main__":
    main() 