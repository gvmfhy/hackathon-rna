#!/usr/bin/env python3
"""
NASH-HCC vs Normal Liver miRNA Expression Analysis Workflow

This script provides an end-to-end workflow for acquiring and analyzing miRNA expression data
for NASH-HCC vs Normal Liver comparison. It supports multiple data sources and analysis methods.

Usage:
    python run_nash_hcc_analysis.py --source xena --analysis python

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import time
import shutil

# Add src directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run end-to-end NASH-HCC vs Normal Liver miRNA expression analysis")
    
    # Data source parameters
    parser.add_argument("--source", type=str, required=True, choices=["xena", "geo", "custom"],
                      help="Data source: xena (TCGA via UCSC Xena), geo (Gene Expression Omnibus), or custom (local data)")
    parser.add_argument("--geo-id", type=str, default=None,
                      help="GEO dataset ID (e.g., GSE12345) when source is 'geo'")
    parser.add_argument("--custom-file", type=str, default=None,
                      help="Custom input file path when source is 'custom'")
    
    # Analysis parameters
    parser.add_argument("--analysis", type=str, required=True, choices=["python", "r"],
                      help="Analysis method: python (t-test/Wilcoxon) or r (DESeq2)")
    parser.add_argument("--test", type=str, default="t-test", choices=["t-test", "wilcoxon"],
                      help="Statistical test for Python analysis: t-test or wilcoxon (default: t-test)")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="results",
                      help="Output directory for results (default: results)")
    parser.add_argument("--data-dir", type=str, default="data",
                      help="Directory for intermediate data files (default: data)")
    parser.add_argument("--prefix", type=str, default="NASH-HCC_miRNA",
                      help="Prefix for output files (default: NASH-HCC_miRNA)")
    
    # Filtering parameters
    parser.add_argument("--fdr", type=float, default=0.05,
                      help="FDR threshold for significance (default: 0.05)")
    parser.add_argument("--fc-threshold", type=float, default=1.5,
                      help="Fold change threshold (default: 1.5)")
    
    return parser.parse_args()

def run_xena_acquisition(data_dir, timestamp):
    """
    Run Xena data acquisition.
    
    Args:
        data_dir: Data directory
        timestamp: Timestamp for output files
    
    Returns:
        Path to processed data file
    """
    print("\n=== Acquiring data from UCSC Xena Browser (TCGA LIHC) ===\n")
    
    cmd = [
        "python", 
        str(Path(__file__).parent / "fetch_xena_data.py"),
        "--output-dir", str(data_dir),
        "--processed-filename", f"TCGA_LIHC_miRNA_processed_{timestamp}.csv"
    ]
    
    subprocess.run(cmd, check=True)
    
    return data_dir / f"TCGA_LIHC_miRNA_processed_{timestamp}.csv"

def run_geo_acquisition(geo_id, data_dir, timestamp):
    """
    Run GEO data acquisition.
    
    Args:
        geo_id: GEO dataset ID
        data_dir: Data directory
        timestamp: Timestamp for output files
    
    Returns:
        Path to processed data file
    """
    print(f"\n=== Acquiring data from GEO (dataset: {geo_id}) ===\n")
    
    if not geo_id:
        print("Error: GEO ID is required when source is 'geo'")
        sys.exit(1)
    
    cmd = [
        "python", 
        str(Path(__file__).parent / "fetch_geo_data.py"),
        "--geo-id", geo_id,
        "--output-dir", str(data_dir),
        "--processed-filename", f"{geo_id}_processed_{timestamp}.csv"
    ]
    
    subprocess.run(cmd, check=True)
    
    return data_dir / f"{geo_id}_processed_{timestamp}.csv"

def run_python_analysis(input_file, output_dir, prefix, test_method, fdr, fc_threshold):
    """
    Run Python-based differential expression analysis.
    
    Args:
        input_file: Input data file
        output_dir: Output directory
        prefix: Prefix for output files
        test_method: Statistical test method
        fdr: FDR threshold
        fc_threshold: Fold change threshold
    """
    print("\n=== Running Python-based differential expression analysis ===\n")
    
    cmd = [
        "python",
        str(Path(__file__).parent.parent / "target_prioritizer.py"),
        "--input", str(input_file),
        "--output-dir", str(output_dir),
        "--prefix", prefix,
        "--test", test_method,
        "--fdr", str(fdr),
        "--fc-threshold", str(fc_threshold)
    ]
    
    subprocess.run(cmd, check=True)

def run_r_analysis(input_file, output_dir, timestamp, fdr, log2fc_threshold):
    """
    Run R-based differential expression analysis (DESeq2).
    
    Args:
        input_file: Input data file
        output_dir: Output directory
        timestamp: Timestamp for output files
        fdr: FDR threshold
        log2fc_threshold: Log2 fold change threshold
    """
    print("\n=== Running R-based differential expression analysis (DESeq2) ===\n")
    
    plots_dir = output_dir / "plots"
    output_file = output_dir / f"deseq2_results_{timestamp}.csv"
    
    cmd = [
        "Rscript",
        str(Path(__file__).parent / "deseq2_analysis.R"),
        "--input", str(input_file),
        "--output", str(output_file),
        "--plots-dir", str(plots_dir),
        "--padj-cutoff", str(fdr),
        "--log2fc-cutoff", str(log2fc_threshold)
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print("Error: R analysis failed. Make sure R and required packages are installed.")
        print("Install R packages with:")
        print('  R -e "if (!require(\'BiocManager\', quietly = TRUE)) install.packages(\'BiocManager\'); BiocManager::install(c(\'DESeq2\', \'dplyr\', \'readr\', \'argparse\', \'ggplot2\', \'pheatmap\'))"')
        sys.exit(1)
    except FileNotFoundError:
        print("Error: Rscript not found. Make sure R is installed and in your PATH.")
        sys.exit(1)

def main():
    """Main function."""
    args = parse_arguments()
    
    # Generate timestamp for file naming
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    
    # Create directories
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    
    data_dir.mkdir(exist_ok=True, parents=True)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Acquire data
    if args.source == "xena":
        input_file = run_xena_acquisition(data_dir, timestamp)
    elif args.source == "geo":
        input_file = run_geo_acquisition(args.geo_id, data_dir, timestamp)
    elif args.source == "custom":
        if not args.custom_file:
            print("Error: Custom file path is required when source is 'custom'")
            sys.exit(1)
        input_file = Path(args.custom_file)
        if not input_file.exists():
            print(f"Error: Custom file not found: {input_file}")
            sys.exit(1)
    
    # Run analysis
    analysis_prefix = f"{args.prefix}_{timestamp}"
    
    if args.analysis == "python":
        run_python_analysis(
            input_file, 
            output_dir, 
            analysis_prefix, 
            args.test, 
            args.fdr, 
            args.fc_threshold
        )
    elif args.analysis == "r":
        # For R analysis, calculate log2 fold change threshold
        log2fc_threshold = abs(float(args.fc_threshold))
        
        run_r_analysis(
            input_file, 
            output_dir, 
            timestamp, 
            args.fdr, 
            log2fc_threshold
        )
    
    print("\n=== Analysis Complete ===\n")
    print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main() 