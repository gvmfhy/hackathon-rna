#!/usr/bin/env python3
"""
Target Prioritizer for miRNA-Based Therapeutics

This tool identifies and ranks miRNAs differentially expressed between 
NASH-HCC and Normal Liver tissue to select optimal therapeutic targets.

Usage:
    python target_prioritizer.py --input <expression_data.csv> [options]

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Add src directory to path to import our modules
sys.path.append(str(Path(__file__).parent.parent))

# Set plot style
plt.style.use('seaborn-whitegrid')
sns.set_context("talk")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Identify and rank miRNAs for NASH-HCC targeting")
    
    # Input parameters
    parser.add_argument("--input", type=str, required=True,
                      help="Path to miRNA expression data (CSV/TSV)")
    parser.add_argument("--format", type=str, default="csv",
                      help="Input file format: csv or tsv (default: csv)")
    parser.add_argument("--id-col", type=str, default="miRNA",
                      help="Column name containing miRNA IDs")
    parser.add_argument("--group-col", type=str, default="Group",
                      help="Column name for sample groups")
    parser.add_argument("--normal-label", type=str, default="Normal",
                      help="Label for normal liver samples in group column")
    parser.add_argument("--disease-label", type=str, default="NASH-HCC",
                      help="Label for NASH-HCC samples in group column")
    
    # Analysis parameters
    parser.add_argument("--test", type=str, default="t-test", choices=["t-test", "wilcoxon"],
                      help="Statistical test: t-test or wilcoxon (default: t-test)")
    parser.add_argument("--fdr", type=float, default=0.05,
                      help="FDR threshold for significance (default: 0.05)")
    parser.add_argument("--fc-threshold", type=float, default=1.5,
                      help="Fold change threshold (default: 1.5)")
    parser.add_argument("--top-n", type=int, default=10,
                      help="Number of top miRNAs to return (default: 10)")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="results",
                      help="Output directory for results (default: results)")
    parser.add_argument("--prefix", type=str, default="NASH-HCC_miRNA",
                      help="Prefix for output files (default: NASH-HCC_miRNA)")
    
    return parser.parse_args()

def load_data(input_file, file_format='csv', **kwargs):
    """
    Load miRNA expression data from file.
    
    Args:
        input_file: Path to input file
        file_format: File format (csv or tsv)
        **kwargs: Additional arguments for pandas read_csv
        
    Returns:
        DataFrame with miRNA expression data
    """
    try:
        if file_format.lower() == 'tsv':
            data = pd.read_csv(input_file, sep='\t', **kwargs)
        else:
            data = pd.read_csv(input_file, **kwargs)
        
        print(f"Loaded data with {data.shape[0]} miRNAs and {data.shape[1]} samples/columns")
        return data
    except Exception as e:
        print(f"Error loading data: {str(e)}")
        sys.exit(1)

def preprocess_data(data, id_col, group_col, normal_label, disease_label):
    """
    Preprocess miRNA expression data.
    
    Args:
        data: Input DataFrame
        id_col: Column name containing miRNA IDs
        group_col: Column name for sample groups
        normal_label: Label for normal samples
        disease_label: Label for disease samples
        
    Returns:
        Tuple of (processed_data, normal_samples, disease_samples)
    """
    # Check if data has expected columns
    if group_col not in data.columns:
        print(f"Error: Group column '{group_col}' not found in data.")
        print(f"Available columns: {', '.join(data.columns)}")
        sys.exit(1)
    
    # Identify sample groups
    normal_samples = data[data[group_col] == normal_label].drop([id_col, group_col], axis=1, errors='ignore')
    disease_samples = data[data[group_col] == disease_label].drop([id_col, group_col], axis=1, errors='ignore')
    
    if normal_samples.empty or disease_samples.empty:
        print(f"Error: Could not find samples with labels '{normal_label}' and/or '{disease_label}'")
        print(f"Available group labels: {data[group_col].unique()}")
        sys.exit(1)
    
    print(f"Found {normal_samples.shape[0]} normal samples and {disease_samples.shape[0]} disease samples")
    
    # Return processed data
    processed_data = data.set_index(id_col)
    return processed_data, normal_samples, disease_samples

def calculate_differential_expression(data, normal_samples, disease_samples, test_method='t-test'):
    """
    Calculate differential expression statistics.
    
    Args:
        data: Processed data with miRNA IDs as index
        normal_samples: Normal sample expression values
        disease_samples: Disease sample expression values
        test_method: Statistical test method ('t-test' or 'wilcoxon')
        
    Returns:
        DataFrame with differential expression results
    """
    results = pd.DataFrame(index=data.index)
    
    # Calculate log2 fold change
    normal_means = normal_samples.mean(axis=1)
    disease_means = disease_samples.mean(axis=1)
    results['log2FC'] = np.log2(disease_means / normal_means)
    results['FC'] = disease_means / normal_means
    
    # Calculate p-values
    pvals = []
    for idx in data.index:
        normal_expr = normal_samples.loc[idx].values
        disease_expr = disease_samples.loc[idx].values
        
        if test_method == 'wilcoxon':
            # Wilcoxon rank-sum test
            stat, pval = stats.ranksums(disease_expr, normal_expr)
        else:
            # Student's t-test
            stat, pval = stats.ttest_ind(disease_expr, normal_expr, equal_var=False)
        
        pvals.append(pval)
    
    results['pvalue'] = pvals
    
    # Apply FDR correction
    results['padj'] = multipletests(results['pvalue'], method='fdr_bh')[1]
    
    return results

def rank_mirnas(results, fdr_threshold=0.05, fc_threshold=1.5, top_n=10):
    """
    Rank miRNAs based on significance and fold change.
    
    Args:
        results: DataFrame with differential expression results
        fdr_threshold: FDR threshold for significance
        fc_threshold: Fold change threshold
        top_n: Number of top miRNAs to return
        
    Returns:
        DataFrame with ranked miRNAs and their statistics
    """
    # Filter by significance and fold change
    sig_results = results[(results['padj'] < fdr_threshold) & 
                          ((results['FC'] > fc_threshold) | (results['FC'] < 1/fc_threshold))]
    
    # Sort by absolute log2FC
    sig_results['abs_log2FC'] = sig_results['log2FC'].abs()
    ranked_results = sig_results.sort_values('abs_log2FC', ascending=False)
    
    # Add regulation status
    ranked_results['regulation'] = ['Upregulated' if fc > 0 else 'Downregulated' 
                                     for fc in ranked_results['log2FC']]
    
    # Return top N miRNAs
    if len(ranked_results) > top_n:
        return ranked_results.iloc[:top_n]
    return ranked_results

def plot_volcano(results, fdr_threshold=0.05, fc_threshold=1.5, 
                top_n=10, output_dir='results', filename='volcano.png'):
    """
    Generate volcano plot for miRNA differential expression.
    
    Args:
        results: DataFrame with differential expression results
        fdr_threshold: FDR threshold for significance
        fc_threshold: Fold change threshold
        top_n: Number of top miRNAs to label
        output_dir: Output directory
        filename: Output filename
    """
    plt.figure(figsize=(12, 10))
    
    # Define colors for plot
    colors = ['gray'] * len(results)
    
    # Calculate threshold values
    log2_fc_thresh = np.log2(fc_threshold)
    
    # Color significant points
    for i, (idx, row) in enumerate(results.iterrows()):
        if row['padj'] < fdr_threshold:
            if row['log2FC'] > log2_fc_thresh:
                colors[i] = 'red'  # Upregulated
            elif row['log2FC'] < -log2_fc_thresh:
                colors[i] = 'blue'  # Downregulated
    
    # Create scatter plot
    plt.scatter(results['log2FC'], -np.log10(results['pvalue']), c=colors, alpha=0.6)
    
    # Add threshold lines
    plt.axvline(x=log2_fc_thresh, color='gray', linestyle='--')
    plt.axvline(x=-log2_fc_thresh, color='gray', linestyle='--')
    plt.axhline(y=-np.log10(fdr_threshold), color='gray', linestyle='--')
    
    # Label top differentially expressed miRNAs
    sig_results = results[(results['padj'] < fdr_threshold) & 
                          ((results['log2FC'] > log2_fc_thresh) | 
                           (results['log2FC'] < -log2_fc_thresh))]
    
    sig_results = sig_results.sort_values('padj').head(top_n)
    
    for idx, row in sig_results.iterrows():
        plt.annotate(idx, 
                     xy=(row['log2FC'], -np.log10(row['pvalue'])),
                     xytext=(5, 5), textcoords='offset points',
                     fontsize=10, fontweight='bold')
    
    # Set axis labels and title
    plt.xlabel('log2 Fold Change (NASH-HCC / Normal)')
    plt.ylabel('-log10 p-value')
    plt.title('Differential miRNA Expression: NASH-HCC vs Normal Liver')
    
    # Add legend
    plt.plot([], [], 'o', color='red', label='Upregulated')
    plt.plot([], [], 'o', color='blue', label='Downregulated')
    plt.plot([], [], 'o', color='gray', label='Not significant')
    plt.legend()
    
    # Save figure
    output_path = os.path.join(output_dir, filename)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Volcano plot saved to {output_path}")
    plt.close()

def plot_top_mirnas(results, top_n=10, output_dir='results', filename='top_mirnas.png'):
    """
    Generate bar plot for top differentially expressed miRNAs.
    
    Args:
        results: DataFrame with ranked miRNA results
        top_n: Number of top miRNAs to show
        output_dir: Output directory
        filename: Output filename
    """
    if len(results) == 0:
        print("No significant miRNAs found for plotting.")
        return
    
    # Get top N miRNAs
    top_results = results.head(min(top_n, len(results)))
    
    plt.figure(figsize=(12, 8))
    
    # Create color map
    colors = ['red' if fc > 0 else 'blue' for fc in top_results['log2FC']]
    
    # Create bar plot
    bars = plt.barh(top_results.index, top_results['log2FC'], color=colors)
    
    # Add p-value annotations
    for i, (idx, row) in enumerate(top_results.iterrows()):
        plt.text(0, i, f"p-adj = {row['padj']:.2e}", 
                 ha='center', va='center', color='white', fontweight='bold')
    
    # Set axis labels and title
    plt.xlabel('log2 Fold Change (NASH-HCC / Normal)')
    plt.ylabel('miRNA')
    plt.title('Top Differentially Expressed miRNAs: NASH-HCC vs Normal Liver')
    
    # Add legend
    plt.plot([], [], 'o', color='red', label='Upregulated in NASH-HCC')
    plt.plot([], [], 'o', color='blue', label='Downregulated in NASH-HCC')
    plt.legend()
    
    # Save figure
    output_path = os.path.join(output_dir, filename)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Top miRNAs plot saved to {output_path}")
    plt.close()

def save_results(results, output_dir='results', prefix='NASH-HCC_miRNA'):
    """
    Save differential expression results to CSV and JSON.
    
    Args:
        results: DataFrame with differential expression results
        output_dir: Output directory
        prefix: Prefix for output files
    """
    # Save full results to CSV
    csv_path = os.path.join(output_dir, f"{prefix}_differential_expression.csv")
    results.to_csv(csv_path)
    print(f"Full results saved to {csv_path}")
    
    # Save top results to JSON
    json_path = os.path.join(output_dir, f"{prefix}_top_targets.json")
    
    # Get only significant results
    sig_results = results[results['padj'] < 0.05]
    
    # Convert to dictionary for JSON
    result_dict = {
        'analysis_info': {
            'date': pd.Timestamp.now().strftime('%Y-%m-%d'),
            'comparison': 'NASH-HCC vs Normal Liver',
            'total_mirnas_analyzed': len(results),
            'significant_mirnas': len(sig_results)
        },
        'upregulated': {},
        'downregulated': {}
    }
    
    # Add up/down regulated miRNAs
    for idx, row in sig_results.iterrows():
        mirna_data = {
            'log2FC': row['log2FC'],
            'FC': row['FC'],
            'pvalue': row['pvalue'],
            'padj': row['padj']
        }
        
        if row['log2FC'] > 0:
            result_dict['upregulated'][idx] = mirna_data
        else:
            result_dict['downregulated'][idx] = mirna_data
    
    # Save to JSON
    with open(json_path, 'w') as f:
        json.dump(result_dict, f, indent=2)
    
    print(f"Top targets saved to {json_path}")

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Load and preprocess data
    data = load_data(
        args.input, 
        file_format=args.format
    )
    
    processed_data, normal_samples, disease_samples = preprocess_data(
        data,
        id_col=args.id_col,
        group_col=args.group_col,
        normal_label=args.normal_label,
        disease_label=args.disease_label
    )
    
    # Calculate differential expression
    results = calculate_differential_expression(
        processed_data,
        normal_samples,
        disease_samples,
        test_method=args.test
    )
    
    # Rank miRNAs
    ranked_results = rank_mirnas(
        results,
        fdr_threshold=args.fdr,
        fc_threshold=args.fc_threshold,
        top_n=args.top_n
    )
    
    # Generate visualizations
    plot_volcano(
        results,
        fdr_threshold=args.fdr,
        fc_threshold=args.fc_threshold,
        top_n=args.top_n,
        output_dir=args.output_dir,
        filename=f"{args.prefix}_volcano.png"
    )
    
    plot_top_mirnas(
        ranked_results,
        top_n=args.top_n,
        output_dir=args.output_dir,
        filename=f"{args.prefix}_top_mirnas.png"
    )
    
    # Save results
    save_results(
        results,
        output_dir=args.output_dir,
        prefix=args.prefix
    )
    
    # Print summary
    print("\nAnalysis Summary:")
    print(f"Total miRNAs analyzed: {len(results)}")
    print(f"Significant miRNAs (FDR < {args.fdr}): {sum(results['padj'] < args.fdr)}")
    print(f"Upregulated: {sum((results['padj'] < args.fdr) & (results['log2FC'] > 0))}")
    print(f"Downregulated: {sum((results['padj'] < args.fdr) & (results['log2FC'] < 0))}")
    
    if len(ranked_results) > 0:
        print("\nTop differentially expressed miRNAs:")
        for idx, row in ranked_results.iterrows():
            direction = "↑" if row['log2FC'] > 0 else "↓"
            print(f"{idx} {direction} (FC: {row['FC']:.2f}, FDR: {row['padj']:.2e})")
    else:
        print("\nNo significant differentially expressed miRNAs found.")

if __name__ == "__main__":
    main() 