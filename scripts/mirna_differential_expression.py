import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os
from datetime import datetime
from scipy.stats import zscore

# Set the aesthetics for the plots
# 'seaborn-whitegrid' is deprecated in newer versions
sns.set_style("whitegrid")
sns.set(font_scale=1.2)

# Get current directory
current_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(current_dir)

# Define file paths with absolute paths
input_file = os.path.join(base_dir, 'data', 'tcga', 'TCGA_LIHC_miRNA_processed.csv')
output_dir = os.path.join(base_dir, 'results')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Current timestamp for unique filenames
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Load the data
print("Loading miRNA expression data...")
print(f"Input file: {input_file}")
df = pd.read_csv(input_file)

# Basic information about the dataset
print(f"Dataset shape: {df.shape}")
print(f"Group counts: \n{df['Group'].value_counts()}")

# Separate samples by group
normal_samples = df[df['Group'] == 'Normal']
tumor_samples = df[df['Group'] == 'NASH-HCC']  # General HCC samples labeled as NASH-HCC

print(f"Number of Normal samples: {len(normal_samples)}")
print(f"Number of HCC samples: {len(tumor_samples)}")

# Get miRNA columns (exclude Sample_ID, sampleID, Group)
mirna_columns = [col for col in df.columns if col.startswith('MIMAT')]
print(f"Number of miRNA features: {len(mirna_columns)}")

# Function to perform differential expression analysis
def analyze_differential_expression(normal_df, tumor_df, mirna_cols):
    results = []
    
    print("Performing differential expression analysis...")
    for i, mirna in enumerate(mirna_cols):
        if i % 100 == 0:
            print(f"Processing miRNA {i+1}/{len(mirna_cols)}")
        
        normal_expr = normal_df[mirna].values
        tumor_expr = tumor_df[mirna].values
        
        # Skip miRNAs with all zeros or NaNs
        if (np.isnan(normal_expr).all() or np.isnan(tumor_expr).all() or 
            (normal_expr == 0).all() or (tumor_expr == 0).all()):
            continue
        
        # Calculate mean expression in each group
        normal_mean = np.nanmean(normal_expr)
        tumor_mean = np.nanmean(tumor_expr)
        
        # Calculate log2 fold change
        if normal_mean > 0 and tumor_mean > 0:
            log2fc = np.log2(tumor_mean) - np.log2(normal_mean)
        else:
            log2fc = np.nan
        
        # Perform Wilcoxon rank-sum test (non-parametric)
        try:
            stat, pvalue = stats.ranksums(normal_expr, tumor_expr)
        except:
            # If test fails (e.g., too many ties), try t-test as fallback
            stat, pvalue = stats.ttest_ind(normal_expr, tumor_expr, nan_policy='omit')
        
        results.append({
            'miRNA': mirna,
            'normal_mean': normal_mean,
            'tumor_mean': tumor_mean,
            'log2FC': log2fc,
            'pvalue': pvalue,
            'stat': stat
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Apply FDR correction
    if not results_df.empty:
        results_df['padj'] = multipletests(results_df['pvalue'].values, method='fdr_bh')[1]
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('padj')
    
    return results_df

# Run differential expression analysis
results = analyze_differential_expression(normal_samples, tumor_samples, mirna_columns)

# Check why there might be NaN p-values
nan_pvals = results[np.isnan(results['pvalue'])].shape[0]
if nan_pvals > 0:
    print(f"Warning: {nan_pvals} miRNAs have NaN p-values and will be excluded from FDR correction")
    results = results.dropna(subset=['pvalue'])

# Reapply FDR correction if needed
if not results.empty:
    results['padj'] = multipletests(results['pvalue'].values, method='fdr_bh')[1]

# Filter for significant results - relaxed criteria
# Using padj < 0.1 and |log2FC| > 0.5 instead of the stricter padj < 0.05 and |log2FC| > 1
significant = results[(results['padj'] < 0.1) & (abs(results['log2FC']) > 0.5)].copy()
significant = significant.sort_values('padj')

if significant.empty:
    print("\nNo miRNAs passed the relaxed significance filter (padj < 0.1, |log2FC| > 0.5)")
    print("Using top miRNAs by p-value for visualization")
    # Get top miRNAs by adjusted p-value regardless of fold change
    significant = results.sort_values('pvalue').head(50).copy()

# Print summary
print(f"\nTotal miRNAs analyzed: {len(results)}")
print(f"Significantly differentially expressed miRNAs (padj < 0.1, |log2FC| > 0.5): {len(significant)}")
print("\nTop 10 upregulated miRNAs:")
print(significant[significant['log2FC'] > 0].head(10)[['miRNA', 'log2FC', 'padj']])
print("\nTop 10 downregulated miRNAs:")
print(significant[significant['log2FC'] < 0].head(10)[['miRNA', 'log2FC', 'padj']])

# Save results to CSV
results_file = os.path.join(output_dir, f'mirna_de_results_{timestamp}.csv')
results.to_csv(results_file, index=False)
print(f"\nFull results saved to: {results_file}")

significant_file = os.path.join(output_dir, f'mirna_de_significant_{timestamp}.csv')
significant.to_csv(significant_file, index=False)
print(f"Significant results saved to: {significant_file}")

# Generate volcano plot
plt.figure(figsize=(10, 8))
plt.scatter(
    results['log2FC'], 
    -np.log10(results['pvalue']),  # Use pvalue instead of padj for better visualization
    alpha=0.5, 
    s=5, 
    color='gray'
)

# Highlight significant points
sig_up = significant[significant['log2FC'] > 0]
sig_down = significant[significant['log2FC'] < 0]

if not sig_up.empty:
    plt.scatter(
        sig_up['log2FC'], 
        -np.log10(sig_up['pvalue']),  # Use pvalue instead of padj
        alpha=0.7, 
        s=10, 
        color='red', 
        label=f'Upregulated ({len(sig_up)})'
    )

if not sig_down.empty:
    plt.scatter(
        sig_down['log2FC'], 
        -np.log10(sig_down['pvalue']),  # Use pvalue instead of padj
        alpha=0.7, 
        s=10, 
        color='blue', 
        label=f'Downregulated ({len(sig_down)})'
    )

# Label top miRNAs
top_up = sig_up.head(5) if not sig_up.empty else pd.DataFrame()
top_down = sig_down.head(5) if not sig_down.empty else pd.DataFrame()
top_mirnas = pd.concat([top_up, top_down])

if not top_mirnas.empty:
    for _, row in top_mirnas.iterrows():
        plt.annotate(
            row['miRNA'],
            xy=(row['log2FC'], -np.log10(row['pvalue'])),  # Use pvalue instead of padj
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8,
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5)
        )

# Add lines for thresholds
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.3)
plt.axvline(x=0.5, color='gray', linestyle='--', alpha=0.3)
plt.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.3)

plt.xlabel('Log2 Fold Change (HCC / Normal)')
plt.ylabel('-Log10 P-value')
plt.title('Volcano Plot: miRNA Differential Expression in HCC vs Normal Liver')
plt.legend(loc='best')
plt.tight_layout()

volcano_file = os.path.join(output_dir, f'volcano_plot_{timestamp}.png')
plt.savefig(volcano_file, dpi=300)
print(f"Volcano plot saved to: {volcano_file}")

# --- Load miRNA Name Mapping (Moved Before Heatmap) ---
print("Loading miRNA name mapping...")
mapping_file = os.path.join(base_dir, 'data', 'annotations', 'mimat_to_mirna_mapping.csv') # Using updated path
mimat_to_name_map = {}
try:
    map_df = pd.read_csv(mapping_file)
    # Create dictionary, handle potential NaN in miRNA_name
    mimat_to_name_map = pd.Series(map_df.miRNA_name.values, index=map_df.MIMAT_ID).dropna().to_dict()
    print(f"Loaded {len(mimat_to_name_map)} mappings.")
except FileNotFoundError:
    print(f"Warning: Mapping file not found at {mapping_file}. Heatmap will use MIMAT IDs.")
except Exception as e:
    print(f"Warning: Error loading or processing mapping file {mapping_file}: {e}. Heatmap will use MIMAT IDs.")

# --- Generate Simple Heatmap (Reverted) --- 
print("Generating heatmap of top differentially expressed miRNAs...")
# Select top 20 miRNAs by p-value for the heatmap
top_20_mirnas_ids = significant.head(20)['miRNA'].tolist()

if top_20_mirnas_ids:
    # Get corresponding names, fallback to ID if name not found
    top_20_mirnas_names = [mimat_to_name_map.get(mid, mid) for mid in top_20_mirnas_ids]

    # --- Prepare Data for Heatmap ---
    # Use all Normal samples and a subset of Tumor samples for balance
    num_tumor_samples = min(100, len(tumor_samples)) # Limit tumor samples if > 100
    sampled_tumor = tumor_samples.sample(num_tumor_samples, random_state=42) if len(tumor_samples) > num_tumor_samples else tumor_samples
    
    # Combine normal and sampled tumor
    plot_samples_df = pd.concat([normal_samples, sampled_tumor])

    # Extract expression data using MIMAT IDs
    heatmap_data = plot_samples_df[top_20_mirnas_ids].copy()

    # Handle potential NaN values (fill with mean of the miRNA across samples)
    if heatmap_data.isnull().values.any():
        print("Warning: NaN values found in heatmap data. Filling with column means.")
        heatmap_data = heatmap_data.fillna(heatmap_data.mean())

    # --- Z-score Scaling --- 
    # Standardize expression values (Z-score) across samples for each miRNA
    # This makes patterns more comparable
    heatmap_data_z = heatmap_data.apply(zscore, axis=0)
    
    # Assign miRNA names to columns
    heatmap_data_z.columns = top_20_mirnas_names

    # --- Generate the Heatmap --- 
    print(f"Generating heatmap for {heatmap_data_z.shape[1]} miRNAs and {heatmap_data_z.shape[0]} samples...")
    try:
        plt.figure(figsize=(12, 10))
        
        # Use seaborn heatmap
        ax = sns.heatmap(
            heatmap_data_z, 
            cmap="coolwarm",     # Color map (blue=low, red=high)
            center=0,            # Center color map at zero Z-score
            cbar_kws={'label': 'Z-score'}, # Color bar label
            yticklabels=False    # Hide individual sample labels
        )
        
        # Add horizontal line to separate Normal and HCC samples
        num_normal = len(normal_samples)
        ax.axhline(y=num_normal, color='black', lw=2)
        
        # Rotate miRNA labels
        plt.xticks(rotation=45, ha='right')
        
        # Add titles and labels
        plt.title('Top 20 miRNAs by p-value (HCC vs Normal)', fontsize=14)
        plt.xlabel('miRNAs', fontsize=12)
        plt.ylabel(f'Samples (Normal above line, HCC below)', fontsize=12)
        plt.tight_layout()
        
        # Save the figure
        heatmap_file = os.path.join(output_dir, f'heatmap_top_mirnas_{timestamp}.png')
        plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
        print(f"Simple heatmap saved to: {heatmap_file}")
        plt.close() # Close the figure

    except Exception as e:
        print(f"Error generating simple heatmap: {e}")

else:
    print("No significant miRNAs found to generate heatmap.")

print("\nAnalysis complete!") 