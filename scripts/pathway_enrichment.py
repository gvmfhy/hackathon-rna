import pandas as pd
import os
import numpy as np
from gprofiler import GProfiler
from datetime import datetime

# Get current directory and set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(current_dir)
results_dir = os.path.join(base_dir, 'results')
data_dir = os.path.join(base_dir, 'data')

# Set file paths
mirna_significant_file = os.path.join(results_dir, 'mirna_de_significant_with_names.csv')
mirtarbase_file = os.path.join(data_dir, 'mirnet', 'miRNet-mir-gene-hsa-mirtarbase.csv')

# Create timestamp for output filenames
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Load miRNA significant results file
print("Loading significantly differentially expressed miRNAs...")
df_mirna = pd.read_csv(mirna_significant_file)

# Verify if file exists
if not os.path.exists(mirna_significant_file):
    print(f"Error: File not found - {mirna_significant_file}")
    exit(1)

# Load the miRTarBase file
print("Loading miRTarBase miRNA-target interactions...")
df_targets = pd.read_csv(mirtarbase_file)

# Verify if file exists
if not os.path.exists(mirtarbase_file):
    print(f"Error: File not found - {mirtarbase_file}")
    exit(1)

# Print basic info about the files
print(f"Number of significant miRNAs: {len(df_mirna)}")
print(f"Number of miRNA-target interactions in miRTarBase: {len(df_targets)}")

# Match significant miRNAs to targets, considering naming differences (miR- vs mir-)
# Create a standardized lowercase version for matching
df_mirna['mirna_lower'] = df_mirna['miRNA_name'].str.lower().str.replace('mir-', 'mir-')
df_targets['mirna_lower'] = df_targets['mir_id'].str.lower().str.replace('mir-', 'mir-')

# Get lists of upregulated and downregulated miRNAs
up_mirnas = df_mirna[df_mirna['log2FC'] > 0]['mirna_lower'].tolist()
down_mirnas = df_mirna[df_mirna['log2FC'] < 0]['mirna_lower'].tolist()

print(f"Upregulated miRNAs: {up_mirnas}")
print(f"Downregulated miRNAs: {down_mirnas}")

# Filter for targets of upregulated miRNAs
up_targets_df = df_targets[df_targets['mirna_lower'].isin(up_mirnas)]
up_target_genes = up_targets_df['symbol'].unique().tolist()

# Filter for targets of downregulated miRNAs
down_targets_df = df_targets[df_targets['mirna_lower'].isin(down_mirnas)]
down_target_genes = down_targets_df['symbol'].unique().tolist()

print(f"Number of target genes for upregulated miRNAs: {len(up_target_genes)}")
print(f"Number of target genes for downregulated miRNAs: {len(down_target_genes)}")

# Save target gene lists to files for reference
with open(os.path.join(results_dir, f'up_mirna_target_genes_{timestamp}.txt'), 'w') as f:
    for gene in sorted(up_target_genes):
        f.write(f"{gene}\n")

with open(os.path.join(results_dir, f'down_mirna_target_genes_{timestamp}.txt'), 'w') as f:
    for gene in sorted(down_target_genes):
        f.write(f"{gene}\n")

# Initialize g:Profiler
gp = GProfiler(return_dataframe=True)

# Function to perform enrichment analysis and save results
def perform_enrichment(gene_list, file_prefix, description):
    if not gene_list:
        print(f"No genes in {description} list, skipping enrichment")
        return None
    
    print(f"Performing pathway enrichment for {description} ({len(gene_list)} genes)...")
    
    # Perform enrichment analysis with g:Profiler
    # Query KEGG, Reactome, and GO Biological Process
    enrichment_results = gp.profile(
        organism='hsapiens',
        query=gene_list,
        sources=['KEGG', 'REAC', 'GO:BP'],
        user_threshold=0.05,  # FDR threshold
        all_results=True,
        no_evidences=False,
        no_iea=True  # Exclude electronically inferred annotations
    )
    
    # Check if we got results
    if enrichment_results is None or len(enrichment_results) == 0:
        print(f"No significant enrichment found for {description}")
        return None
    
    # Save results to CSV
    output_file = os.path.join(results_dir, f'{file_prefix}_pathway_enrichment_{timestamp}.csv')
    enrichment_results.to_csv(output_file, index=False)
    print(f"Saved enrichment results to: {output_file}")
    
    # Print top 10 pathways
    print(f"\nTop 10 enriched pathways for {description}:")
    top_results = enrichment_results.sort_values('p_value').head(10)
    for _, row in top_results.iterrows():
        print(f"  - {row['source']}: {row['name']} (p = {row['p_value']:.2e}, genes: {row['intersection_size']})")
    
    return enrichment_results

# Perform enrichment for targets of upregulated miRNAs
# These genes are potentially suppressed in HCC
up_mirna_enrichment = perform_enrichment(
    up_target_genes, 
    'up_mirna_targets', 
    'targets of upregulated miRNAs (potentially suppressed in HCC)'
)

# Perform enrichment for targets of downregulated miRNAs
# These genes are potentially more expressed in HCC
down_mirna_enrichment = perform_enrichment(
    down_target_genes, 
    'down_mirna_targets', 
    'targets of downregulated miRNAs (potentially overexpressed in HCC)'
)

print("\nPathway enrichment analysis complete!")

# Create a summary file
summary_file = os.path.join(results_dir, f'pathway_analysis_summary_{timestamp}.txt')
with open(summary_file, 'w') as f:
    f.write("# miRNA Target Pathway Enrichment Analysis\n\n")
    
    f.write("## Analysis Overview\n")
    f.write(f"- Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"- Total significant miRNAs: {len(df_mirna)}\n")
    f.write(f"- Upregulated miRNAs: {len(up_mirnas)}\n")
    f.write(f"- Downregulated miRNAs: {len(down_mirnas)}\n")
    f.write(f"- Target genes of upregulated miRNAs: {len(up_target_genes)}\n")
    f.write(f"- Target genes of downregulated miRNAs: {len(down_target_genes)}\n\n")
    
    f.write("## Interpretation\n")
    f.write("- Target genes of upregulated miRNAs: These genes are likely suppressed in HCC compared to normal tissue\n")
    f.write("- Target genes of downregulated miRNAs: These genes are likely more expressed in HCC compared to normal tissue\n\n")
    
    if up_mirna_enrichment is not None and len(up_mirna_enrichment) > 0:
        f.write("## Top Pathways for Targets of Upregulated miRNAs (likely suppressed in HCC)\n")
        top_results = up_mirna_enrichment.sort_values('p_value').head(10)
        for i, (_, row) in enumerate(top_results.iterrows(), 1):
            f.write(f"{i}. {row['source']}: {row['name']} (p = {row['p_value']:.2e}, genes: {row['intersection_size']})\n")
    else:
        f.write("## No significant pathways found for itargets of upregulated miRNAs\n")
    
    f.write("\n")
    
    if down_mirna_enrichment is not None and len(down_mirna_enrichment) > 0:
        f.write("## Top Pathways for Targets of Downregulated miRNAs (likely overexpressed in HCC)\n")
        top_results = down_mirna_enrichment.sort_values('p_value').head(10)
        for i, (_, row) in enumerate(top_results.iterrows(), 1):
            f.write(f"{i}. {row['source']}: {row['name']} (p = {row['p_value']:.2e}, genes: {row['intersection_size']})\n")
    else:
        f.write("## No significant pathways found for targets of downregulated miRNAs\n")

print(f"Summary saved to: {summary_file}")
print("Analysis complete!") 