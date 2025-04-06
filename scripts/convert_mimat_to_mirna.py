import pandas as pd
import os
import requests
import time
import json

# Get the current directory and set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(current_dir)
results_dir = os.path.join(base_dir, 'results')

# Find the latest significant results file
def find_latest_file(directory, pattern):
    files = [f for f in os.listdir(directory) if f.startswith(pattern)]
    if not files:
        return None
    return max(files, key=lambda x: os.path.getmtime(os.path.join(directory, x)))

sig_file = find_latest_file(results_dir, 'mirna_de_significant_')
results_file = find_latest_file(results_dir, 'mirna_de_results_')

if not sig_file:
    print("No significant results file found.")
    exit(1)

print(f"Processing file: {sig_file}")
sig_path = os.path.join(results_dir, sig_file)
df = pd.read_csv(sig_path)

# Manual mapping for our top differentially expressed miRNAs
# Based on miRBase and literature searches
mimat_to_mirna = {
    'MIMAT0004552': 'hsa-miR-421',
    'MIMAT0000436': 'hsa-miR-133a-3p',
    'MIMAT0001545': 'hsa-miR-483-5p',
    'MIMAT0002808': 'hsa-miR-511-5p',
    'MIMAT0000765': 'hsa-miR-335-5p',
    'MIMAT0006790': 'hsa-miR-1257',
    'MIMAT0003879': 'hsa-miR-551b-3p',
    'MIMAT0000281': 'hsa-miR-223-3p',
    'MIMAT0005797': 'hsa-miR-1261'
}

# Function to query miRBase REST API
def get_mirna_name(mimat_id):
    if mimat_id in mimat_to_mirna:
        return mimat_to_mirna[mimat_id]

    try:
        # Use miRBase REST API
        url = f"https://mirbase.org/api/v1/accession/{mimat_id}"
        
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            if 'id' in data:
                return data['id']
        
        return f"Unknown ({mimat_id})"
    except Exception as e:
        print(f"Error looking up {mimat_id}: {e}")
        return f"Error ({mimat_id})"

print(f"Converting {len(df)} MIMAT IDs to miRNA names...")

# For each row in the dataframe, add the miRNA name
df['miRNA_name'] = df['miRNA'].apply(get_mirna_name)

# Reorder columns to put miRNA name after MIMAT ID
cols = df.columns.tolist()
cols.insert(1, cols.pop(cols.index('miRNA_name')))
df = df[cols]

# Save the enriched results
output_file = os.path.join(results_dir, 'mirna_de_significant_with_names.csv')
df.to_csv(output_file, index=False)
print(f"Results with miRNA names saved to: {output_file}")

# Create a mapping file for future reference
mapping_df = pd.DataFrame({
    'MIMAT_ID': list(mimat_to_mirna.keys()),
    'miRNA_name': list(mimat_to_mirna.values())
})
mapping_file = os.path.join(base_dir, 'data', 'mimat_to_mirna_mapping.csv')
mapping_df.to_csv(mapping_file, index=False)
print(f"MIMAT to miRNA name mapping saved to: {mapping_file}")

# If we also have the full results file, create an enriched version of that too
if results_file:
    print(f"Also processing full results file: {results_file}")
    results_path = os.path.join(results_dir, results_file)
    full_df = pd.read_csv(results_path)
    
    # Add miRNA names to all results (only for the top 50 by p-value to keep it manageable)
    top_50 = full_df.sort_values('pvalue').head(50)
    
    for mimat_id in top_50['miRNA'].unique():
        if mimat_id not in mimat_to_mirna:
            mirna_name = get_mirna_name(mimat_id)
            mimat_to_mirna[mimat_id] = mirna_name
    
    # Update the mapping file
    mapping_df = pd.DataFrame({
        'MIMAT_ID': list(mimat_to_mirna.keys()),
        'miRNA_name': list(mimat_to_mirna.values())
    })
    mapping_df.to_csv(mapping_file, index=False)
    
    # Add miRNA names to the full results where we have them
    full_df['miRNA_name'] = full_df['miRNA'].map(mimat_to_mirna).fillna("Not mapped")
    
    # Create a smaller version with just the top 50 miRNAs by p-value
    top_50_df = full_df.sort_values('pvalue').head(50).copy()
    
    # Reorder columns
    cols = top_50_df.columns.tolist()
    cols.insert(1, cols.pop(cols.index('miRNA_name')))
    top_50_df = top_50_df[cols]
    
    # Save the top 50 results
    output_file = os.path.join(results_dir, 'mirna_de_top50_with_names.csv')
    top_50_df.to_csv(output_file, index=False)
    print(f"Top 50 results with miRNA names saved to: {output_file}")

print("Conversion complete!") 