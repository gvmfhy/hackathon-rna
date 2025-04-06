import pandas as pd
import pyarrow.parquet as pq
import os

# Get the directory containing the script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Construct the path to the Open Targets data directory relative to the script
ot_dir = os.path.join(script_dir, 'data', 'open_targets')

# Get a list of parquet files
parquet_files = [f for f in os.listdir(ot_dir) if f.endswith('.parquet')]

# Select a sample file to examine
sample_file = os.path.join(ot_dir, parquet_files[0])
print(f"Examining file: {sample_file}")

# Read the parquet file
table = pq.read_table(sample_file)
df = table.to_pandas()

# Check if the tractability column exists
if 'tractability' in df.columns:
    print("\nSample tractability data:")
    for index, row in df[['approvedSymbol', 'tractability']].head(10).iterrows():
        print(f"\n{row['approvedSymbol']}:")
        print(f"Type: {type(row['tractability'])}")
        print(f"Value: {row['tractability']}")
else:
    print("No 'tractability' column found")

# Count genes with tractability information
tractable_count = 0
small_molecule_count = 0
antibody_count = 0
other_count = 0

for _, row in df.iterrows():
    if 'tractability' in df.columns and isinstance(row['tractability'], list) and len(row['tractability']) > 0:
        has_tractability = False
        has_small_molecule = False
        has_antibody = False
        has_other = False
        
        for item in row['tractability']:
            if isinstance(item, dict) and 'value' in item and item['value'] == True:
                has_tractability = True
                modality = item.get('modality', '').lower()
                if 'small' in modality:
                    has_small_molecule = True
                elif 'antibody' in modality:
                    has_antibody = True
                else:
                    has_other = True
        
        if has_tractability:
            tractable_count += 1
            if has_small_molecule:
                small_molecule_count += 1
            if has_antibody:
                antibody_count += 1
            if has_other:
                other_count += 1

print(f"\nStats from this sample file:")
print(f"Total genes: {len(df)}")
print(f"Tractable genes: {tractable_count}")
print(f"Small molecule targets: {small_molecule_count}")
print(f"Antibody targets: {antibody_count}")
print(f"Other modality targets: {other_count}")

# Find and display some examples of tractable genes
print("\nExamples of tractable genes:")
for index, row in df.iterrows():
    if 'tractability' in df.columns and isinstance(row['tractability'], list) and len(row['tractability']) > 0:
        for item in row['tractability']:
            if isinstance(item, dict) and 'value' in item and item['value'] == True:
                print(f"\nGene: {row['approvedSymbol']}")
                print(f"Tractability details: {row['tractability']}")
                break 