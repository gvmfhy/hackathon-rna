import pandas as pd
import requests
import os
import json

# Get the current directory and set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(current_dir)
results_dir = os.path.join(base_dir, 'results')
data_dir = os.path.join(base_dir, 'data')

# Input file with significant miRNAs
sig_file = os.path.join(results_dir, 'mirna_de_significant_with_names.csv')

# Output files for gene lists
up_targets_file = os.path.join(data_dir, 'up_mir_targets.txt')
down_targets_file = os.path.join(data_dir, 'down_mir_targets.txt')

# Check if input file exists
if not os.path.exists(sig_file):
    print(f"Error: Significant miRNA file not found at {sig_file}")
    exit(1)

# Load significant miRNAs
df_sig = pd.read_csv(sig_file)

# Separate upregulated and downregulated miRNAs
up_mirnas = df_sig[df_sig['log2FC'] > 0]['miRNA_name'].tolist()
down_mirnas = df_sig[df_sig['log2FC'] < 0]['miRNA_name'].tolist()

print(f"Upregulated miRNAs: {up_mirnas}")
print(f"Downregulated miRNAs: {down_mirnas}")

# miRNet API endpoint
mirnet_api_url = "http://api.mirnet.ca/table/mir"
headers = {
    'Content-Type': 'application/json'
}

# Function to fetch targets from miRNet API
def fetch_mirnet_targets(mirna_list):
    if not mirna_list:
        return set() # Return empty set if list is empty

    # Format miRNA list for API
    query_list = ";".join(mirna_list)
    
    # Updated payload without selSource
    payload_dict = {
        "org": "hsa",
        "idOpt": "mir_id", 
        # "selSource": "Liver", # Removed based on troubleshooting step 1
        "targetOpt": "gene",
        "myList": query_list
    }
    
    # Convert payload to JSON string
    payload_json = json.dumps(payload_dict)
    
    print(f"Querying miRNet for targets of: {query_list}") # Updated print statement
    try:
        # Send request with data parameter as JSON string
        response = requests.post(mirnet_api_url, headers=headers, data=payload_json)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
        
        data = response.json()
        
        # Extract target gene symbols - they seem to be in the 'Target' column
        target_genes = set()
        if isinstance(data, list) and len(data) > 0:
             # Assuming the first dictionary contains the key that holds the list of results
             results_key = list(data[0].keys())[0]
             if results_key in data[0] and isinstance(data[0][results_key], list):
                 for item in data[0][results_key]:
                    if 'Target' in item:
                        target_genes.add(item['Target'])
             else:
                 print(f"Warning: Unexpected JSON structure for key {results_key}. Data: {data[0]}")
        elif isinstance(data, dict) and 'msg' in data:
            print(f"API returned a message (possibly no results): {data['msg']}")
        else:
            print(f"Warning: Received unexpected data format from API: {data}")
            
        print(f"Found {len(target_genes)} unique target genes.")
        return target_genes
        
    except requests.exceptions.RequestException as e:
        print(f"Error querying miRNet API: {e}")
        return set() # Return empty set on error
    except json.JSONDecodeError:
        print("Error decoding JSON response from miRNet API.")
        print(f"Response text: {response.text}") # Print response text for debugging
        return set()
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return set()

# Fetch targets for upregulated miRNAs
up_targets = fetch_mirnet_targets(up_mirnas)

# Fetch targets for downregulated miRNAs
down_targets = fetch_mirnet_targets(down_mirnas)

# Save target gene lists
if up_targets:
    with open(up_targets_file, 'w') as f:
        for gene in sorted(list(up_targets)):
            f.write(f"{gene}\n")
    print(f"Upregulated miRNA target genes saved to: {up_targets_file}")
else:
    print("No target genes found or saved for upregulated miRNAs.")

if down_targets:
    with open(down_targets_file, 'w') as f:
        for gene in sorted(list(down_targets)):
            f.write(f"{gene}\n")
    print(f"Downregulated miRNA target genes saved to: {down_targets_file}")
else:
    print("No target genes found or saved for downregulated miRNAs.")

print("miRNet target fetching complete.") 