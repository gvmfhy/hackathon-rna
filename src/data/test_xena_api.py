#!/usr/bin/env python3
"""
Test script to list datasets on the Xena Hub and find relevant LIHC miRNA datasets.
"""

import pandas as pd
import re

# Import our SSL fix first
try:
    import xena_ssl_fix
except ImportError:
    print("Warning: xena_ssl_fix module not found. SSL certificate errors may occur.")

# Import xenaPython after the SSL fix
try:
    import xenaPython as xena
    # Print available functions/attributes to debug
    print("\nAvailable attributes in xenaPython:")
    print(dir(xena))
    print("-"*30 + "\n")
except ImportError:
    print("Error: xenaPython library not found. Please install it.")
    exit(1)

# Constants
XENA_HUB = "https://tcga.xenahubs.net"

def find_lihc_mirna_datasets(hub):
    """
    Connects to the Xena hub, lists all datasets using available functions,
    and filters for LIHC miRNA data.
    
    Args:
        hub: Xena hub URL
        
    Returns:
        List of potential matching dataset names.
    """
    print(f"Connecting to Xena Hub: {hub}")
    
    try:
        # Step 1: Get all cohorts
        print("Fetching all cohorts...")
        all_cohorts_list = xena.all_cohorts(hub, [])
        if not all_cohorts_list:
            print("Error: Could not retrieve cohort list.")
            return None
        print(f"Found {len(all_cohorts_list)} cohorts.")

        # Step 2: Get dataset metadata for all cohorts
        print("Fetching dataset list for all cohorts... (This might take a while)")
        all_datasets_metadata = xena.dataset_list(hub, all_cohorts_list)
        if not all_datasets_metadata:
             print("Error: Could not retrieve dataset list from cohorts.")
             return None
        print(f"Found metadata for {len(all_datasets_metadata)} datasets.")

        # Extract dataset names from metadata
        all_dataset_names = [ds['name'] for ds in all_datasets_metadata if 'name' in ds]
        print(f"Extracted {len(all_dataset_names)} dataset names.")

        # Step 3: Filter dataset names for LIHC miRNA
        lihc_mirna_datasets = []
        print("\nFiltering for LIHC miRNA datasets...")
        pattern = re.compile(r'TCGA.*LIHC.*(miRNA|HiSeq|miR)', re.IGNORECASE)
        
        for dataset_name in all_dataset_names:
            if pattern.search(dataset_name):
                lihc_mirna_datasets.append(dataset_name)
                
        return lihc_mirna_datasets
        
    except Exception as e:
        print(f"\nError interacting with Xena Hub: {str(e)}")
        # Print available attributes again if there's an error during execution
        try:
            print("\nAvailable attributes in xenaPython:")
            print(dir(xena))
        except:
            pass
        return None

if __name__ == "__main__":
    matching_datasets = find_lihc_mirna_datasets(XENA_HUB)
    
    if matching_datasets is not None:
        if matching_datasets:
            print("\nFound potential LIHC miRNA datasets:")
            for ds in matching_datasets:
                print(f"- {ds}")
        else:
            print("\nNo datasets matching 'LIHC' and 'miRNA/HiSeq/miR' found.")
    else:
        print("\nCould not retrieve dataset list.") 