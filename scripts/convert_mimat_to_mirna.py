import pandas as pd
import os
import argparse
import requests # Re-add requests for API fallback
import time     # For potential delays in API calls

# --- Function to query miRBase REST API (Fallback) ---
def get_mirna_name_from_api(mimat_id):
    """Queries miRBase API to get miRNA name for a MIMAT ID."""
    print(f"  Attempting API lookup for {mimat_id}...")
    try:
        # Use miRBase REST API
        url = f"https://mirbase.org/api/v1/accession/{mimat_id}"
        # Add a small delay to avoid overwhelming the API
        time.sleep(0.2)
        
        response = requests.get(url, timeout=10) # Added timeout
        response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
        
        data = response.json()
        if data and 'id' in data:
            print(f"    API Success: Found name '{data['id']}'")
            return data['id']
        else:
            print(f"    API Warning: No 'id' field in response for {mimat_id}")
            return None # Indicate failure to find name via API
            
    except requests.exceptions.RequestException as e:
        print(f"    API Error looking up {mimat_id}: {e}")
        return None # Indicate failure
    except Exception as e:
        print(f"    Unexpected Error during API lookup for {mimat_id}: {e}")
        return None # Indicate failure

def main():
    parser = argparse.ArgumentParser(description="Convert MIMAT IDs to miRNA names in a results file.")
    parser.add_argument("--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("--output", required=True, help="Path for the output CSV file with names.")
    parser.add_argument("--map", required=True, help="Path to the MIMAT ID to miRNA name mapping CSV file.")
    parser.add_argument("--col", required=True, help="Name of the column containing MIMAT IDs in the input file.")
    args = parser.parse_args()

    # --- Load Input File ---
    print(f"Loading input file: {args.input}")
    try:
        df = pd.read_csv(args.input)
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input}")
        exit(1)
    except Exception as e:
        print(f"Error reading input file {args.input}: {e}")
        exit(1)

    # Check if specified column exists
    if args.col not in df.columns:
        print(f"Error: Column '{args.col}' not found in the input file.")
        exit(1)

    # --- Load Mapping File ---
    print(f"Loading mapping file: {args.map}")
    mimat_to_name_map = {}
    try:
        map_df = pd.read_csv(args.map)
        # Ensure required columns exist in map file
        if 'MIMAT_ID' not in map_df.columns or 'miRNA_name' not in map_df.columns:
            print(f"Error: Mapping file must contain 'MIMAT_ID' and 'miRNA_name' columns.")
            exit(1)
        # Create dictionary, handle potential NaN in miRNA_name
        mimat_to_name_map = pd.Series(map_df.miRNA_name.values, index=map_df.MIMAT_ID).dropna().to_dict()
        print(f"Loaded {len(mimat_to_name_map)} mappings from file.")
    except FileNotFoundError:
        print(f"Warning: Mapping file not found at {args.map}. Will rely heavily on API lookup.")
    except Exception as e:
        print(f"Warning: Error reading mapping file {args.map}: {e}. Will rely heavily on API lookup.")
    
    if not mimat_to_name_map:
        print("Warning: No mappings loaded from the map file.")

    # --- Define Conversion Function with Fallback ---
    def get_name_with_fallback(mimat_id):
        # 1. Check local map first
        name = mimat_to_name_map.get(mimat_id)
        if name and name != mimat_id: # Ensure it's not just the ID itself if map fails
            return name
        
        # 2. If not in local map, try API
        api_name = get_mirna_name_from_api(mimat_id)
        if api_name:
            # Optionally update the map in memory for this run? (No, avoid side effects)
            return api_name
            
        # 3. If API fails or no name found, return original ID
        print(f"  Fallback: Using original ID '{mimat_id}' as name.")
        return mimat_id # Fallback to original ID

    # --- Add miRNA Name Column ---
    print(f"Adding 'miRNA_name' column based on '{args.col}' with API fallback...")
    df['miRNA_name'] = df[args.col].apply(get_name_with_fallback)
    
    # Reorder columns to put miRNA_name after the original ID column
    try:
        cols = df.columns.tolist()
        # Find the index of the original column
        id_col_index = cols.index(args.col)
        # Insert the new column right after it
        cols.insert(id_col_index + 1, cols.pop(cols.index('miRNA_name')))
        df = df[cols]
    except ValueError:
        # If column index logic fails for some reason, just proceed without reordering
        print(f"Warning: Could not reorder columns. 'miRNA_name' will be appended at the end.")

    # --- Save Output File ---
    print(f"Saving results with miRNA names to: {args.output}")
    try:
        df.to_csv(args.output, index=False)
    except Exception as e:
        print(f"Error writing output file {args.output}: {e}")
        exit(1)

    print("Conversion complete!")

if __name__ == "__main__":
    main() 