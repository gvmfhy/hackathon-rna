import os
import subprocess
import pandas as pd
import glob
import time # Use time instead of datetime for timestamp generation if datetime not needed otherwise
import sys

# Try importing Biopython, needed for FASTA parsing
try:
    from Bio import SeqIO
except ImportError:
    print("Error: Biopython not installed. Please install it: pip install biopython")
    sys.exit(1)

# Define base directory and script paths
script_dir = os.path.dirname(os.path.abspath(__file__))
# base_dir = os.path.dirname(script_dir) # Original line - likely incorrect
base_dir = script_dir # Assume this orchestrator script is in the hackathonbro/ directory

de_script_path = os.path.join(base_dir, 'scripts', 'mirna_differential_expression.py')
accessibility_script_path = os.path.join(base_dir, 'src', 'accessibility_predictor.py')
prioritization_script_path = os.path.join(base_dir, 'scripts', 'target_prioritization_analysis.py')
off_target_script_path = os.path.join(base_dir, 'scripts', 'off_target_check.py')
results_dir = os.path.join(base_dir, 'results')
mature_fasta_path = os.path.join(base_dir, 'data', 'annotations', 'mature.fa') # Updated path

# Import the default UTR sequences - use RNA versions for computational analysis
# Correct path to UTR sequences needs base_dir context
sys.path.append(os.path.join(base_dir, 'data', 'utr_sequences'))
from hbb_utr_sequences import HBB_3UTR_RNA_SHORT # Import relative to the updated sys.path

def parse_mature_fasta(fasta_file):
    """Parses mature.fa to create a MIMAT ID to sequence mapping for HUMAN miRNAs only."""
    mimat_to_seq = {}
    if not os.path.exists(fasta_file):
        print(f"Warning: miRBase mature sequence file not found at {fasta_file}")
        print("Accessibility predictor step will likely fail without sequences.")
        return mimat_to_seq
    
    print(f"Parsing HUMAN miRNA sequences from: {fasta_file}")
    human_mirna_count = 0
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Check if the ID starts with 'hsa-' for Homo sapiens
            if record.id.startswith("hsa-"):
                human_mirna_count += 1
                # Header format usually >hsa-miR-XXX MIMATYYYYYYY ...
                parts = record.description.split()
                mimat_id = None
                for part in parts:
                    if part.startswith("MIMAT"):
                        mimat_id = part
                        break
                if mimat_id:
                    # Convert sequence to string and ensure it's RNA (U instead of T)
                    sequence = str(record.seq).upper().replace("T", "U")
                    mimat_to_seq[mimat_id] = sequence
                else:
                    # Handle cases where MIMAT ID might be missing in header for an hsa entry
                    print(f"  Warning: Found human miRNA '{record.id}' but no MIMAT ID in header: {record.description}")
    except Exception as e:
        print(f"Error parsing {fasta_file}: {e}")
        # Return potentially partially filled dict
    
    print(f"Found {human_mirna_count} human miRNA entries.")
    print(f"Loaded sequences for {len(mimat_to_seq)} MIMAT IDs.")
    return mimat_to_seq

# --- Step 1: Run miRNA Differential Expression ---
print("--- Step 1: Running miRNA Differential Expression Analysis ---")
print(f"Executing: {de_script_path}...")

# Ensure results directory exists
os.makedirs(results_dir, exist_ok=True)

# Execute the DE script
# Run scripts from the base directory context
process_de = subprocess.run(['python', de_script_path], capture_output=True, text=True, cwd=base_dir)
print("STDOUT:")
print(process_de.stdout)
print("STDERR:")
print(process_de.stderr)

if process_de.returncode != 0:
    print("\nError running miRNA Differential Expression script. Exiting.")
    exit(1)

print("\nDifferential Expression analysis complete.")

# --- Identify Downregulated Candidate miRNAs ---
print("\n--- Identifying Downregulated Candidate miRNAs ---")

# Find the latest significant results CSV file
# Add a small delay to ensure the file system registers the new file completely
time.sleep(1) 
significant_files = glob.glob(os.path.join(results_dir, 'mirna_de_significant_*.csv'))
if not significant_files:
    print("Error: No significant miRNA DE results file found. Exiting.")
    exit(1)

downregulated_mirna_ids = []
try:
    latest_sig_file = max(significant_files, key=os.path.getctime)
    print(f"Reading candidates from: {latest_sig_file}")
    sig_df = pd.read_csv(latest_sig_file)
    if 'log2FC' not in sig_df.columns or 'miRNA' not in sig_df.columns:
        print(f"Error: 'log2FC' or 'miRNA' column not found in {latest_sig_file}. Exiting.")
        exit(1)
    downregulated_mirna_ids = sig_df[sig_df['log2FC'] < 0]['miRNA'].tolist()
except Exception as e:
    print(f"Error reading or processing significant miRNA file: {e}")
    exit(1)

if not downregulated_mirna_ids:
    print("No downregulated miRNAs found in the significant results. Skipping Accessibility and Off-Target steps.")
    exit(0) 

print(f"Found {len(downregulated_mirna_ids)} downregulated miRNAs (MIMAT IDs) as candidates:")
for mirna_id in downregulated_mirna_ids:
    print(f" - {mirna_id}")

# --- NEW STEP: Convert MIMAT IDs to Names in DE Results ---
print("\n--- Preparing DE Results with miRNA Names for Downstream Analysis ---")
conversion_script_path = os.path.join(base_dir, 'scripts', 'convert_mimat_to_mirna.py')
mapping_file_path = os.path.join(base_dir, 'data', 'annotations', 'mimat_to_mirna_mapping.csv')
input_de_file = latest_sig_file # Use the latest significant file found earlier
output_de_named_file = os.path.join(results_dir, 'mirna_de_significant_with_names.csv')

if not os.path.exists(conversion_script_path):
    print(f"Error: Conversion script not found at {conversion_script_path}. Cannot create named DE file.")
    exit(1)
if not os.path.exists(mapping_file_path):
    print(f"Error: MIMAT mapping file not found at {mapping_file_path} for conversion script.")
    exit(1)

print(f"Running conversion script: {conversion_script_path}")
print(f"  Input DE file: {os.path.basename(input_de_file)}")
print(f"  Output named file: {os.path.basename(output_de_named_file)}")

cmd_convert = [
    'python',
    conversion_script_path,
    '--input', input_de_file,
    '--output', output_de_named_file,
    '--map', mapping_file_path,
    '--col', 'miRNA' # Specify the column containing MIMAT IDs
]

process_convert = subprocess.run(cmd_convert, capture_output=True, text=True, cwd=base_dir)
print(f"STDOUT: {process_convert.stdout.strip()}")
print(f"STDERR: {process_convert.stderr.strip()}")

if process_convert.returncode != 0 or not os.path.exists(output_de_named_file):
    print(f"\nError running conversion script or output file not generated. Check logs. Exiting.")
    exit(1)

print("Successfully generated DE results file with miRNA names.")

# --- NEW STEP: Fetch miRNA Sequences ---
print("\n--- Fetching miRNA Sequences for Candidates ---")
mimat_to_sequence_map = parse_mature_fasta(mature_fasta_path)

candidate_sequences = {}
missing_sequences = []
for mirna_id in downregulated_mirna_ids:
    sequence = mimat_to_sequence_map.get(mirna_id)
    if sequence:
        candidate_sequences[mirna_id] = sequence
    else:
        print(f"  Warning: Sequence not found in mature.fa for {mirna_id}")
        missing_sequences.append(mirna_id)

print(f"Found sequences for {len(candidate_sequences)} out of {len(downregulated_mirna_ids)} candidates.")
if not candidate_sequences:
    print("Error: No sequences found for any candidate miRNAs. Cannot proceed with Accessibility Prediction.")
    exit(1)

# --- Step 2 (Now Step 3): Accessibility Prediction (Using HBB 3' UTR and Sequences) ---
print("\n--- Step 3: Running Accessibility Prediction (Using HBB 3' UTR and miRNA Sequences) ---")
print(f"Using HBB 3' UTR RNA sequence: {HBB_3UTR_RNA_SHORT[:30]}...")
print(f"Executing: {accessibility_script_path} for each candidate miRNA sequence...")

accessibility_results_dir = os.path.join(results_dir, "accessibility_analysis")
os.makedirs(accessibility_results_dir, exist_ok=True)

accessibility_success_count = 0
# Loop through the candidates for which we found sequences
for mirna_id, mirna_sequence in candidate_sequences.items():
    print(f"  Running for miRNA: {mirna_id} (Sequence: {mirna_sequence[:15]}...)")
    
    # Check how accessibility_predictor.py expects sequence input.
    # Assuming it uses --mirna for either name or sequence, we pass the sequence.
    cmd = [
        'python',
        accessibility_script_path,
        '--utr', HBB_3UTR_RNA_SHORT,
        '--mirna', mirna_sequence, # Pass the actual RNA sequence
        '--output-dir', accessibility_results_dir,
        # '--visualize' # Optional: uncomment to generate plots
    ]

    # Run the accessibility predictor script
    process_acc = subprocess.run(cmd, capture_output=True, text=True, cwd=base_dir)

    # Limit output printing to avoid clutter
    stdout_snippet = process_acc.stdout.strip().replace('\n', ' ')[-200:]
    stderr_snippet = process_acc.stderr.strip().replace('\n', ' ')

    print(f"    STDOUT (sample): ...{stdout_snippet}" )
    if stderr_snippet:
        print(f"    STDERR: {stderr_snippet}")

    if process_acc.returncode == 0:
        print(f"    Accessibility analysis successful for {mirna_id}.")
        accessibility_success_count += 1
    else:
        # Error is less likely now, but could still happen (e.g., RNAfold issues)
        print(f"    Warning: Accessibility analysis failed or produced errors for {mirna_id} (Return Code: {process_acc.returncode}). Check STDERR.")

print(f"\nAccessibility analysis attempted for {len(candidate_sequences)} miRNAs, {accessibility_success_count} succeeded.")

# --- Step 3 (Now Step 4): Off-Target Context Analysis ---
print("\n--- Step 4: Running Off-Target Context Analysis ---")
print(f"Executing: {prioritization_script_path}...")

# Execute the target prioritization script
process_prio = subprocess.run(['python', prioritization_script_path], capture_output=True, text=True, cwd=base_dir)
print("STDOUT:")
print(process_prio.stdout)
print("STDERR:")
print(process_prio.stderr)

if process_prio.returncode != 0:
    print("\nError running Target Prioritization script. Exiting.")
    exit(1)

print("\nOff-target context analysis complete.")

# --- Step 4 (Now Step 5): Safety Check with Contextual Off-Target Scoring ---
print("\n--- Step 5: Running Safety Check with Contextual Off-Target Scoring ---")
print(f"Executing: {off_target_script_path}...")

# Execute the off-target check script
cmd = [
    'python',
    off_target_script_path,
    '--utr', HBB_3UTR_RNA_SHORT,  # Use the HBB 3' UTR RNA sequence
    '--output-dir', results_dir,
    '--detailed'  # Get detailed output
]

process_safety = subprocess.run(cmd, capture_output=True, text=True, cwd=base_dir)
print("STDOUT:")
print(process_safety.stdout)
print("STDERR:")
print(process_safety.stderr)

if process_safety.returncode != 0:
    print("\nError running Off-Target Check script. Exiting.")
    exit(1)

print("\nSafety check with contextual off-target scoring complete.")

# --- Workflow Complete ---
print("\n--- Complete miRNA Switch Platform Workflow Integrated ---")
print("Successfully ran DE analysis, accessibility prediction, off-target context analysis, and safety check.")
print("All components for miRNA switch design and safety assessment are now integrated.")
print(f"All results are available in: {results_dir}") 