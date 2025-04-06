#!/usr/bin/env python3
"""
Off-Target Check with Contextual Scoring for miRNA Switch Technology

This script analyzes a therapeutic mRNA's UTR for potential binding sites of high-risk off-target miRNAs.
It provides a safety assessment by scoring the risk of unintended miRNA binding based on:
1. Seed match type
2. Site accessibility (RNA structure)
3. miRNA abundance context

Usage:
    python off_target_check.py [options]

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import json
import tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd

# Add parent directory to sys.path to import our modules
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

# Import our utility functions and modules
from src.models.rna_structure import RNAStructurePrediction
from src.utils.mirna_utils import (
    reverse_complement,
    get_seed_sequence,
    validate_rna_sequence,
    find_seed_matches
)

# Import the default UTR sequences - use RNA version for analysis
sys.path.append(os.path.join(parent_dir, 'data'))
from utr_sequences.hbb_utr_sequences import HBB_3UTR_RNA_SHORT

# High-risk miRNAs (highly abundant in healthy liver)
# These are miRNAs we want to avoid binding to our therapeutic mRNA
HIGH_RISK_MIRNAS = {
    "hsa-miR-122-5p": "UGGAGUGUGACAAUGGUGUUUG",  # liver-specific, highly abundant
    "hsa-miR-192-5p": "CUGACCUAUGAAUUGACAGCC"    # also abundant in healthy liver
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Off-Target Check with Contextual Scoring")
    
    # Input parameters
    parser.add_argument("--utr", type=str, 
                      help="Target UTR sequence or path to FASTA file (default: HBB 3' UTR RNA sequence)")
    parser.add_argument("--mirnas", type=str, 
                      help="Path to file with high-risk miRNAs (default: built-in liver miRNAs)")
    
    # Analysis parameters
    parser.add_argument("--seed-only", action="store_true",
                      help="Only consider seed region matches (positions 2-8)")
    parser.add_argument("--accessibility-threshold", type=float, default=0.6,
                      help="Threshold for considering a site accessible (default: 0.6)")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="../results",
                      help="Output directory for results (default: ../results)")
    parser.add_argument("--detailed", action="store_true",
                      help="Generate detailed report with all potential sites")
    
    return parser.parse_args()

def read_sequence(input_str):
    """
    Read RNA sequence from string or FASTA file.
    
    Args:
        input_str: Input string that could be a sequence or file path
        
    Returns:
        RNA sequence string
    """
    # Check if it's a FASTA file
    if input_str and os.path.exists(input_str):
        with open(input_str, 'r') as f:
            lines = f.readlines()
            # Skip header lines starting with '>'
            seq_lines = [line.strip() for line in lines if not line.startswith('>')]
            return ''.join(seq_lines).upper()
    
    # If not a file or no input provided, return the input string or default HBB 3' UTR
    return input_str.strip().upper() if input_str else HBB_3UTR_RNA_SHORT

def read_mirnas(file_path):
    """
    Read miRNAs from a file or use the default high-risk ones.
    
    Args:
        file_path: Path to file with miRNA sequences
        
    Returns:
        Dictionary of miRNA names to sequences
    """
    if file_path and os.path.exists(file_path):
        mirnas = {}
        with open(file_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        mirnas[parts[0]] = parts[1].upper()
        return mirnas
    
    # Default to built-in high-risk miRNAs
    return HIGH_RISK_MIRNAS

def calculate_accessibility(dot_bracket, position, window_size=7):
    """
    Calculate accessibility score for a potential binding site.
    
    Args:
        dot_bracket: RNA structure in dot-bracket notation
        position: Starting position for binding site
        window_size: Size of binding site (typically 7 for seed region)
        
    Returns:
        Accessibility score (0-1) where 1 is fully accessible
    """
    if position + window_size > len(dot_bracket):
        return 0.0
    
    # Parse dot-bracket to identify unpaired nucleotides (represented by ".")
    unpaired_regions = [char == '.' for char in dot_bracket]
    
    # Count unpaired nucleotides in the window
    site_window = unpaired_regions[position:position+window_size]
    accessibility = sum(site_window) / len(site_window)
    
    return accessibility

def analyze_seed_match(match_type):
    """
    Analyze seed match type and assign a base risk score.
    
    Args:
        match_type: Type of seed match (e.g., '7mer-A1', '7mer-m8', '8mer')
        
    Returns:
        Base risk score (0-1)
    """
    # For now, simple scoring based on match type
    if match_type == '8mer':
        return 1.0  # Highest risk
    elif match_type == '7mer-m8':
        return 0.9
    elif match_type == '7mer-A1':
        return 0.8
    elif match_type == '6mer':
        return 0.6
    else:
        return 0.4  # Other match types

def determine_match_type(mrna_seq, mirna_seq, position):
    """
    Determine the type of seed match at a given position.
    
    Args:
        mrna_seq: Target mRNA sequence
        mirna_seq: miRNA sequence
        position: Position of the match in the mRNA
        
    Returns:
        Match type description (e.g., '8mer', '7mer-m8', etc.)
    """
    seed = get_seed_sequence(mirna_seq)  # positions 2-8
    seed_rc = reverse_complement(seed)
    
    # Check if the position is valid for checking
    if position + len(seed_rc) > len(mrna_seq):
        return "Invalid position"
    
    # Get the sequence at the position
    site_seq = mrna_seq[position:position+len(seed_rc)]
    
    # Check if there's a perfect seed match
    if site_seq == seed_rc:
        # Check for 8mer (seed match + position 8 match + A at position 1)
        if (position > 0 and 
            mrna_seq[position-1] == 'A' and 
            position + len(seed_rc) < len(mrna_seq) and
            mirna_seq[8] and 
            mrna_seq[position+len(seed_rc)] == reverse_complement(mirna_seq[8])):
            return "8mer"
        # Check for 7mer-m8 (seed match + position 8 match)
        elif (position + len(seed_rc) < len(mrna_seq) and
             mirna_seq[8] and 
             mrna_seq[position+len(seed_rc)] == reverse_complement(mirna_seq[8])):
            return "7mer-m8"
        # Check for 7mer-A1 (seed match + A at position 1)
        elif position > 0 and mrna_seq[position-1] == 'A':
            return "7mer-A1"
        else:
            return "7mer"
    
    # Check for 6mer (positions 2-7)
    elif site_seq[:-1] == seed_rc[:-1]:
        return "6mer"
    
    # Partial matches or other types
    else:
        match_count = sum(1 for a, b in zip(site_seq, seed_rc) if a == b)
        return f"Partial ({match_count}/{len(seed_rc)} matches)"

def main():
    """Main function."""
    args = parse_arguments()
    
    # Set up output directory
    if args.output_dir.startswith(".."):
        # If relative path, make it relative to the script directory
        output_dir = os.path.join(script_dir, args.output_dir)
    else:
        output_dir = args.output_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get current timestamp for file naming
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Read target UTR sequence
    utr_seq = read_sequence(args.utr)
    print(f"Using UTR sequence of length {len(utr_seq)} nt")
    
    # Validate sequence
    if not validate_rna_sequence(utr_seq):
        print(f"Error: Invalid UTR RNA sequence. Must contain only A, U, G, C.")
        return 1
        
    # Read high-risk miRNAs
    mirnas = read_mirnas(args.mirnas)
    print(f"Analyzing {len(mirnas)} high-risk miRNAs")
    
    # Predict UTR structure
    print("Predicting UTR secondary structure...")
    rna_predictor = RNAStructurePrediction()
    structure_result = rna_predictor.predict_structure(
        utr_seq,
        temperature=37.0,
        avoid_lonely_pairs=True,
        partition_function=True
    )
    
    # Extract structure in dot-bracket notation
    if 'mfe_structure' not in structure_result:
        print("Error: Unable to predict structure.")
        return 1
    
    dot_bracket = structure_result['mfe_structure']
    
    # Initialize results data structure
    results = {
        "utr_sequence": utr_seq,
        "utr_length": len(utr_seq),
        "structure_energy": structure_result.get('mfe_energy', 0),
        "high_risk_mirnas": list(mirnas.keys()),
        "potential_off_target_sites": []
    }
    
    # Initialize summary statistics
    high_risk_count = 0
    medium_risk_count = 0
    low_risk_count = 0
    
    # Check for potential binding sites for each high-risk miRNA
    print("\nScanning for potential off-target binding sites...")
    for mirna_name, mirna_seq in mirnas.items():
        print(f"\nAnalyzing miRNA: {mirna_name}")
        
        # Find seed matches
        matches = find_seed_matches(utr_seq, mirna_seq)
        
        if not matches:
            print(f"  No seed matches found for {mirna_name}")
            continue
            
        print(f"  Found {len(matches)} potential binding sites")
        
        # Analyze each match
        for position in matches:
            # Determine match type
            match_type = determine_match_type(utr_seq, mirna_seq, position)
            
            # Calculate accessibility
            accessibility = calculate_accessibility(dot_bracket, position, window_size=7)
            
            # Calculate risk score based on match type and accessibility
            # Base score from match type
            base_score = analyze_seed_match(match_type)
            
            # Adjust for accessibility (high accessibility increases risk)
            if accessibility >= args.accessibility_threshold:
                # High accessibility - full risk
                risk_score = base_score
                accessibility_impact = "High"
            elif accessibility >= args.accessibility_threshold / 2:
                # Medium accessibility - reduced risk
                risk_score = base_score * 0.7
                accessibility_impact = "Medium"
            else:
                # Low accessibility - greatly reduced risk
                risk_score = base_score * 0.3
                accessibility_impact = "Low"
            
            # Determine risk category
            if risk_score >= 0.7:
                risk_category = "High"
                high_risk_count += 1
            elif risk_score >= 0.4:
                risk_category = "Medium"
                medium_risk_count += 1
            else:
                risk_category = "Low"
                low_risk_count += 1
            
            # Get local sequence context
            context_start = max(0, position - 5)
            context_end = min(len(utr_seq), position + 12)
            sequence_context = utr_seq[context_start:context_end]
            
            # Store the result
            site_info = {
                "mirna": mirna_name,
                "position": position,
                "match_type": match_type,
                "accessibility": accessibility,
                "accessibility_impact": accessibility_impact,
                "risk_score": risk_score,
                "risk_category": risk_category,
                "sequence_context": sequence_context,
                "local_structure": dot_bracket[context_start:context_end]
            }
            
            results["potential_off_target_sites"].append(site_info)
            
            # Print finding if it's high risk or detailed output is requested
            if risk_category == "High" or args.detailed:
                print(f"  Site at position {position}:")
                print(f"    Match type: {match_type}")
                print(f"    Accessibility: {accessibility:.2f} ({accessibility_impact})")
                print(f"    Risk score: {risk_score:.2f} ({risk_category} Risk)")
                print(f"    Sequence: {sequence_context}")
                print(f"    Structure: {dot_bracket[context_start:context_end]}")
    
    # Add summary statistics
    results["summary"] = {
        "total_potential_sites": len(results["potential_off_target_sites"]),
        "high_risk_sites": high_risk_count,
        "medium_risk_sites": medium_risk_count,
        "low_risk_sites": low_risk_count
    }
    
    # Sort by risk score (highest first)
    results["potential_off_target_sites"].sort(key=lambda x: x["risk_score"], reverse=True)
    
    # Save results to JSON
    output_file = os.path.join(output_dir, f"off_target_analysis_{timestamp}.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    # Generate a summary report
    summary_file = os.path.join(output_dir, f"off_target_summary_{timestamp}.txt")
    with open(summary_file, "w") as f:
        f.write("# Off-Target Binding Site Analysis\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"UTR sequence length: {len(utr_seq)} nt\n")
        f.write(f"miRNAs analyzed: {', '.join(mirnas.keys())}\n\n")
        
        f.write("## Summary\n")
        f.write(f"Total potential binding sites: {len(results['potential_off_target_sites'])}\n")
        f.write(f"High risk sites: {high_risk_count}\n")
        f.write(f"Medium risk sites: {medium_risk_count}\n")
        f.write(f"Low risk sites: {low_risk_count}\n\n")
        
        f.write("## High Risk Sites\n")
        high_risk_sites = [site for site in results["potential_off_target_sites"] 
                           if site["risk_category"] == "High"]
        
        if high_risk_sites:
            for i, site in enumerate(high_risk_sites, 1):
                f.write(f"{i}. miRNA: {site['mirna']}, Position: {site['position']}, "
                       f"Match: {site['match_type']}, "
                       f"Accessibility: {site['accessibility']:.2f}\n")
                f.write(f"   Sequence context: {site['sequence_context']}\n")
                f.write(f"   Structure context: {site['local_structure']}\n\n")
        else:
            f.write("No high risk sites identified.\n\n")
        
        f.write("## Safety Assessment\n")
        if high_risk_count > 0:
            f.write("⚠️ POTENTIAL SAFETY CONCERN: High-risk binding sites detected for liver-abundant miRNAs.\n")
            f.write("Recommendation: Consider modifying the UTR design to disrupt these sites.\n")
        elif medium_risk_count > 0:
            f.write("⚠️ MONITOR: Medium-risk binding sites detected.\n")
            f.write("Recommendation: Monitor during preclinical testing. May not require redesign.\n")
        else:
            f.write("✅ LOW RISK: No high or medium risk binding sites detected for liver-abundant miRNAs.\n")
            f.write("Recommendation: Proceed with current design.\n")

    # Print summary to console
    print("\n--- Off-Target Analysis Summary ---")
    print(f"Total potential binding sites: {len(results['potential_off_target_sites'])}")
    print(f"High risk sites: {high_risk_count}")
    print(f"Medium risk sites: {medium_risk_count}")
    print(f"Low risk sites: {low_risk_count}")
    
    print(f"\nDetailed results saved to: {output_file}")
    print(f"Summary report saved to: {summary_file}")
    
    # Safety assessment
    print("\n--- Safety Assessment ---")
    if high_risk_count > 0:
        print("⚠️ POTENTIAL SAFETY CONCERN: High-risk binding sites detected for liver-abundant miRNAs.")
        print("Recommendation: Consider modifying the UTR design to disrupt these sites.")
    elif medium_risk_count > 0:
        print("⚠️ MONITOR: Medium-risk binding sites detected.")
        print("Recommendation: Monitor during preclinical testing. May not require redesign.")
    else:
        print("✅ LOW RISK: No high or medium risk binding sites detected for liver-abundant miRNAs.")
        print("Recommendation: Proceed with current design.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 