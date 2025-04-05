#!/usr/bin/env python3
"""
Accessibility Predictor for miRNA Binding Sites

This tool predicts optimal binding site locations for miRNAs (particularly miR-122)
within target mRNA sequences by analyzing RNA secondary structure and accessibility.

Usage:
    python accessibility_predictor.py --utr <utr_sequence_or_file> --mirna <mirna_sequence> [options]

Author: Hackathon-Bio Team
"""

import os
import sys
import argparse
import json
import tempfile
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Add src directory to path to import our modules
sys.path.append(str(Path(__file__).parent.parent))

from src.models.rna_structure import RNAStructurePrediction
from src.utils.mirna_utils import (
    reverse_complement,
    get_seed_sequence,
    validate_rna_sequence,
    design_binding_site
)
from src.visualization.structure_viz import (
    visualize_rna_structure,
    plot_binding_site_accessibility,
    plot_mirna_mrna_interaction
)

# Default miRNA sequences (focusing on NASH-HCC relevant miRNAs)
DEFAULT_MIRNAS = {
    "miR-122-5p": "UGGAGUGUGACAAUGGUGUUUG",  # liver-specific, downregulated in NASH-HCC
    "miR-21-5p": "UAGCUUAUCAGACUGAUGUUGA",   # pro-oncogenic, upregulated in NASH-HCC
    "miR-34a-5p": "UGGCAGUGUCUUAGCUGGUUGU",  # complex role in NASH-HCC
    "miR-26a-5p": "UUCAAGUAAUCCAGGAUAGGCU"   # tumor suppressor, downregulated in NASH-HCC
}

# Example UTR sequences (for testing)
EXAMPLE_UTRS = {
    "test_utr": "GGGCUACUUUAAAGUGCUGCUAUAGAUACAAUGAGGAUGAGUACCUGAUCUAAAUAGGUACUUAA",
    "albumin_partial": "AUCUUACUGAGGUGAGAACUGUCCUUGCUCAGCUUGGAGAAGACAUUUAGCUAGCUUACUCAAUAAAUGCUUUUUAAGCA",
    "luciferase_partial": "UGAAGGGCUGGAAAAAUGCUUUAGCUGAAGCUGUAUGUGAGAUAUCUUUAAUAAAGAUCUUUGCCAAGAGGCA"
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Predict optimal miRNA binding sites based on accessibility")
    
    # Input parameters
    parser.add_argument("--utr", type=str, 
                      help="Target UTR sequence or path to FASTA file")
    parser.add_argument("--mirna", type=str, default="miR-122-5p",
                      help=f"miRNA name or sequence (default: miR-122-5p). Options: {', '.join(DEFAULT_MIRNAS.keys())}")
    
    # Analysis parameters
    parser.add_argument("--window-size", type=int, default=80,
                      help="Window size for local folding (default: 80)")
    parser.add_argument("--max-span", type=int, default=40,
                      help="Maximum base pair span (default: 40)")
    parser.add_argument("--temperature", type=float, default=37.0,
                      help="Temperature for folding (default: 37.0)")
    parser.add_argument("--num-sites", type=int, default=3,
                      help="Number of top binding sites to report (default: 3)")
    
    # Output parameters
    parser.add_argument("--output-dir", type=str, default="results",
                      help="Output directory for results (default: results)")
    parser.add_argument("--visualize", action="store_true",
                      help="Generate visualizations")
    parser.add_argument("--example", action="store_true",
                      help=f"Use example UTR (options: {', '.join(EXAMPLE_UTRS.keys())})")
    parser.add_argument("--example-name", type=str, default="test_utr",
                      help="Name of example UTR to use")
    
    return parser.parse_args()

def read_sequence(input_str):
    """
    Read RNA sequence from string, FASTA file, or example.
    
    Args:
        input_str: Input string that could be a sequence, file path, or example name
        
    Returns:
        RNA sequence string
    """
    # Check if it's a FASTA file
    if os.path.exists(input_str):
        with open(input_str, 'r') as f:
            lines = f.readlines()
            # Skip header lines starting with '>'
            seq_lines = [line.strip() for line in lines if not line.startswith('>')]
            return ''.join(seq_lines).upper()
    
    # Check if it's an example name
    if input_str in EXAMPLE_UTRS:
        return EXAMPLE_UTRS[input_str]
    
    # Assume it's a sequence
    return input_str.strip().upper()

def parse_dotbracket(dot_bracket):
    """
    Parse dot-bracket notation to identify paired and unpaired regions.
    
    Args:
        dot_bracket: RNA structure in dot-bracket notation
        
    Returns:
        List of booleans where True indicates unpaired nucleotides
    """
    return [char == '.' for char in dot_bracket]

def calculate_accessibility(unpaired_regions, position, window_size=7):
    """
    Calculate accessibility score for a potential binding site.
    
    Args:
        unpaired_regions: List of booleans where True indicates unpaired nucleotides
        position: Starting position for binding site
        window_size: Size of binding site (typically 7 for seed region)
        
    Returns:
        Accessibility score (0-1) where 1 is fully accessible
    """
    if position + window_size > len(unpaired_regions):
        return 0.0
    
    # Count unpaired nucleotides in the window
    site_window = unpaired_regions[position:position+window_size]
    accessibility = sum(site_window) / len(site_window)
    
    # Add bonus for unpaired flanking regions (2nt on each side)
    flank_start = max(0, position - 2)
    flank_end = min(len(unpaired_regions), position + window_size + 2)
    
    flank_left = unpaired_regions[flank_start:position] if position > flank_start else []
    flank_right = unpaired_regions[position+window_size:flank_end] if position+window_size < flank_end else []
    
    flank_accessibility = (sum(flank_left) + sum(flank_right)) / max(1, len(flank_left) + len(flank_right))
    
    # Weighted combination of site and flank accessibility
    return 0.8 * accessibility + 0.2 * flank_accessibility

def find_optimal_binding_sites(utr_seq, mirna_seq, num_sites=3):
    """
    Find optimal binding sites for a miRNA in a UTR sequence.
    
    Args:
        utr_seq: Target UTR sequence
        mirna_seq: miRNA sequence
        num_sites: Number of top sites to return
        
    Returns:
        List of dictionaries with site information
    """
    # Initialize RNA structure prediction
    rna_predictor = RNAStructurePrediction()
    
    # Predict UTR structure
    structure_result = rna_predictor.predict_structure(
        utr_seq,
        temperature=37.0,
        avoid_lonely_pairs=True,
        partition_function=True
    )
    
    # Extract structure in dot-bracket notation
    if 'mfe_structure' not in structure_result:
        print("Error: Unable to predict structure.")
        return []
    
    dot_bracket = structure_result['mfe_structure']
    unpaired_regions = parse_dotbracket(dot_bracket)
    
    # Get seed sequence from miRNA
    seed_seq = get_seed_sequence(mirna_seq)
    seed_length = len(seed_seq)
    
    # Calculate accessibility scores for all possible positions
    scores = []
    for pos in range(len(utr_seq) - seed_length + 1):
        accessibility = calculate_accessibility(unpaired_regions, pos, seed_length)
        
        # Simple sequence match bonus (prefer seed complementary regions)
        seed_complement = reverse_complement(seed_seq)
        site_seq = utr_seq[pos:pos+seed_length]
        
        # Calculate sequence match score
        match_score = sum(a == b for a, b in zip(site_seq, seed_complement)) / seed_length
        
        # Combined score (weighted accessibility and sequence match)
        combined_score = 0.7 * accessibility + 0.3 * match_score
        
        scores.append({
            'position': pos,
            'accessibility': accessibility,
            'sequence_match': match_score,
            'combined_score': combined_score,
            'site_sequence': site_seq,
            'local_structure': dot_bracket[max(0, pos-5):min(len(dot_bracket), pos+seed_length+5)]
        })
    
    # Sort by combined score and return top sites
    sorted_scores = sorted(scores, key=lambda x: x['combined_score'], reverse=True)
    return sorted_scores[:num_sites], structure_result

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Handle input sequences
    utr_seq = None
    mirna_seq = None
    
    # Get UTR sequence
    if args.example:
        utr_seq = EXAMPLE_UTRS.get(args.example_name, EXAMPLE_UTRS['test_utr'])
        print(f"Using example UTR: {args.example_name}")
    elif args.utr:
        utr_seq = read_sequence(args.utr)
    else:
        print("Error: No UTR sequence provided. Use --utr or --example.")
        return 1
    
    # Get miRNA sequence
    if args.mirna in DEFAULT_MIRNAS:
        mirna_seq = DEFAULT_MIRNAS[args.mirna]
        mirna_name = args.mirna
        print(f"Using {mirna_name}: {mirna_seq}")
    else:
        mirna_seq = args.mirna
        mirna_name = "custom_mirna"
        print(f"Using custom miRNA: {mirna_seq}")
    
    # Validate sequences
    if not validate_rna_sequence(utr_seq):
        print(f"Error: Invalid UTR RNA sequence. Must contain only A, U, G, C.")
        return 1
        
    if not validate_rna_sequence(mirna_seq):
        print(f"Error: Invalid miRNA RNA sequence. Must contain only A, U, G, C.")
        return 1
    
    print(f"UTR length: {len(utr_seq)} nt")
    print(f"miRNA length: {len(mirna_seq)} nt")
    
    # Find optimal binding sites
    print(f"Finding optimal binding sites for {mirna_name} in target UTR...")
    optimal_sites, structure_result = find_optimal_binding_sites(
        utr_seq, mirna_seq, num_sites=args.num_sites
    )
    
    # Display results
    print("\nTop binding sites:")
    for i, site in enumerate(optimal_sites):
        print(f"\nSite {i+1} (position {site['position']})")
        print(f"  Sequence: {site['site_sequence']}")
        print(f"  Accessibility score: {site['accessibility']:.3f}")
        print(f"  Sequence match score: {site['sequence_match']:.3f}")
        print(f"  Combined score: {site['combined_score']:.3f}")
        print(f"  Local structure: {site['local_structure']}")
    
    # Generate accessibility profile
    accessibility_profile = [
        calculate_accessibility(parse_dotbracket(structure_result['mfe_structure']), pos, 7)
        for pos in range(len(utr_seq) - 7 + 1)
    ]
    
    # Save results to JSON
    results = {
        "utr_sequence": utr_seq,
        "mirna_name": mirna_name,
        "mirna_sequence": mirna_seq,
        "structure": structure_result['mfe_structure'],
        "mfe_energy": structure_result.get('mfe_energy', 0),
        "optimal_sites": optimal_sites,
        "accessibility_profile": accessibility_profile
    }
    
    with open(output_dir / f"{mirna_name}_binding_analysis.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {output_dir / f'{mirna_name}_binding_analysis.json'}")
    
    # Generate visualizations if requested
    if args.visualize:
        # Plot accessibility profile
        site_positions = [site['position'] for site in optimal_sites]
        plot_binding_site_accessibility(
            utr_seq,
            accessibility_profile + [0] * 6,  # Pad with zeros to match UTR length
            binding_sites=site_positions,
            title=f"{mirna_name} Binding Site Accessibility",
            output_file=str(output_dir / f"{mirna_name}_accessibility.png")
        )
        
        # Plot miRNA-mRNA interaction for top site
        if optimal_sites:
            plot_mirna_mrna_interaction(
                mirna_seq,
                utr_seq,
                optimal_sites[0]['position'],
                title=f"{mirna_name} Binding to Target UTR",
                output_file=str(output_dir / f"{mirna_name}_interaction.png")
            )
        
        # Plot RNA structure
        try:
            visualize_rna_structure(
                utr_seq,
                structure_result['mfe_structure'],
                title="UTR Secondary Structure",
                output_file=str(output_dir / "utr_structure.svg"),
                show=False
            )
        except Exception as e:
            print(f"Warning: Could not generate structure visualization: {e}")
        
        print(f"Visualizations saved to {output_dir}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 