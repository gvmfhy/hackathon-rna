#!/usr/bin/env python3
"""
Example script demonstrating miRNA target site analysis and optimization.

This script shows how to use the RNA structure prediction and miRNA utility
modules to analyze and optimize miRNA binding sites in mRNA sequences.
"""

import os
import sys
import json
from pathlib import Path

# Add the src directory to the path so we can import our modules
sys.path.append(str(Path(__file__).parent.parent))

from src.models.rna_structure import RNAStructurePrediction
from src.utils.mirna_utils import (
    reverse_complement, 
    get_seed_sequence,
    find_seed_matches,
    evaluate_binding_site,
    design_binding_site,
    validate_rna_sequence
)

# Try to import matplotlib visualization if available
try:
    import matplotlib
    import numpy as np
    from src.visualization.structure_viz import (
        visualize_rna_structure,
        plot_binding_site_accessibility,
        plot_mirna_mrna_interaction
    )
    visualization_available = True
except ImportError:
    visualization_available = False
    print("Matplotlib not available. Visualization will be skipped.")

# Create output directory if it doesn't exist
output_dir = Path(__file__).parent / "data" / "output"
output_dir.mkdir(exist_ok=True, parents=True)

def main():
    print("miRNA Target Site Analysis and Optimization Example")
    print("-" * 60)
    
    # Example miRNA and mRNA sequences
    mirna_sequence = "UAGCUUAUCAGACUGAUGUUGA"  # Example: miR-122-5p (liver-specific)
    mrna_sequence = "GGGCUACUUUAAAGUGCUGCUAUAGAUACAAUGAGGAUGAGUACCUGAUCUAAAUAGGUACUUAA"
    
    # Validate sequences
    if not validate_rna_sequence(mirna_sequence):
        print(f"Invalid miRNA sequence: {mirna_sequence}")
        return
        
    if not validate_rna_sequence(mrna_sequence):
        print(f"Invalid mRNA sequence: {mrna_sequence}")
        return
        
    print(f"miRNA: {mirna_sequence}")
    print(f"mRNA: {mrna_sequence}")
    print()
    
    # 1. Basic sequence analysis
    seed = get_seed_sequence(mirna_sequence)
    print(f"miRNA seed sequence: {seed}")
    print(f"Reverse complement: {reverse_complement(mirna_sequence)}")
    print()
    
    # 2. Find potential binding sites
    binding_sites = find_seed_matches(mrna_sequence, mirna_sequence)
    print(f"Found {len(binding_sites)} potential binding sites at positions: {binding_sites}")
    print()
    
    # 3. RNA structure prediction
    rna_structure = RNAStructurePrediction()
    structure_result = rna_structure.predict_structure(
        mrna_sequence,
        temperature=37.0,
        avoid_lonely_pairs=True,
        partition_function=True
    )
    
    print("RNA Structure Prediction Results:")
    print(f"MFE Structure: {structure_result.get('mfe_structure', 'N/A')}")
    print(f"MFE Energy: {structure_result.get('mfe_energy', 'N/A')} kcal/mol")
    
    if 'ensemble_energy' in structure_result:
        print(f"Ensemble Energy: {structure_result.get('ensemble_energy')} kcal/mol")
    
    if 'mfe_frequency' in structure_result:
        print(f"MFE Frequency: {structure_result.get('mfe_frequency')}")
        
    if 'ensemble_diversity' in structure_result:
        print(f"Ensemble Diversity: {structure_result.get('ensemble_diversity')}")
    print()
    
    # 4. Design optimal binding sites
    print("Designing optimal binding sites:")
    perfect_site = design_binding_site(mirna_sequence, perfect_complement=True)
    seed_site = design_binding_site(mirna_sequence, perfect_complement=False)
    
    print(f"Perfect complement site: {perfect_site}")
    print(f"Seed-based site: {seed_site}")
    print()
    
    # 5. Evaluate binding sites
    if binding_sites:
        print("Evaluating binding sites:")
        for pos in binding_sites:
            evaluation = evaluate_binding_site(mrna_sequence, mirna_sequence, pos)
            print(f"Site at position {pos}:")
            print(f"  Sequence with flanks: {evaluation['site_with_flanks']}")
            print(f"  Binding score: {evaluation['binding_score']}")
            print(f"  Normalized score: {evaluation['normalized_score']:.2f}")
    print()
    
    # 6. Visualization (if available)
    if visualization_available and 'mfe_structure' in structure_result:
        try:
            # Visualize RNA structure
            svg_file = visualize_rna_structure(
                mrna_sequence, 
                structure_result['mfe_structure'],
                title="mRNA Secondary Structure",
                output_file=str(output_dir / "mrna_structure.svg"),
                show=False
            )
            print(f"Structure visualization saved to: {svg_file}")
            
            # Visualize miRNA-mRNA interaction if binding sites exist
            if binding_sites:
                plot_mirna_mrna_interaction(
                    mirna_sequence,
                    mrna_sequence,
                    binding_sites[0],
                    title="miRNA-mRNA Interaction",
                    output_file=str(output_dir / "mirna_mrna_interaction.png")
                )
                print(f"Interaction visualization saved to: {output_dir / 'mirna_mrna_interaction.png'}")
                
            # Create a simple accessibility profile (mock data for demonstration)
            accessibility = [0.5 + 0.3 * (i % 3 - 1) for i in range(len(mrna_sequence))]
            plot_binding_site_accessibility(
                mrna_sequence,
                accessibility,
                binding_sites=binding_sites,
                title="mRNA Accessibility Profile",
                output_file=str(output_dir / "accessibility_profile.png")
            )
            print(f"Accessibility profile saved to: {output_dir / 'accessibility_profile.png'}")
        except Exception as e:
            print(f"Error during visualization: {e}")
    
    # 7. Save results
    result_data = {
        "mirna": mirna_sequence,
        "mrna": mrna_sequence,
        "structure_prediction": structure_result,
        "binding_sites": binding_sites,
        "designed_sites": {
            "perfect_complement": perfect_site,
            "seed_based": seed_site
        }
    }
    
    # Save to JSON
    with open(output_dir / "example_results.json", "w") as f:
        json.dump(result_data, f, indent=2)
    
    print(f"Results saved to {output_dir / 'example_results.json'}")
    
if __name__ == "__main__":
    main() 