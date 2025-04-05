"""
Visualization module for RNA structures and miRNA binding sites.

This module provides functions to visualize RNA secondary structures
and miRNA binding sites, making it easier to analyze and interpret results.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
import subprocess
import tempfile
import os
from pathlib import Path

def visualize_rna_structure(sequence: str, structure: str, 
                          title: str = "RNA Secondary Structure",
                          output_file: Optional[str] = None,
                          show: bool = True) -> str:
    """
    Visualize RNA secondary structure using ViennaRNA's RNAplot.
    
    Args:
        sequence: RNA sequence
        structure: Secondary structure in dot-bracket notation
        title: Title for the plot
        output_file: Path to save the output file (SVG format)
        show: Whether to display the plot
        
    Returns:
        Path to the generated visualization file
    """
    # Create temporary input file
    fd, input_file = tempfile.mkstemp(suffix=".fa")
    try:
        with os.fdopen(fd, 'w') as tmp:
            tmp.write(f">{title}\n{sequence}\n{structure}\n")
            
        # Determine output file
        if output_file is None:
            output_file = f"{title.replace(' ', '_')}.svg"
            
        # Run RNAplot to generate visualization
        cmd = ["RNAplot", "--output-format", "svg"]
        subprocess.run(
            cmd,
            input=f">{title}\n{sequence}\n{structure}\n".encode(),
            capture_output=True,
            check=True
        )
        
        # RNAplot generates a file named "rna.svg" in the current directory
        # Rename it to the desired output file
        if os.path.exists("rna.svg"):
            os.rename("rna.svg", output_file)
        
        # Open the file if show is True
        if show:
            # On macOS, use 'open' command to display the SVG
            if os.uname().sysname == "Darwin":
                subprocess.run(["open", output_file])
                
        return output_file
        
    finally:
        # Clean up
        os.unlink(input_file)
        
def plot_binding_site_accessibility(sequence: str, 
                                  accessibility_profile: List[float],
                                  binding_sites: Optional[List[int]] = None,
                                  window_size: int = 10,
                                  title: str = "Binding Site Accessibility",
                                  output_file: Optional[str] = None) -> None:
    """
    Plot the accessibility of potential miRNA binding sites.
    
    Args:
        sequence: RNA sequence
        accessibility_profile: List of accessibility values for each position
        binding_sites: List of binding site positions to highlight
        window_size: Window size for moving average smoothing
        title: Title for the plot
        output_file: Path to save the output file
    """
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Plot accessibility profile
    positions = list(range(len(accessibility_profile)))
    plt.plot(positions, accessibility_profile, color='blue', alpha=0.7)
    
    # Apply moving average for smoother visualization
    if window_size > 0 and window_size < len(accessibility_profile):
        smooth_data = np.convolve(
            accessibility_profile, 
            np.ones(window_size)/window_size, 
            mode='valid'
        )
        smooth_positions = positions[window_size-1:len(smooth_data)+window_size-1]
        plt.plot(smooth_positions, smooth_data, color='red', linewidth=2)
    
    # Highlight binding sites if provided
    if binding_sites:
        for pos in binding_sites:
            plt.axvspan(pos, pos + 7, color='green', alpha=0.3)
    
    # Labels and title
    plt.xlabel('Position')
    plt.ylabel('Accessibility')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    # Save if output file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    plt.show()
    
def plot_mirna_mrna_interaction(mirna: str, mrna: str, 
                              binding_position: int,
                              title: str = "miRNA-mRNA Interaction",
                              output_file: Optional[str] = None) -> None:
    """
    Visualize the interaction between a miRNA and its target site on an mRNA.
    
    Args:
        mirna: miRNA sequence
        mrna: mRNA sequence
        binding_position: Position of the binding site in the mRNA
        title: Title for the plot
        output_file: Path to save the output file
    """
    # Create figure
    plt.figure(figsize=(10, 4))
    
    # Format sequences for visualization
    mirna_vis = list(mirna[::-1])  # Reverse the miRNA for 5'->3' orientation
    
    # Extract the binding region from mRNA
    binding_end = min(binding_position + len(mirna), len(mrna))
    mrna_binding = mrna[binding_position:binding_end]
    
    # Extend with dots if necessary
    if len(mrna_binding) < len(mirna):
        mrna_binding += '.' * (len(mirna) - len(mrna_binding))
    
    # Plot sequences
    for i, base in enumerate(mirna_vis):
        plt.text(i, 0, base, ha='center', va='center', fontsize=12, 
                color='blue', fontweight='bold')
        
        if i < len(mrna_binding):
            plt.text(i, 1, mrna_binding[i], ha='center', va='center', 
                    fontsize=12, color='red', fontweight='bold')
            
            # Draw connection line for complementary bases
            if (mirna_vis[i] == 'A' and mrna_binding[i] == 'U') or \
               (mirna_vis[i] == 'U' and mrna_binding[i] == 'A') or \
               (mirna_vis[i] == 'G' and mrna_binding[i] == 'C') or \
               (mirna_vis[i] == 'C' and mrna_binding[i] == 'G'):
                plt.plot([i, i], [0.2, 0.8], 'k-', linewidth=1)
    
    # Labels and title
    plt.text(-1, 0, "miRNA 5'", ha='right', va='center', fontsize=10)
    plt.text(len(mirna), 0, "3'", ha='left', va='center', fontsize=10)
    plt.text(-1, 1, "mRNA 3'", ha='right', va='center', fontsize=10)
    plt.text(len(mrna_binding), 1, "5'", ha='left', va='center', fontsize=10)
    
    plt.title(title)
    plt.axis('off')
    
    # Save if output file is provided
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    plt.show() 