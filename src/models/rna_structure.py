"""
RNA Secondary Structure Prediction Module.

This module provides functions to predict RNA secondary structures using the ViennaRNA package.
It serves as a wrapper around the RNAfold command-line tool.
"""

import subprocess
import os
import tempfile
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Union

class RNAStructurePrediction:
    """Class for RNA secondary structure prediction using ViennaRNA package."""
    
    def __init__(self, temp_dir: Optional[str] = None):
        """
        Initialize the RNA structure prediction class.
        
        Args:
            temp_dir: Optional directory to store temporary files. If None, system temp directory is used.
        """
        self.temp_dir = temp_dir if temp_dir else tempfile.gettempdir()
        
    def predict_structure(self, sequence: str, 
                         temperature: float = 37.0,
                         avoid_lonely_pairs: bool = True,
                         partition_function: bool = True) -> Dict:
        """
        Predict RNA secondary structure using RNAfold.
        
        Args:
            sequence: RNA sequence string
            temperature: Temperature in Celsius for folding calculation
            avoid_lonely_pairs: If True, avoid lonely pairs in structure prediction
            partition_function: If True, calculate partition function and base pair probabilities
            
        Returns:
            Dictionary containing structure prediction results
        """
        # Create command with appropriate parameters
        cmd = ["RNAfold"]
        cmd.extend(["-T", str(temperature)])
        
        if avoid_lonely_pairs:
            cmd.append("--noLP")
            
        if partition_function:
            cmd.append("-p")
        
        # Create temporary input file
        fd, input_file = tempfile.mkstemp(suffix=".fa", dir=self.temp_dir)
        try:
            with os.fdopen(fd, 'w') as tmp:
                tmp.write(f"{sequence}\n")
                
            # Run RNAfold using the file instead of piping the input
            with open(input_file, 'r') as infile:
                process = subprocess.run(
                    cmd,
                    stdin=infile,
                    capture_output=True,
                    text=True,
                    check=True
                )
            
            # Parse output
            result = self._parse_rnafold_output(process.stdout, sequence)
            return result
            
        finally:
            # Clean up
            os.unlink(input_file)
            # Also clean up any dot.ps files created in current directory
            ps_file = Path("dot.ps")
            if ps_file.exists():
                ps_file.unlink()
                
    def _parse_rnafold_output(self, output: str, sequence: str) -> Dict:
        """
        Parse the output from RNAfold.
        
        Args:
            output: String output from RNAfold
            sequence: Original RNA sequence
            
        Returns:
            Dictionary with parsed results
        """
        lines = output.strip().split('\n')
        result = {
            "sequence": sequence,
            "length": len(sequence)
        }
        
        # Extract MFE structure and energy
        for i, line in enumerate(lines):
            if i > 0 and '(' in line:  # Skip first line which is the sequence
                structure, energy_str = line.split(' (')
                result["mfe_structure"] = structure
                result["mfe_energy"] = float(energy_str.replace(')', ''))
                break
                
        # Extract ensemble information if available
        if len(lines) > 2 and '[' in lines[2]:
            ensemble_info = lines[2].split(' [')[1]
            result["ensemble_energy"] = float(ensemble_info.replace(']', ''))
            
        if len(lines) > 3 and 'frequency' in lines[3]:
            freq_info = lines[3].split('frequency of mfe structure in ensemble ')[1]
            freq_parts = freq_info.split(';')
            result["mfe_frequency"] = float(freq_parts[0])
            
            if 'ensemble diversity' in freq_parts[1]:
                diversity = freq_parts[1].strip().split('ensemble diversity ')[1]
                result["ensemble_diversity"] = float(diversity)
                
        return result
    
    def analyze_mirna_binding_site(self, mrna_seq: str, mirna_seq: str) -> Dict:
        """
        Analyze potential miRNA binding site in an mRNA sequence.
        
        Args:
            mrna_seq: Target mRNA sequence
            mirna_seq: microRNA sequence
            
        Returns:
            Dictionary with binding analysis results
        """
        # Use RNAcofold for this analysis if available
        # For now, this is a placeholder
        return {
            "mrna_seq": mrna_seq,
            "mirna_seq": mirna_seq,
            "binding_energy": 0.0,  # Placeholder
            "binding_position": 0,  # Placeholder
        }
    
    def analyze_site_accessibility(self, sequence: str, window_size: int = 80, 
                                 max_span: int = 40, unpaired_window: int = 10) -> Dict:
        """
        Analyze the accessibility of potential binding sites using RNAplfold.
        
        Args:
            sequence: RNA sequence
            window_size: Window size for local folding
            max_span: Maximum base pair span
            unpaired_window: Size of window for unpaired region
            
        Returns:
            Dictionary with accessibility analysis results
        """
        # This would use RNAplfold in a real implementation
        # For now, return placeholder
        return {
            "sequence": sequence,
            "accessibility_profile": [0.5] * len(sequence),  # Placeholder
        } 