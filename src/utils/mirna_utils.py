"""
Utility functions for miRNA sequence handling and binding site analysis.
"""

from typing import List, Dict, Tuple, Optional
import re

def reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of an RNA sequence.
    
    Args:
        sequence: RNA sequence string
        
    Returns:
        Reverse complement of the sequence
    """
    complement_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                      'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
    return ''.join([complement_dict.get(base, base) for base in reversed(sequence)])

def get_seed_sequence(mirna: str) -> str:
    """
    Extract the seed sequence (positions 2-8) from a miRNA.
    
    Args:
        mirna: miRNA sequence
        
    Returns:
        Seed sequence of the miRNA
    """
    if len(mirna) < 8:
        raise ValueError(f"miRNA sequence too short: {mirna}")
    return mirna[1:8]

def find_seed_matches(mrna: str, mirna: str) -> List[int]:
    """
    Find potential miRNA binding sites based on seed sequence matches.
    
    Args:
        mrna: Target mRNA sequence
        mirna: miRNA sequence
        
    Returns:
        List of starting positions of potential binding sites
    """
    seed = get_seed_sequence(mirna)
    seed_rc = reverse_complement(seed)
    
    matches = []
    for i in range(len(mrna) - len(seed_rc) + 1):
        if mrna[i:i+len(seed_rc)] == seed_rc:
            matches.append(i)
    
    return matches

def evaluate_binding_site(mrna: str, mirna: str, position: int, 
                         flank_size: int = 15) -> Dict:
    """
    Evaluate a potential miRNA binding site in the mRNA.
    
    Args:
        mrna: Target mRNA sequence
        mirna: miRNA sequence
        position: Starting position of the potential binding site
        flank_size: Size of flanking regions to include in evaluation
        
    Returns:
        Dictionary with binding site evaluation metrics
    """
    # Extract site with flanking regions
    site_start = max(0, position - flank_size)
    site_end = min(len(mrna), position + len(mirna) + flank_size)
    site_with_flanks = mrna[site_start:site_end]
    
    # Simple scoring (placeholder)
    # In practice, would use ViennaRNA's RNAduplex or RNAcofold
    binding_score = 0
    mirna_rc = reverse_complement(mirna)
    
    # Count matching bases as a simple score
    for i in range(min(len(mirna), len(mrna) - position)):
        if position + i < len(mrna) and i < len(mirna_rc):
            if mrna[position + i] == mirna_rc[i]:
                binding_score += 1
    
    # Return evaluation metrics
    return {
        "position": position,
        "site_with_flanks": site_with_flanks,
        "binding_score": binding_score,
        "normalized_score": binding_score / len(mirna) if len(mirna) > 0 else 0
    }

def design_binding_site(mirna: str, 
                       perfect_complement: bool = False,
                       include_flanking: bool = True,
                       flank_sequence: str = "NNNN") -> str:
    """
    Design an optimal binding site for a miRNA.
    
    Args:
        mirna: miRNA sequence
        perfect_complement: If True, create a perfect complement site
        include_flanking: If True, add flanking sequences
        flank_sequence: Sequence to use for flanking regions
        
    Returns:
        Designed binding site sequence
    """
    if perfect_complement:
        # Create perfect complement
        binding_site = reverse_complement(mirna)
    else:
        # Create seed-based binding site (positions 2-8)
        seed = get_seed_sequence(mirna)
        binding_site = reverse_complement(seed)
        
        # Add random matches for non-seed region
        non_seed = reverse_complement(mirna[8:])
        # In practice, optimize this region based on accessibility
        binding_site += non_seed
    
    # Add flanking regions if requested
    if include_flanking:
        binding_site = flank_sequence + binding_site + flank_sequence
        
    return binding_site

def validate_rna_sequence(sequence: str) -> bool:
    """
    Validate that a string is a valid RNA sequence.
    
    Args:
        sequence: String to validate
        
    Returns:
        True if the sequence is valid RNA, False otherwise
    """
    if not sequence:
        return False
    
    # Check for valid RNA characters (allow lowercase too)
    pattern = re.compile(r'^[AUGCaugc]+$')
    return bool(pattern.match(sequence)) 