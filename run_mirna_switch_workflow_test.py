#!/usr/bin/env python3
"""
Test script to verify the UTR sequence imports.
"""

import os
import sys

# Define base directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Import the UTR sequences
sys.path.append(os.path.join(script_dir, 'data'))

try:
    from utr_sequences.hbb_utr_sequences import (
        HBB_5UTR_RNA, HBB_3UTR_RNA, HBB_3UTR_RNA_SHORT,
        HBB_5UTR_DNA, HBB_3UTR_DNA, HBB_3UTR_DNA_SHORT,
        HBB_5UTR, HBB_3UTR, HBB_3UTR_SHORT  # Aliases
    )
    
    print("Successfully imported UTR sequences!")
    print(f"HBB_3UTR_RNA_SHORT length: {len(HBB_3UTR_RNA_SHORT)}")
    print(f"First 30 nucleotides: {HBB_3UTR_RNA_SHORT[:30]}")
    print(f"RNA version has U's: {'U' in HBB_3UTR_RNA_SHORT}")
    print(f"DNA version has T's: {'T' in HBB_3UTR_DNA_SHORT}")
    
    # Check alias mapping
    print("\nVerifying aliases:")
    print(f"HBB_3UTR_SHORT is HBB_3UTR_RNA_SHORT: {HBB_3UTR_SHORT is HBB_3UTR_RNA_SHORT}")
    
except ImportError as e:
    print(f"Error importing UTR sequences: {e}")
    print(f"Current path: {sys.path}")
    sys.exit(1)

print("\nAll imports successful!") 