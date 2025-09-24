#!/usr/bin/env python3
"""
Helper script to copy example files to inputs directory for testing
"""

import shutil
import glob
import os
from pathlib import Path

def copy_examples_to_inputs():
    """Copy all example bedgraph files to inputs directory"""
    
    script_dir = Path(__file__).parent
    examples_dir = script_dir / "workflows" / "examples"
    inputs_dir = script_dir / "workflows" / "inputs"
    
    # Create inputs directory if it doesn't exist
    inputs_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bedgraph files in examples
    example_files = list(examples_dir.glob("*.bedgraph"))
    
    if not example_files:
        print(f"No .bedgraph files found in {examples_dir}")
        return False
    
    # Copy files
    copied_files = []
    for example_file in example_files:
        target_file = inputs_dir / example_file.name
        shutil.copy2(example_file, target_file)
        copied_files.append(example_file.name)
    
    print(f"Copied {len(copied_files)} example files to workflows/inputs/:")
    for filename in copied_files:
        print(f"  - {filename}")
    
    print(f"\nNow you can run: python run_analysis.py")
    return True

if __name__ == "__main__":
    copy_examples_to_inputs()