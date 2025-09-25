#!/usr/bin/env python3
"""
Example usage of the liftover_annotations.py script.
This demonstrates various ways to use the liftover functionality.
"""

import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd):
    """Run a command and print it."""
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print("✓ Success!")
        if result.stdout.strip():
            print(result.stdout)
    else:
        print("✗ Failed!")
        if result.stderr.strip():
            print(result.stderr)
    print("-" * 50)
    return result.returncode == 0

def main():
    print("GLOseq Liftover Examples")
    print("=" * 50)
    
    # Check if pyliftover is installed
    try:
        import pyliftover
        print("✓ pyliftover is installed")
    except ImportError:
        print("✗ pyliftover not installed. Run: pip install pyliftover")
        return
    
    # Check for chain files
    chain_dir = Path("misc")
    if not chain_dir.exists():
        print("✗ misc directory not found. Run scripts/download_chain_files.py first")
        return
    
    chain_files = list(chain_dir.glob("*.chain.gz"))
    if not chain_files:
        print("✗ No chain files found. Run scripts/download_chain_files.py first")
        return
    
    print(f"✓ Found {len(chain_files)} chain files:")
    for cf in chain_files:
        print(f"  - {cf.name}")
    print()
    
    # Example BED files (use results from GLOseq analysis)
    results_dir = Path("workflows/results")
    bed_files = list(results_dir.glob("*_critical_points.bed")) + list(results_dir.glob("*_zones.bed"))
    
    if not bed_files:
        print("✗ No BED files found in workflows/results/")
        print("Run the GLOseq analysis first to generate BED files")
        return
    
    print(f"✓ Found {len(bed_files)} BED files for liftover:")
    for bf in bed_files[:3]:  # Show first 3
        print(f"  - {bf.name}")
    if len(bed_files) > 3:
        print(f"  ... and {len(bed_files) - 3} more")
    print()
    
    # Example 1: Single file liftover
    print("Example 1: Single file liftover (hg38 → hg19)")
    input_file = bed_files[0]
    output_file = input_file.parent / f"{input_file.stem}_hg19.bed"
    
    cmd = [
        "python", "scripts/liftover_annotations.py",
        str(input_file),
        "--from", "hg38",
        "--to", "hg19", 
        "--output", str(output_file)
    ]
    
    success = run_command(cmd)
    
    if success and output_file.exists():
        print(f"✓ Output saved to: {output_file}")
        # Show first few lines
        with open(output_file) as f:
            lines = f.readlines()[:3]
            for line in lines:
                print(f"  {line.strip()}")
    print()
    
    # Example 2: Batch processing
    print("Example 2: Batch processing (hg38 → hg19)")
    output_dir = results_dir / "hg19_lifted"
    
    cmd = [
        "python", "scripts/liftover_annotations.py"
    ] + [str(bf) for bf in bed_files[:2]]  + [
        "--from", "hg38",
        "--to", "hg19",
        "--output-dir", str(output_dir),
        "--suffix", "_hg19"
    ]
    
    success = run_command(cmd)
    
    if success:
        output_files = list(output_dir.glob("*_hg19.bed"))
        print(f"✓ {len(output_files)} files created in {output_dir}")
        for of in output_files:
            print(f"  - {of.name}")
    print()
    
    # Example 3: Using custom chain file
    if len(chain_files) > 0:
        print("Example 3: Using custom chain file")
        chain_file = chain_files[0]
        
        # Determine assemblies from chain file name
        if "hg19ToHg38" in chain_file.name:
            from_asm, to_asm = "hg19", "hg38"
        elif "hg38ToHg19" in chain_file.name:
            from_asm, to_asm = "hg38", "hg19"
        else:
            from_asm, to_asm = "hg38", "hg19"  # default
        
        output_file = input_file.parent / f"{input_file.stem}_custom.bed"
        
        cmd = [
            "python", "scripts/liftover_annotations.py",
            str(input_file),
            "--from", from_asm,
            "--to", to_asm,
            "--chain-file", str(chain_file),
            "--output", str(output_file)
        ]
        
        success = run_command(cmd)
        
        if success and output_file.exists():
            print(f"✓ Custom liftover completed: {output_file}")
    
    print("\nLiftover examples completed!")
    print("\nFor more options, run:")
    print("python scripts/liftover_annotations.py --help")

if __name__ == "__main__":
    main()