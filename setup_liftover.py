#!/usr/bin/env python3
"""
Setup script for liftover chain files in the workflows/misc directory.
This script downloads and organizes chain files for the Snakemake liftover pipeline.
"""

import os
import shutil
import urllib.request
from pathlib import Path

def download_file(url: str, output_path: str):
    """Download a file from URL to output path."""
    print(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"✓ Downloaded: {output_path}")
        return True
    except Exception as e:
        print(f"✗ Failed to download {url}: {e}")
        return False

def copy_existing_chain_files():
    """Copy existing chain files from misc/ to workflows/misc/"""
    misc_dir = Path("misc")
    workflows_misc_dir = Path("workflows/misc")
    workflows_misc_dir.mkdir(exist_ok=True)
    
    if misc_dir.exists():
        chain_files = list(misc_dir.glob("*.chain.gz"))
        if chain_files:
            print(f"Copying {len(chain_files)} chain files from misc/ to workflows/misc/")
            for chain_file in chain_files:
                dest = workflows_misc_dir / chain_file.name
                shutil.copy2(chain_file, dest)
                print(f"✓ Copied: {chain_file.name}")
            return True
    return False

def main():
    print("Setting up liftover chain files for Snakemake pipeline")
    print("=" * 55)
    
    # Create workflows/misc directory
    workflows_misc_dir = Path("workflows/misc")
    workflows_misc_dir.mkdir(exist_ok=True)
    print(f"✓ Created directory: {workflows_misc_dir}")
    
    # Try to copy existing files first
    copied_existing = copy_existing_chain_files()
    
    # Chain file URLs and paths
    chain_files = {
        "hs1ToHg19.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hs1/liftOver/hs1ToHg19.over.chain.gz",
        "hs1ToHg38.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hs1/liftOver/hs1ToHg38.over.chain.gz",
        "hg19ToHg38.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "hg38ToHg19.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    }
    
    # Download missing chain files
    for filename, url in chain_files.items():
        output_path = workflows_misc_dir / filename
        if output_path.exists():
            print(f"✓ Already exists: {filename}")
        else:
            success = download_file(url, str(output_path))
            if not success:
                print(f"✗ Failed to download {filename}")
    
    # Check for T2T chain files
    t2t_files = ["hg38ToHs1.over.chain.gz", "hs1ToHg38.over.chain.gz"]
    missing_t2t = []
    
    for t2t_file in t2t_files:
        if not (workflows_misc_dir / t2t_file).exists():
            missing_t2t.append(t2t_file)
    
    if missing_t2t:
        print(f"\n⚠️  Missing T2T-CHM13 chain files:")
        for file in missing_t2t:
            print(f"   - {file}")
        print("\nDownload T2T chain files from:")
        print("https://github.com/marbl/CHM13-liftover")
        print("https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain_files/")
    
    # List final contents
    print(f"\nFinal chain files in {workflows_misc_dir}:")
    chain_files = list(workflows_misc_dir.glob("*.chain.gz"))
    for cf in sorted(chain_files):
        size_mb = cf.stat().st_size / (1024 * 1024)
        print(f"  ✓ {cf.name} ({size_mb:.1f} MB)")
    
    if not chain_files:
        print("  (No chain files found)")
        
    print(f"\nSetup complete! Chain files are in: {workflows_misc_dir.absolute()}")
    
    # Instructions for usage
    print("\nUsage:")
    print("1. To enable liftover in Snakemake, set 'liftover_to_hs1: true' in config.yaml")
    print("2. Run full pipeline with liftover: snakemake --cores 1")
    print("3. Run only liftover: snakemake liftover_all --cores 1")

if __name__ == "__main__":
    main()