#!/usr/bin/env python3
"""
Download chain files for genome liftover.
This script downloads commonly needed chain files for the GLOseq pipeline.
"""

import os
import urllib.request
import sys
from pathlib import Path

def download_file(url: str, output_path: str):
    """Download a file from URL to output path."""
    print(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"Downloaded: {output_path}")
        return True
    except Exception as e:
        print(f"Failed to download {url}: {e}")
        return False

def main():
    # Create misc directory
    misc_dir = Path("misc")
    misc_dir.mkdir(exist_ok=True)
    
    # Chain file URLs
    chain_files = {
        "hg19ToHg38.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "hg38ToHg19.over.chain.gz": "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz",
    }
    
    print("Downloading UCSC chain files...")
    for filename, url in chain_files.items():
        output_path = misc_dir / filename
        if output_path.exists():
            print(f"Already exists: {output_path}")
        else:
            success = download_file(url, str(output_path))
            if not success:
                print(f"Failed to download {filename}")
    
    print("\nNote: T2T-CHM13 (hs1) chain files need to be downloaded manually from:")
    print("https://github.com/marbl/CHM13-liftover")
    print("https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain_files/")
    
    print(f"\nChain files are saved in: {misc_dir.absolute()}")

if __name__ == "__main__":
    main()