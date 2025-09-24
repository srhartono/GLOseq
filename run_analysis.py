#!/usr/bin/env python3
"""
Wrapper script for RFD Analysis Snakemake workflow with sensible defaults
Usage: python run_analysis.py [additional snakemake arguments]
"""

import sys
import subprocess
import os
from pathlib import Path

def main():
    # Get the directory containing this script
    parent_dir = Path(__file__).parent.absolute()
    # print(parent_dir)
    # Default arguments
    default_args = [
        "snakemake",
        "--cores", "1",  # Default to single core
        "--dag",  # Print the DAG of jobs
        "--use-conda",   # Use conda environments if available
        "--printshellcmds"  # Print shell commands for transparency
    ]
    
    # Add any additional arguments passed by user
    user_args = sys.argv[1:]
    
    # Check if user specified cores, if so remove our default
    cores_specified = any("--cores" in arg or "-c" in arg for arg in user_args)
    if cores_specified:
        default_args = [arg for arg in default_args if arg not in ["--cores", "1"]]
    
    # Combine arguments
    all_args = default_args + user_args
    
    # Print command for transparency
    print(f"Running: {' '.join(all_args)}")
    print("-" * 50)
    
    # Execute snakemake
    try:
        # Change to script directory to ensure Snakefile is found
        os.chdir(parent_dir)
        dag_filename = parent_dir/"dag.dot"
        with open(dag_filename, "w") as dag_file:
            result = subprocess.run(all_args, check=True, stdout=dag_file, stderr=subprocess.PIPE, text=True)

        from scripts import dot2pdf
        print(dag_filename)
        dot2pdf.main(dag_filename)
        
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error: Snakemake failed with exit code {e.returncode}")
        return e.returncode
    except FileNotFoundError:
        print("Error: snakemake not found. Please ensure it's installed and in your PATH.")
        print("Install with: pip install snakemake")
        return 1

if __name__ == "__main__":
    sys.exit(main())