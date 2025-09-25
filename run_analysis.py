#!/usr/bin/env python3
"""
Wrapper script for RFD Analysis Snakemake workflow with sensible defaults
Usage: python run_analysis.py [additional snakemake arguments]
"""

import sys
import subprocess
import os
from pathlib import Path
import argparse as args
def main():
    # Get the directory containing this script
    parent_dir = Path(__file__).parent.absolute()
    # print(parent_dir)
    # Default arguments
    # "--dag",  # Print the DAG of jobs
    default_args = [
        "snakemake",
        "--cores", "1",  # Default to single core
        # "--use-conda",   # Use conda environments if available
        "--printshellcmds"  # Print shell commands for transparency
    ]
    VERSION = (parent_dir/"VERSION.md").read_text().strip()
    # Filter out handled arguments from user_args to avoid duplication
    handled_args = ["--dryrun", "-0", "-c", "--clean", "-d", "--dag", "-t", "--trc", "-T", "--trc-complete", "-A", "--trc-all-assemblies", "-g", "--gene-file"]
    
    # More sophisticated filtering to handle arguments with values
    filtered_user_args = []
    skip_next = False
    for i, arg in enumerate(sys.argv[1:]):
        if skip_next:
            skip_next = False
            continue
        if any(arg.startswith(h) for h in handled_args):
            # Check if this argument takes a value
            if arg in ["-g", "--gene-file"] and i + 1 < len(sys.argv[1:]):
                skip_next = True  # Skip the next argument (the value)
            continue
        filtered_user_args.append(arg)
    
    user_args = filtered_user_args 
    myargs = args.ArgumentParser(description="RFD Analysis Snakemake Wrapper with TRC Analysis Options")
    myargs.add_argument("-j", "--cores", type=int, help="Number of CPU cores to use [1]")
    #myargs.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    myargs.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}", help="Version Information")
    myargs.add_argument("-m","--snakefile", type=str, default=str(parent_dir/"Snakefile"), help="Path to the Snakefile [default: %(default)s]")
    myargs.add_argument("-0","--dryrun", action="store_true", help="Perform a dry run without executing commands")
    myargs.add_argument("-d","--dag", action="store_true", help="Generate and save the DAG of jobs to dag.dot and dag.pdf")
    myargs.add_argument("-c", "--clean", action="store_true", help="Remove all results and temporary files")
    myargs.add_argument("-t", "--trc", action="store_true", help="Run transcription-replication conflict (TRC) analysis (hg19 only)")
    myargs.add_argument("-T", "--trc-complete", action="store_true", help="Run complete TRC analysis including peak distribution with liftover")
    myargs.add_argument("-A", "--trc-all-assemblies", action="store_true", help="Run TRC analysis with liftover to all assemblies")
    myargs.add_argument("-g", "--gene-file", type=str, help="Path to gene coordinates file for TRC analysis")

    myargs2 = myargs.parse_args()
    # Check if user specified cores, if so remove our default
    cores_specified = myargs2.cores #("--cores" in arg or "-c" in arg for arg in user_args)
    if cores_specified:
        default_args = [arg for arg in default_args if arg not in ["--cores", "1"]]
    
    # Combine arguments
    all_args = default_args + user_args
    # Handle special targets
    target_specified = False
    
    if (myargs2.dag):
        all_args += ["--dag"]
    if (myargs2.clean):
        all_args += ["clean"]
        target_specified = True
    # Set target based on TRC options
    target = None
    if (myargs2.trc):
        target = "trc_all"
        target_specified = True
    if (myargs2.trc_complete):
        target = "trc_complete"
        target_specified = True
    if (myargs2.trc_all_assemblies):
        target = "trc_all_assemblies"
        target_specified = True
        
    # Add target to arguments
    if target:
        all_args += [target]
    elif not target_specified and not myargs2.dag:
        all_args += ["all"]
        
    # Handle gene file override
    if myargs2.gene_file:
        # Update config dynamically by creating a config override
        all_args += ["--config", f"gene_coordinates_file={myargs2.gene_file}"]
    
    # Print analysis information
    if myargs2.trc:
        print("🧬 Running TRC Analysis (Transcription-Replication Conflicts)")
        print("   Analysis performed in hg19 coordinates only")
    elif myargs2.trc_complete:
        print("🧬 Running Complete TRC Analysis with Peak Distribution")
        print("   Analysis in hg19, with liftover to other assemblies if enabled")
    elif myargs2.trc_all_assemblies:
        print("🧬 Running TRC Analysis for All Genome Assemblies")
        print("   Analysis in hg19, with liftover to hg38 and hs1")
    elif myargs2.clean:
        print("🧹 Cleaning analysis results...")
    elif myargs2.dag:
        print("📊 Generating DAG visualization...")
    else:
        print("🔬 Running Standard RFD Analysis")
        
    if myargs2.gene_file:
        print(f"📋 Using custom gene coordinates: {myargs2.gene_file}")
    elif myargs2.trc or myargs2.trc_complete or myargs2.trc_all_assemblies:
        print(f"📋 Using default hg19 gene coordinates from config.yaml")
    
    # Print command for transparency
    print(f"Running: {' '.join(all_args)}")
    print("-" * 50)
    if myargs2.dryrun:
        print("Debug ran successfully. Exiting without executing snakemake.\n")
        sys.exit(0)
    # Execute snakemake
    try:
        # Change to script directory to ensure Snakefile is found
        os.chdir(parent_dir)
        dag_filename = parent_dir/"dag.dot"
        
        if myargs2.dag:
            with open(dag_filename, "w") as dag_file:
                result = subprocess.run(all_args, check=True, stdout=dag_file, stderr=subprocess.PIPE, text=True)
            print(f"Generating DAG visualization to dag.pdf...")
            from scripts import dot2pdf
            print(dag_filename)
            dot2pdf.main(dag_filename)
        else:
            result = subprocess.run(all_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(f"Snakemake output:\n{result.stdout}")
            if result.stderr:
                print(f"Snakemake errors:\n{result.stderr}")
            
            # Print completion message
            if myargs2.trc:
                print("\n✅ TRC Analysis completed!")
                print("📁 hg19 Results: workflows/results/*_trc_analysis.tsv")
                print("🎨 hg19 BED files: *_head_on_genes.bed, *_co_directional_genes.bed")
            elif myargs2.trc_all_assemblies:
                print("\n✅ TRC Analysis for All Assemblies completed!")
                print("📁 hg19 Results: workflows/results/*_trc_analysis.tsv")
                print("🎨 hg19 BED files: *_head_on_genes.bed, *_co_directional_genes.bed")
                print("🔄 hs1 Results: workflows/results/hs1/*_head_on_genes.bed, *_co_directional_genes.bed")
            elif myargs2.trc_complete:
                print("\n✅ Complete TRC Analysis with Peak Distribution completed!")
                print("📁 Results in workflows/results/")
                print("   hg19 results:")
                print("   - *_trc_analysis.tsv (basic TRC results in hg19)")
                print("   - *_enhanced_trc_results.tsv (with IS/TS distances)")
                print("   - *_peak_analysis_report.txt (statistical summary)")
                print("   - *_peak_distributions.png (visualization)")
                print("   - *_head_on_genes.bed, *_co_directional_genes.bed (strand-colored)")
                print("   Lifted over results in: workflows/results/hs1/ (if enabled)")
            elif not myargs2.clean and not myargs2.dag:
                print("\n✅ RFD Analysis completed!")
                print("📁 Results in workflows/results/")

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