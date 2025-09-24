# Snakemake workflow for RFD (Replication Fork Directionality) analysis
# 
# This workflow processes OK-seq bedgraph files to identify:
# - Initiation Sites (IS): local maxima
# - Termination Sites (TS): local minima  
# - Initiation Zones (IZ): regions from TS to IS
# - Termination Zones (TZ): regions from IS to TS

import os
import glob

# Configuration
configfile: "config.yaml"

# Directory configuration
INPUTS_DIR = "workflows/inputs"      # Main input directory for real data
EXAMPLES_DIR = "workflows/examples"  # Example data directory
RESULTS_DIR = "workflows/results"
SCRIPTS_DIR = "scripts"

# Find all .bedgraph files in both inputs and examples directories
# Prioritize inputs folder, fall back to examples if inputs is empty
inputs_files = glob.glob(f"{INPUTS_DIR}/*.bedgraph")
examples_files = glob.glob(f"{EXAMPLES_DIR}/*.bedgraph")

if inputs_files:
    # Use real data from inputs folder
    BEDGRAPH_FILES = inputs_files
    INPUT_DIR = INPUTS_DIR
    print(f"Using {len(inputs_files)} input files from {INPUTS_DIR}/")
else:
    # Fall back to example data
    BEDGRAPH_FILES = examples_files
    INPUT_DIR = EXAMPLES_DIR
    print(f"No files in {INPUTS_DIR}/, using {len(examples_files)} example files from {EXAMPLES_DIR}/")

SAMPLES = [os.path.basename(f).replace(".bedgraph", "") for f in BEDGRAPH_FILES]

if not SAMPLES:
    print(f"Warning: No .bedgraph files found in {INPUTS_DIR}/ or {EXAMPLES_DIR}/")
    print("Please add your bedgraph files to one of these directories.")
    SAMPLES = []

# Default rule - run all analyses
rule all:
    input:
        # Critical points (IS/TS) files
        expand(f"{RESULTS_DIR}/{{sample}}_critical_points.bed", sample=SAMPLES),
        # Zone classification files
        expand(f"{RESULTS_DIR}/{{sample}}_zones.bed", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv", sample=SAMPLES),
        # Summary reports
        expand(f"{RESULTS_DIR}/{{sample}}_summary.txt", sample=SAMPLES),
        # Plots (if enabled)
        expand(f"{RESULTS_DIR}/{{sample}}_plot.png", sample=SAMPLES) if config.get("generate_plots", True) else [],
        # Combined summary report
        f"{RESULTS_DIR}/combined_analysis_summary.tsv"

# Rule to run RFD analysis on each bedgraph file
rule rfd_analysis:
    input:
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph"
    output:
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/{{sample}}_zones.bed", 
        zones_detailed = f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv"
    params:
        output_prefix = f"{RESULTS_DIR}/{{sample}}",
        min_distance = config.get("min_peak_distance", 10),
        min_prominence = config.get("min_peak_prominence", 0.1),
        chromosome = config.get("chromosome", ""),
        plot_flag = "--plot" if config.get("generate_plots", True) else ""
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_analysis.log"
    shell:
        """
        python {SCRIPTS_DIR}/rfd_analyzer.py {input.bedgraph} \
            -o {params.output_prefix} \
            -d {params.min_distance} \
            -p {params.min_prominence} \
            {params.chromosome} \
            {params.plot_flag} \
            > {log} 2>&1
        """

# Rule to generate individual summary reports
rule generate_summary:
    input:
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph",
        zones_detailed = f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv"
    output:
        summary = f"{RESULTS_DIR}/{{sample}}_summary.txt"
    params:
        sample_name = "{sample}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_summary.log"
    run:
        import pandas as pd
        import numpy as np
        
        # Read zones data
        zones_df = pd.read_csv(input.zones_detailed, sep='\t')
        
        # Read original bedgraph for basic stats
        bedgraph_df = pd.read_csv(input.bedgraph, sep='\t', 
                                names=['chr', 'start', 'end', 'RFD'])
        
        # Calculate statistics
        iz_zones = zones_df[zones_df['zone_type'] == 'IZ']
        tz_zones = zones_df[zones_df['zone_type'] == 'TZ']
        
        # Write summary
        with open(output.summary, 'w') as f:
            f.write(f"RFD Analysis Summary: {params.sample_name}\n")
            f.write("=" * 50 + "\n\n")
            
            # Basic data info
            f.write("Dataset Information:\n")
            f.write(f"  Total genomic bins: {len(bedgraph_df):,}\n")
            f.write(f"  RFD value range: {bedgraph_df['RFD'].min():.3f} to {bedgraph_df['RFD'].max():.3f}\n")
            f.write(f"  Chromosomes: {', '.join(sorted(bedgraph_df['chr'].unique()))}\n\n")
            
            # Zone statistics
            f.write("Zone Analysis:\n")
            f.write(f"  Total zones identified: {len(zones_df)}\n")
            f.write(f"  Initiation Zones (IZ): {len(iz_zones)}\n")
            f.write(f"  Termination Zones (TZ): {len(tz_zones)}\n\n")
            
            if len(iz_zones) > 0:
                f.write(f"  IZ Statistics:\n")
                f.write(f"    Mean length: {iz_zones['length'].mean():,.0f} bp\n")
                f.write(f"    Median length: {iz_zones['length'].median():,.0f} bp\n")
                f.write(f"    Total coverage: {iz_zones['length'].sum():,.0f} bp\n\n")
            
            if len(tz_zones) > 0:
                f.write(f"  TZ Statistics:\n")
                f.write(f"    Mean length: {tz_zones['length'].mean():,.0f} bp\n")
                f.write(f"    Median length: {tz_zones['length'].median():,.0f} bp\n")
                f.write(f"    Total coverage: {tz_zones['length'].sum():,.0f} bp\n\n")

# Rule to generate plots (conditional)
rule generate_plot:
    input:
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph",
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed",
        zones_detailed = f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv"
    output:
        plot = f"{RESULTS_DIR}/{{sample}}_plot.png"
    params:
        sample_name = "{sample}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_plot.log"
    run:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        import pandas as pd
        import numpy as np
        
        # Read data
        bedgraph_df = pd.read_csv(input.bedgraph, sep='\t', 
                                names=['chr', 'start', 'end', 'RFD'])
        critical_points_df = pd.read_csv(input.critical_points, sep='\t',
                                       names=['chr', 'start', 'end', 'type', 'value'])
        zones_df = pd.read_csv(input.zones_detailed, sep='\t')
        
        # Create plot for first chromosome
        first_chr = bedgraph_df['chr'].iloc[0]
        chr_data = bedgraph_df[bedgraph_df['chr'] == first_chr]
        chr_critical = critical_points_df[critical_points_df['chr'] == first_chr]
        chr_zones = zones_df[zones_df['chr'] == first_chr]
        
        plt.figure(figsize=(15, 6))
        
        # Plot RFD signal
        plt.plot(chr_data['start'], chr_data['RFD'], 'b-', linewidth=1, alpha=0.7, label='RFD')
        
        # Plot critical points
        is_points = chr_critical[chr_critical['type'] == 'IS']
        ts_points = chr_critical[chr_critical['type'] == 'TS']
        
        if len(is_points) > 0:
            plt.scatter(is_points['start'], is_points['value'], 
                       color='red', s=50, marker='^', 
                       label=f'IS (n={len(is_points)})', zorder=5)
        
        if len(ts_points) > 0:
            plt.scatter(ts_points['start'], ts_points['value'], 
                       color='blue', s=50, marker='v', 
                       label=f'TS (n={len(ts_points)})', zorder=5)
        
        # Add zone shading
        for _, zone in chr_zones.iterrows():
            if zone['zone_type'] == 'IZ':
                plt.axvspan(zone['start'], zone['end'], alpha=0.2, color='green')
            elif zone['zone_type'] == 'TZ':
                plt.axvspan(zone['start'], zone['end'], alpha=0.2, color='orange')
        
        plt.xlabel('Genomic Position')
        plt.ylabel('RFD Value')
        plt.title(f'RFD Analysis - {params.sample_name} ({first_chr})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output.plot, dpi=300, bbox_inches='tight')
        plt.close()

# Rule to combine all summaries into one report
rule combine_summaries:
    input:
        summaries = expand(f"{RESULTS_DIR}/{{sample}}_summary.txt", sample=SAMPLES),
        zones_files = expand(f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv", sample=SAMPLES)
    output:
        combined = f"{RESULTS_DIR}/combined_analysis_summary.tsv"
    log:
        f"{RESULTS_DIR}/logs/combined_summary.log"
    run:
        import pandas as pd
        
        combined_data = []
        
        for i, zones_file in enumerate(input.zones_files):
            sample_name = SAMPLES[i]
            
            # Read zones data
            zones_df = pd.read_csv(zones_file, sep='\t')
            
            # Calculate basic statistics
            iz_count = len(zones_df[zones_df['zone_type'] == 'IZ'])
            tz_count = len(zones_df[zones_df['zone_type'] == 'TZ'])
            total_zones = len(zones_df)
            
            iz_zones = zones_df[zones_df['zone_type'] == 'IZ']
            tz_zones = zones_df[zones_df['zone_type'] == 'TZ']
            
            iz_mean_length = iz_zones['length'].mean() if len(iz_zones) > 0 else 0
            tz_mean_length = tz_zones['length'].mean() if len(tz_zones) > 0 else 0
            iz_total_coverage = iz_zones['length'].sum() if len(iz_zones) > 0 else 0
            tz_total_coverage = tz_zones['length'].sum() if len(tz_zones) > 0 else 0
            
            combined_data.append({
                'sample': sample_name,
                'total_zones': total_zones,
                'iz_count': iz_count,
                'tz_count': tz_count,
                'iz_mean_length': iz_mean_length,
                'tz_mean_length': tz_mean_length,
                'iz_total_coverage': iz_total_coverage,
                'tz_total_coverage': tz_total_coverage,
                'chromosomes': len(zones_df['chr'].unique())
            })
        
        # Create combined DataFrame and save
        combined_df = pd.DataFrame(combined_data)
        combined_df.to_csv(output.combined, sep='\t', index=False)

# Rule to create logs directory
rule create_logs_dir:
    output:
        directory(f"{RESULTS_DIR}/logs")
    shell:
        "mkdir -p {output}"

# Rule to copy example files to inputs directory for testing
rule copy_examples_to_inputs:
    input:
        expand(f"{EXAMPLES_DIR}/{{example}}.bedgraph", 
               example=[os.path.basename(f).replace(".bedgraph", "") 
                       for f in glob.glob(f"{EXAMPLES_DIR}/*.bedgraph")])
    output:
        expand(f"{INPUTS_DIR}/{{example}}.bedgraph",
               example=[os.path.basename(f).replace(".bedgraph", "") 
                       for f in glob.glob(f"{EXAMPLES_DIR}/*.bedgraph")])
    shell:
        """
        mkdir -p {INPUTS_DIR}
        cp {EXAMPLES_DIR}/*.bedgraph {INPUTS_DIR}/
        echo "Copied example files to inputs directory"
        """

# Make sure logs directory exists before running analyses
ruleorder: create_logs_dir > rfd_analysis