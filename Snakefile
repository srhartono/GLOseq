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
print("Running Snakefile", file=sys.stderr)
# Directory configuration
INPUTS_DIR = "workflows/inputs"      # Main input directory for real data
EXAMPLES_DIR = "workflows/examples"  # Example data directory
RESULTS_DIR = "workflows/results"
SCRIPTS_DIR = "scripts"

# Find all .bedgraph.gz files in both inputs and examples directories
# Prioritize inputs folder, fall back to examples if inputs is empty
inputs_files = glob.glob(f"{INPUTS_DIR}/*.bedgraph.gz")
examples_files = glob.glob(f"{EXAMPLES_DIR}/*.bedgraph.gz")

if inputs_files:
    # Use real data from inputs folder
    BEDGRAPH_FILES = inputs_files
    INPUT_DIR = INPUTS_DIR
    print(f"Using {len(inputs_files)} input files from {INPUTS_DIR}/",file=sys.stderr)
else:
    # Fall back to example data
    BEDGRAPH_FILES = examples_files
    INPUT_DIR = EXAMPLES_DIR
    print(f"No files in {INPUT_DIR}/, using {len(examples_files)} example files from {EXAMPLES_DIR}/", file=sys.stderr)

SAMPLES = [os.path.basename(f).replace(".bedgraph.gz", "") for f in BEDGRAPH_FILES]

if not SAMPLES:
    print(f"Warning: No .bedgraph files found in {INPUT_DIR}/ or {EXAMPLES_DIR}/", file=sys.stderr)
    print("Please add your bedgraph files to one of these directories.", file=sys.stderr)
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
        f"{RESULTS_DIR}/combined_analysis_summary.tsv",
        # Liftover results (if enabled)
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_critical_points.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else [],
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_zones.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else []

# Rule to run RFD analysis on each bedgraph.gz file
rule rfd_analysis:
    input:
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph.gz"
    output:
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/{{sample}}_zones.bed", 
        zones_detailed = f"{RESULTS_DIR}/{{sample}}_zones_detailed.tsv"
    params:
        output_prefix = f"{RESULTS_DIR}/{{sample}}",
        min_distance = config.get("min_peak_distance", 100),
        min_prominence = config.get("min_peak_prominence", 1),
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
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph.gz",
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
        import gzip
        with gzip.open(input.bedgraph, 'rt') as f:
            bedgraph_df = pd.read_csv(f, sep='\t', 
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
        bedgraph = f"{INPUT_DIR}/{{sample}}.bedgraph.gz",
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
        import gzip
        with gzip.open(input.bedgraph, 'rt') as f:
            bedgraph_df = pd.read_csv(f, sep='\t', 
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
        
        plt.xlabel('Genomic Poszition')
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

# Liftover rules for genome assembly conversion
# These rules convert coordinates from hg19 -> hg38 -> hs1 (T2T-CHM13)

# Rule to liftover from original assembly (hg19) to hg38
rule liftover_hg19_to_hg38:
    input:
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/{{sample}}_zones.bed"
    output:
        critical_points = f"{RESULTS_DIR}/hg38/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/hg38/{{sample}}_zones.bed"
    params:
        chain_file = "workflows/misc/hg19ToHg38.over.chain.gz",
        sample = "{sample}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_liftover_hg19_to_hg38.log"
    run:
        import os
        import subprocess
        
        # Create output directory
        os.makedirs(os.path.dirname(output.critical_points), exist_ok=True)
        
        # Liftover critical points
        cmd1 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.critical_points,
            "--from", "hg19",
            "--to", "hg38", 
            "--chain-file", params.chain_file,
            "--output", output.critical_points
        ]
        
        with open(log[0], 'w') as logfile:
            result1 = subprocess.run(cmd1, capture_output=True, text=True)
            logfile.write(f"Critical points liftover:\n{result1.stdout}\n{result1.stderr}\n\n")
            if result1.returncode != 0:
                raise Exception(f"Critical points liftover failed: {result1.stderr}")
        
        # Liftover zones
        cmd2 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.zones,
            "--from", "hg19",
            "--to", "hg38",
            "--chain-file", params.chain_file,
            "--output", output.zones
        ]
        
        with open(log[0], 'a') as logfile:
            result2 = subprocess.run(cmd2, capture_output=True, text=True)
            logfile.write(f"Zones liftover:\n{result2.stdout}\n{result2.stderr}\n")
            if result2.returncode != 0:
                raise Exception(f"Zones liftover failed: {result2.stderr}")

# Rule to liftover from hg38 to hs1 (T2T-CHM13)
rule liftover_hg38_to_hs1:
    input:
        critical_points = f"{RESULTS_DIR}/hg38/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/hg38/{{sample}}_zones.bed"
    output:
        critical_points = f"{RESULTS_DIR}/hs1/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/hs1/{{sample}}_zones.bed"
    params:
        chain_file = "workflows/misc/hg38ToHs1.over.chain.gz",
        sample = "{sample}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_liftover_hg38_to_hs1.log"
    run:
        import os
        import subprocess
        
        # Create output directory
        os.makedirs(os.path.dirname(output.critical_points), exist_ok=True)
        
        # Liftover critical points
        cmd1 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.critical_points,
            "--from", "hg38",
            "--to", "hs1",
            "--chain-file", params.chain_file,
            "--output", output.critical_points
        ]
        
        with open(log[0], 'w') as logfile:
            result1 = subprocess.run(cmd1, capture_output=True, text=True)
            logfile.write(f"Critical points liftover:\n{result1.stdout}\n{result1.stderr}\n\n")
            if result1.returncode != 0:
                raise Exception(f"Critical points liftover failed: {result1.stderr}")
        
        # Liftover zones
        cmd2 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.zones,
            "--from", "hg38",
            "--to", "hs1",
            "--chain-file", params.chain_file,
            "--output", output.zones
        ]
        
        with open(log[0], 'a') as logfile:
            result2 = subprocess.run(cmd2, capture_output=True, text=True)
            logfile.write(f"Zones liftover:\n{result2.stdout}\n{result2.stderr}\n")
            if result2.returncode != 0:
                raise Exception(f"Zones liftover failed: {result2.stderr}")

# Combined rule for complete liftover pipeline (hg19 -> hg38 -> hs1)
rule liftover_to_hs1:
    input:
        critical_points = f"{RESULTS_DIR}/hs1/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/hs1/{{sample}}_zones.bed"
    output:
        touch(f"{RESULTS_DIR}/{{sample}}_liftover_complete.flag")
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_liftover_complete.log"
    shell:
        """
        echo "Liftover pipeline complete for sample {wildcards.sample}" > {log}
        echo "Files created:" >> {log}
        echo "  - {input.critical_points}" >> {log} 
        echo "  - {input.zones}" >> {log}
        """

# Rule to create logs directory
rule create_logs_dir:
    output:
        directory(f"{RESULTS_DIR}/logs")
    shell:
        "mkdir -p {output}"

# TRC Liftover rules - liftover TRC results from hg19 to hg38 to hs1
# These rules ensure that TRC analysis is done in hg19, then results are lifted over

# Rule to liftover TRC results from hg19 to hg38
rule trc_liftover_hg19_to_hg38:
    input:
        head_on_genes = f"{RESULTS_DIR}/{{sample}}_head_on_genes.bed",
        co_directional_genes = f"{RESULTS_DIR}/{{sample}}_co_directional_genes.bed"
    output:
        head_on_genes = f"{RESULTS_DIR}/hg38/{{sample}}_head_on_genes.bed",
        co_directional_genes = f"{RESULTS_DIR}/hg38/{{sample}}_co_directional_genes.bed"
    params:
        chain_file = "workflows/misc/hg19ToHg38.over.chain.gz"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_trc_liftover_hg19_to_hg38.log"
    run:
        import os
        import subprocess
        
        # Create output directory
        os.makedirs(os.path.dirname(output.head_on_genes), exist_ok=True)
        
        # Liftover head-on genes
        cmd1 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.head_on_genes,
            "--from", "hg19",
            "--to", "hg38", 
            "--chain-file", params.chain_file,
            "--output", output.head_on_genes
        ]
        
        with open(log[0], 'w') as logfile:
            result1 = subprocess.run(cmd1, capture_output=True, text=True)
            logfile.write(f"Head-on genes liftover:\n{result1.stdout}\n{result1.stderr}\n\n")
            if result1.returncode != 0:
                raise Exception(f"Head-on genes liftover failed: {result1.stderr}")
        
        # Liftover co-directional genes
        cmd2 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.co_directional_genes,
            "--from", "hg19",
            "--to", "hg38",
            "--chain-file", params.chain_file,
            "--output", output.co_directional_genes
        ]
        
        with open(log[0], 'a') as logfile:
            result2 = subprocess.run(cmd2, capture_output=True, text=True)
            logfile.write(f"Co-directional genes liftover:\n{result2.stdout}\n{result2.stderr}\n")
            if result2.returncode != 0:
                raise Exception(f"Co-directional genes liftover failed: {result2.stderr}")

# Rule to liftover TRC results from hg38 to hs1
rule trc_liftover_hg38_to_hs1:
    input:
        head_on_genes = f"{RESULTS_DIR}/hg38/{{sample}}_head_on_genes.bed",
        co_directional_genes = f"{RESULTS_DIR}/hg38/{{sample}}_co_directional_genes.bed"
    output:
        head_on_genes = f"{RESULTS_DIR}/hs1/{{sample}}_head_on_genes.bed",
        co_directional_genes = f"{RESULTS_DIR}/hs1/{{sample}}_co_directional_genes.bed"
    params:
        chain_file = "workflows/misc/hg38ToHs1.over.chain.gz"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_trc_liftover_hg38_to_hs1.log"
    run:
        import os
        import subprocess
        
        # Create output directory
        os.makedirs(os.path.dirname(output.head_on_genes), exist_ok=True)
        
        # Liftover head-on genes
        cmd1 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.head_on_genes,
            "--from", "hg38",
            "--to", "hs1",
            "--chain-file", params.chain_file,
            "--output", output.head_on_genes
        ]
        
        with open(log[0], 'w') as logfile:
            result1 = subprocess.run(cmd1, capture_output=True, text=True)
            logfile.write(f"Head-on genes liftover:\n{result1.stdout}\n{result1.stderr}\n\n")
            if result1.returncode != 0:
                raise Exception(f"Head-on genes liftover failed: {result1.stderr}")
        
        # Liftover co-directional genes
        cmd2 = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe", f"{SCRIPTS_DIR}/liftover_annotations.py",
            input.co_directional_genes,
            "--from", "hg38",
            "--to", "hs1",
            "--chain-file", params.chain_file,
            "--output", output.co_directional_genes
        ]
        
        with open(log[0], 'a') as logfile:
            result2 = subprocess.run(cmd2, capture_output=True, text=True)
            logfile.write(f"Co-directional genes liftover:\n{result2.stdout}\n{result2.stderr}\n")
            if result2.returncode != 0:
                raise Exception(f"Co-directional genes liftover failed: {result2.stderr}")

# Rule to copy example files to inputs directory for testing
rule copy_examples_to_inputs:
    input:
        expand(f"{EXAMPLES_DIR}/{{example}}.bedgraph.gz", 
               example=[os.path.basename(f).replace(".bedgraph.gz", "") 
                       for f in glob.glob(f"{EXAMPLES_DIR}/*.bedgraph.gz")])
    output:
        expand(f"{INPUTS_DIR}/{{example}}.bedgraph.gz",
               example=[os.path.basename(f).replace(".bedgraph.gz", "") 
                       for f in glob.glob(f"{EXAMPLES_DIR}/*.bedgraph.gz")])
    shell:
        """
        mkdir -p {INPUTS_DIR}
        cp {EXAMPLES_DIR}/*.bedgraph.gz {INPUTS_DIR}/
        echo "Copied example files to inputs directory"
        """

# Target rule for liftover only (run with: snakemake liftover_all)
rule liftover_all:
    input:
        # Basic RFD results
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_critical_points.bed", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_zones.bed", sample=SAMPLES),
        # TRC results (if gene coordinates file is provided)
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_head_on_genes.bed", sample=SAMPLES) if config.get("gene_coordinates_file", "") else [],
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_co_directional_genes.bed", sample=SAMPLES) if config.get("gene_coordinates_file", "") else []
    message:
        "Completed liftover for all samples from hg19 -> hg38 -> hs1"

# Make sure logs directory exists before running analyses
ruleorder: create_logs_dir > rfd_analysis

# Rule for transcription-replication conflict analysis
rule trc_analysis:
    input:
        genes = config.get("gene_coordinates_file", ""),  # Gene coordinates file
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed",
        zones = f"{RESULTS_DIR}/{{sample}}_zones.bed",
        rfd_data = f"{INPUT_DIR}/{{sample}}.bedgraph.gz"
    output:
        trc_results = f"{RESULTS_DIR}/{{sample}}_trc_analysis.tsv",
        head_on_genes = f"{RESULTS_DIR}/{{sample}}_head_on_genes.bed",
        co_directional_genes = f"{RESULTS_DIR}/{{sample}}_co_directional_genes.bed"
    params:
        output_prefix = f"{RESULTS_DIR}/{{sample}}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_trc_analysis.log"
    run:
        import os
        import subprocess
        
        # Check if gene coordinates file is provided
        if not input.genes or not os.path.exists(input.genes):
            with open(log[0], 'w') as logfile:
                logfile.write("No gene coordinates file provided or file not found.\n")
                logfile.write("Set 'gene_coordinates_file' in config.yaml to enable TRC analysis.\n")
            
            # Create empty output files
            for output_file in output:
                with open(output_file, 'w') as f:
                    if output_file.endswith('.tsv'):
                        f.write("# No gene coordinates provided - TRC analysis skipped\n")
                    else:
                        f.write("# No gene coordinates provided - TRC analysis skipped\n")
            return
        
        # Run TRC analysis
        cmd = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe",
            f"{SCRIPTS_DIR}/trc_analyzer.py",
            "--genes", input.genes,
            "--rfd-data", input.rfd_data,
            "--critical-points", input.critical_points,
            "--zones", input.zones,
            "--output", params.output_prefix,
            "--formats", "tsv", "bed"
        ]
        
        with open(log[0], 'w') as logfile:
            result = subprocess.run(cmd, capture_output=True, text=True)
            logfile.write(f"TRC analysis command: {' '.join(cmd)}\n\n")
            logfile.write(f"STDOUT:\n{result.stdout}\n\n")
            logfile.write(f"STDERR:\n{result.stderr}\n")
            
            if result.returncode != 0:
                raise Exception(f"TRC analysis failed with return code {result.returncode}")

# Rule for TRC peak distribution analysis
rule trc_peak_analysis:
    input:
        trc_results = f"{RESULTS_DIR}/{{sample}}_trc_analysis.tsv",
        critical_points = f"{RESULTS_DIR}/{{sample}}_critical_points.bed"
    output:
        peak_report = f"{RESULTS_DIR}/{{sample}}_peak_analysis_report.txt",
        peak_plot = f"{RESULTS_DIR}/{{sample}}_peak_distributions.png",
        enhanced_results = f"{RESULTS_DIR}/{{sample}}_enhanced_trc_results.tsv"
    params:
        output_prefix = f"{RESULTS_DIR}/{{sample}}"
    log:
        f"{RESULTS_DIR}/logs/{{sample}}_peak_analysis.log"
    run:
        import subprocess
        import shutil
        import os
        
        # Use temporary prefix to avoid name conflicts
        temp_prefix = f"{RESULTS_DIR}/{wildcards.sample}_temp_peak_analysis"
        
        cmd = [
            "C:/cygwin64/home/srhar/GLOseq/.venv/Scripts/python.exe",
            f"{SCRIPTS_DIR}/trc_peak_analyzer.py",
            "--trc-results", input.trc_results,
            "--critical-points", input.critical_points,
            "--output", temp_prefix,
            "--create-plots"
        ]
        
        # Run command and capture output
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Write log
        with open(log[0], 'w') as logfile:
            logfile.write(f"Peak analysis command: {' '.join(cmd)}\n\n")
            logfile.write(f"STDOUT:\n{result.stdout}\n\n")
            logfile.write(f"STDERR:\n{result.stderr}\n")
            
            if result.returncode != 0:
                raise Exception(f"Peak analysis failed with return code {result.returncode}")
            
            # Move files to expected names
            temp_files = {
                f"{temp_prefix}_peak_analysis_report.txt": output.peak_report,
                f"{temp_prefix}_peak_distributions.png": output.peak_plot,
                f"{temp_prefix}_enhanced_trc_results.tsv": output.enhanced_results
            }
            
            for temp_file, final_file in temp_files.items():
                if os.path.exists(temp_file):
                    shutil.move(temp_file, final_file)
                    logfile.write(f"Moved {temp_file} to {final_file}\n")

# Rule to run TRC analysis for all samples (hg19 only)
rule trc_all:
    input:
        expand(f"{RESULTS_DIR}/{{sample}}_trc_analysis.tsv", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_head_on_genes.bed", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_co_directional_genes.bed", sample=SAMPLES)
    message:
        "Completed transcription-replication conflict analysis for all samples (hg19 coordinates)"

# Rule to run TRC analysis with liftover to all assemblies
rule trc_all_assemblies:
    input:
        # hg19 results (original analysis)
        expand(f"{RESULTS_DIR}/{{sample}}_trc_analysis.tsv", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_head_on_genes.bed", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_co_directional_genes.bed", sample=SAMPLES),
        # Lifted over results (if liftover is enabled)
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_head_on_genes.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else [],
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_co_directional_genes.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else []
    message:
        "Completed TRC analysis for all samples with liftover to all assemblies"

# Rule to run complete TRC analysis including peaks and liftover
rule trc_complete:
    input:
        # Peak analysis results
        expand(f"{RESULTS_DIR}/{{sample}}_peak_analysis_report.txt", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_peak_distributions.png", sample=SAMPLES),
        expand(f"{RESULTS_DIR}/{{sample}}_enhanced_trc_results.tsv", sample=SAMPLES),
        # TRC results with liftover (if enabled)
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_head_on_genes.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else [],
        expand(f"{RESULTS_DIR}/hs1/{{sample}}_co_directional_genes.bed", sample=SAMPLES) if config.get("liftover_to_hs1", False) else []
    message:
        "Completed comprehensive TRC analysis with peak distributions and liftover for all samples"

# Clean rule to remove all result files
rule clean:
    run:
        import os
        import shutil
        
        print("Cleaning all result files...")
        
        # Remove the entire results directory
        if os.path.exists(RESULTS_DIR):
            shutil.rmtree(RESULTS_DIR)
            print(f"Removed directory: {RESULTS_DIR}")
        
        # Remove any temporary files (skip .snakemake logs to avoid permission issues)
        temp_files = [
            "dag.svg", "dag.png", "dag.pdf",  # DAG visualization files
            "*.pyc", "__pycache__"  # Python cache files
        ]
        
        for pattern in temp_files:
            if "*" in pattern:
                import glob
                for file in glob.glob(pattern):
                    try:
                        if os.path.isfile(file):
                            os.remove(file)
                            print(f"Removed file: {file}")
                        elif os.path.isdir(file):
                            shutil.rmtree(file)
                            print(f"Removed directory: {file}")
                    except (OSError, PermissionError) as e:
                        print(f"Skipped {file}: {e}")
            else:
                try:
                    if os.path.exists(pattern):
                        if os.path.isfile(pattern):
                            os.remove(pattern)
                            print(f"Removed file: {pattern}")
                        elif os.path.isdir(pattern):
                            shutil.rmtree(pattern)
                            print(f"Removed directory: {pattern}")
                except (OSError, PermissionError) as e:
                    print(f"Skipped {pattern}: {e}")
        
        print("Clean completed successfully!")
        print("Note: Snakemake logs (.snakemake) are preserved to avoid permission issues")

# Clean rule for results only (keeps Snakemake metadata)
rule clean_results:
    run:
        import os
        import shutil
        
        print("Cleaning result files only...")
        
        # Remove the entire results directory
        if os.path.exists(RESULTS_DIR):
            shutil.rmtree(RESULTS_DIR)
            print(f"Removed directory: {RESULTS_DIR}")
        else:
            print(f"Directory {RESULTS_DIR} does not exist")
        
        print("Results cleanup completed!")

# Clean Snakemake logs and metadata (run this separately when no Snakemake is running)
rule clean_logs:
    run:
        import os
        import shutil
        
        print("Cleaning Snakemake logs and metadata...")
        
        # Remove Snakemake metadata
        if os.path.exists(".snakemake"):
            try:
                shutil.rmtree(".snakemake")
                print("Removed Snakemake metadata: .snakemake")
            except (OSError, PermissionError) as e:
                print(f"Could not remove .snakemake: {e}")
                print("Make sure no Snakemake processes are running and try again")
        else:
            print("No .snakemake directory found")
        
        print("Log cleanup completed!")

# Help rule to show available commands
rule help:
    run:
        print("""
GLOseq Snakemake Workflow - Available Commands:

MAIN ANALYSIS:
  snakemake all                    - Run complete RFD analysis pipeline
  snakemake --cores 4 all          - Run with 4 parallel cores

INDIVIDUAL STEPS:
  snakemake rfd_analysis           - Run RFD analysis only  
  snakemake plot_results           - Generate plots only
  snakemake summary_reports        - Generate summary reports only
  snakemake trc_all                - Run transcription-replication conflict analysis

LIFTOVER (Genome Assembly Conversion):
  snakemake liftover_all           - Convert coordinates: hg19 → hg38 → hs1
  snakemake liftover_hg19_to_hg38  - Convert hg19 → hg38 only
  snakemake liftover_hg38_to_hs1   - Convert hg38 → hs1 (T2T-CHM13) only

UTILITY:
  snakemake clean                  - Remove all results and temporary files
  snakemake clean_results          - Remove results directory only
  snakemake clean_logs             - Remove Snakemake logs (run when no processes active)
  snakemake help                   - Show this help message
  snakemake copy_examples_to_inputs- Copy example data to inputs folder

VISUALIZATION:
  snakemake --dag | dot -Tsvg > dag.svg     - Generate workflow DAG
  snakemake --dry-run              - Show planned jobs without execution

CONFIGURATION:
  - Edit config.yaml for analysis parameters
  - Place .bedgraph.gz files in workflows/inputs/
  - Chain files for liftover in workflows/misc/

For more details, see README.md and WORKFLOW_SUMMARY.md
        """)