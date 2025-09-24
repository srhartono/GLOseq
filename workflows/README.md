# Snakemake Workflow for RFD Analysis

This directory contains a complete Snakemake workflow for analyzing RFD (Replication Fork Directionality) data from OK-seq-like experiments (e.g. GLOe-seq).

## Directory Structure

```
GLOseq/
├── workflows/
│   ├── inputs/             # Your real bedgraph files (place files here)
│   ├── examples/           # Sample/example bedgraph files
│   └── results/            # Analysis outputs
├── scripts/                # Python analysis scripts
├── Snakefile              # Main workflow definition
├── config.yaml            # Configuration parameters
└── README.md              # Documentation
```

## Quick Start

1. **Place your bedgraph files** in `workflows/inputs/` (or use examples in `workflows/examples/`)
2. **Adjust parameters** in `config.yaml` if needed  
3. **Run the workflow** (choose one):
   ```bash
   # Easy way - using wrapper scripts (cores=1 by default)
   python run_analysis.py              # Cross-platform
   .\run_analysis.ps1                  # PowerShell (Windows)
   run_analysis.bat                    # Batch (Windows)
   
   # Direct snakemake (requires --cores specification)
   snakemake --cores 4
   ```

## Configuration

Edit `config.yaml` to customize analysis parameters:

- `min_peak_distance`: Minimum distance between peaks (default: 10 bins)
- `min_peak_prominence`: Minimum peak prominence (default: 0.1)  
- `generate_plots`: Create visualization plots (default: true)
- `analysis_mode`: Choose "standard", "sensitive", or "stringent"

## Input Format

Bedgraph files should be tab-separated with columns:
```
chr    start    end    RFD_value
chr1   1000000  1010000  0.245
chr1   1010000  1020000  0.312
...
```

## Workflow Steps

1. **RFD Analysis** (`rfd_analysis`)
   - Identifies initiation sites (IS) and termination sites (TS)
   - Classifies initiation zones (IZ) and termination zones (TZ)
   - Outputs: BED files with critical points and zones

2. **Summary Generation** (`generate_summary`)
   - Creates detailed statistics for each sample
   - Outputs: Text summary files

3. **Visualization** (`generate_plot`)
   - Creates plots showing RFD signal with annotations
   - Outputs: PNG plot files

4. **Combined Report** (`combine_summaries`)
   - Aggregates statistics across all samples
   - Outputs: TSV file with comparative data

## Output Files

For each input file `sample.bedgraph`, the workflow generates:

- `sample_critical_points.bed`: IS/TS locations
- `sample_zones.bed`: IZ/TZ regions  
- `sample_zones_detailed.tsv`: Detailed zone information
- `sample_summary.txt`: Statistical summary
- `sample_plot.png`: Visualization (optional)

Plus a combined summary: `combined_analysis_summary.tsv`

## Advanced Usage

### Run specific samples only
```bash
snakemake workflows/results/sample_name_critical_points.bed --cores 4
```

### Dry run to check workflow
```bash
snakemake --dry-run --cores 1
```

### Generate workflow report
```bash
snakemake --report report.html --cores 1
```

### Use different parameters
```bash
# Edit config.yaml or override on command line
snakemake --config min_peak_prominence=0.15 --cores 4
```

## Analysis Modes

### Standard Mode (default)
- `min_peak_distance: 10`
- `min_peak_prominence: 0.1`
- Good for typical RFD data

### Sensitive Mode  
- `min_peak_distance: 5`
- `min_peak_prominence: 0.05`
- Detects more peaks, good for weak signals

### Stringent Mode
- `min_peak_distance: 15` 
- `min_peak_prominence: 0.2`
- Fewer, higher-confidence peaks

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure bedgraph files are in `workflows/inputs/` (or `workflows/examples/` for samples)
2. **Permission errors**: Check write permissions for `workflows/results/`
3. **Memory issues**: Reduce parallelism with `--cores 1`
4. **Plot generation fails**: Set `generate_plots: false` in config.yaml

### Checking logs
```bash
# View individual job logs
ls workflows/results/logs/
cat workflows/results/logs/sample_name_analysis.log
```

### Re-running failed jobs
```bash
# Force re-run of specific rule
snakemake --forcerun generate_plot --cores 1
```

## Dependencies

- Python 3.7+
- snakemake
- pandas
- numpy  
- scipy
- matplotlib

## Performance Tips

- Use `--cores N` to parallelize across samples
- For large datasets, consider increasing memory allocation
- Use `--cluster` for HPC environments

## Example Commands

```bash
# Basic run (cores=1 by default)
python run_analysis.py
# or
.\run_analysis.ps1

# Use more cores for faster processing
python run_analysis.py --cores 4

# Sensitive analysis
python run_analysis.py --config analysis_mode=sensitive

# Without plots (faster)
python run_analysis.py --config generate_plots=false

# Specific chromosome only
python run_analysis.py --config chromosome=chr1

# Advanced: direct snakemake usage
snakemake --cores 4 --config analysis_mode=sensitive
```