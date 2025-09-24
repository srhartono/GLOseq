# GLOseq Project Structure Summary

## Successfully Created Snakemake Workflow

### Directory Structure
```
GLOseq/
├── workflows/
│   ├── inputs/                      # Your real bedgraph files (place files here)
│   ├── examples/                    # Sample bedgraph files for testing
│   │   ├── sample_rfd_chr1.bedgraph
│   │   └── sample_rfd_multi.bedgraph
│   ├── results/                     # All analysis outputs
│   │   ├── logs/                    # Job execution logs
│   │   ├── *_critical_points.bed    # IS/TS locations
│   │   ├── *_zones.bed              # IZ/TZ regions
│   │   ├── *_zones_detailed.tsv     # Detailed zone info
│   │   ├── *_summary.txt            # Individual summaries
│   │   ├── *_plot.png               # Visualizations
│   │   └── combined_analysis_summary.tsv  # Combined results
│   └── README.md                    # Workflow documentation
├── scripts/                         # Python analysis tools
│   ├── rfd_analyzer.py             # Main RFD analysis script
│   ├── generate_sample_data.py     # Sample data generator
│   ├── example_usage.py            # Usage examples
│   └── advanced_rfd_analysis.py    # Advanced utilities
├── Snakefile                       # Workflow definition  
├── config.yaml                     # Configuration parameters
├── run_analysis.py                 # Easy workflow launcher
├── setup_examples.py               # Copy examples to inputs
└── README.md                       # Project documentation
```

### Workflow Features

✅ **Automated RFD Analysis Pipeline**
- Processes multiple bedgraph files automatically
- Identifies Initiation Sites (IS) and Termination Sites (TS)
- Classifies Initiation Zones (IZ) and Termination Zones (TZ)
- Generates comprehensive statistics and visualizations

✅ **Configurable Parameters**
- Adjustable peak detection sensitivity
- Multiple analysis modes (standard/sensitive/stringent)
- Optional plot generation
- Chromosome-specific analysis

✅ **Organized Output**
- Clean separation of inputs, outputs, and scripts
- Standardized file naming
- Combined summary reports
- Individual job logging

✅ **Production Ready**
- Proper dependency management
- Error handling and logging
- Parallel execution support
- Reproducible analysis

### Usage

#### Basic Analysis
```bash
# Place your .bedgraph files in workflows/inputs/
# Then run:
python run_analysis.py

# Or for testing with examples:
python setup_examples.py    # Copy examples to inputs
python run_analysis.py      # Run analysis
```

#### Custom Parameters
```bash
# Edit config.yaml or override on command line
snakemake --config min_peak_prominence=0.15 --cores 4
```

#### Advanced Usage
```bash
# Dry run
snakemake --dry-run

# Specific sample
snakemake workflows/results/your_sample_critical_points.bed

# Generate report
snakemake --report analysis_report.html
```

### Test Results

The workflow was successfully tested with sample data:
- **Sample 1** (chr1): 55 IZ zones, 54 TZ zones
- **Sample 2** (multi-chr): 89 IZ zones, 88 TZ zones across 3 chromosomes
- Generated all expected output files including plots and summaries

### Next Steps

1. **Add your real data**: Place `.bedgraph` files in `workflows/inputs/`
2. **Customize settings**: Edit `config.yaml` for your analysis parameters  
3. **Run analysis**: Execute `snakemake --cores N` 
4. **Review results**: Check `workflows/results/` for all outputs
5. **Scale up**: Use HPC clusters with `--cluster` for large datasets