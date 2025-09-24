# Quick Start Guide - RFD Analysis

## Simple 3-Step Process

### 1. Add Your Data
Place your `.bedgraph` files in the `workflows/inputs/` folder
(Sample files are available in `workflows/examples/` for testing)

### 2. Run Analysis
Choose the easiest option for your system:

#### **Recommended: Python Wrapper (Works Everywhere)**
```bash
python run_analysis.py
```

#### **Windows PowerShell**
```powershell
.\run_analysis.ps1
```

#### **Windows Command Prompt**
```cmd
run_analysis.bat
```

### 3. Check Results
Find all outputs in `workflows/results/` folder:
- `*_critical_points.bed` - Initiation/Termination sites
- `*_zones.bed` - IZ/TZ regions
- `*_zones_detailed.tsv` - Detailed statistics
- `*_summary.txt` - Analysis summary
- `*_plot.png` - Visualizations
- `combined_analysis_summary.tsv` - Combined results

## Common Options

### Use More CPU Cores (Faster)
```bash
python run_analysis.py --cores 4
```

### Change Analysis Sensitivity
```bash
python run_analysis.py --config analysis_mode=sensitive    # More peaks
python run_analysis.py --config analysis_mode=stringent    # Fewer peaks
```

### Skip Plots (Faster)
```bash
python run_analysis.py --config generate_plots=false
```

### Analyze Specific Chromosome
```bash
python run_analysis.py --config chromosome=chr1
```

## Troubleshooting

### No Input Files Found
- Make sure `.bedgraph` files are in `workflows/inputs/`
- For testing, copy from `workflows/examples/` to `workflows/inputs/`
- Check file permissions

### Analysis Fails
- Check logs in `workflows/results/logs/`
- Try with `--cores 1` if memory issues
- Verify bedgraph file format (4 columns: chr, start, end, RFD)

### Need Help
- Check `workflows/README.md` for detailed documentation
- View example files in `workflows/examples/`

## File Format
Your bedgraph files should look like:
```
chr1    1000000    1010000    0.245
chr1    1010000    1020000    0.312
chr1    1020000    1030000    -0.156
```

That's it! The analysis runs automatically and generates all outputs.