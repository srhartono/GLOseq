# RFD Analysis Tool for OK-seq Data

This tool analyzes Replication Fork Directionality (RFD) data from OK-seq experiments to identify replication features and zones.

## Features

- **Initiation Sites (IS)**: Local maxima in RFD values indicating where replication initiates
- **Termination Sites (TS)**: Local minima in RFD values indicating where replication terminates  
- **Initiation Zones (IZ)**: Regions from termination sites (TS) to initiation sites (IS)
- **Termination Zones (TZ)**: Regions from initiation sites (IS) to termination sites (TS)

## Input Format

The tool expects a bedgraph file with the following tab-separated columns:
```
chr    start    end    RFD_value
chr1   1000000  1010000  0.245
chr1   1010000  1020000  0.312
...
```

Where:
- `chr`: Chromosome name (e.g., chr1, chr2, chrX)
- `start`: Start position of 10kb bin
- `end`: End position of 10kb bin  
- `RFD_value`: RFD value ranging from -1 to +1

## Usage

### Snakemake Workflow (Recommended)

```bash
# Place your .bedgraph files in workflows/inputs/ then run:
python run_analysis.py              # Easy way (cores=1 default)
python run_analysis.py --cores 4    # Use more cores
```

### Direct Script Usage

```bash
python scripts/rfd_analyzer.py input_file.bedgraph
```

### Advanced Usage

```bash
python scripts/rfd_analyzer.py input_file.bedgraph \
    -o output_prefix \
    --min-distance 10 \
    --min-prominence 0.1 \
    --chromosome chr1 \
    --plot
```

### Parameters

- `-o, --output`: Output prefix for result files (default: 'rfd_analysis')
- `-d, --min-distance`: Minimum distance between peaks in bins (default: 10)
- `-p, --min-prominence`: Minimum peak prominence for detection (default: 0.1)
- `-c, --chromosome`: Analyze specific chromosome only
- `--plot`: Generate plots for visualization

## Output Files

1. **`{prefix}_critical_points.bed`**: BED format file with all initiation and termination sites
   ```
   chr1  1040000  1050000  IS  0.368
   chr1  1150000  1160000  TS  0.334
   ```

2. **`{prefix}_zones.bed`**: BED format file with classified zones
   ```
   chr1  1040000  1150000  TZ
   chr1  1150000  1280000  IZ
   ```

3. **`{prefix}_zones_detailed.tsv`**: Detailed zone information with transitions
   ```
   chr    start    end      zone_type  start_feature  end_feature  length
   chr1   1040000  1150000  TZ         IS             TS           110000
   chr1   1150000  1280000  IZ         TS             IS           130000
   ```

4. **`{prefix}_{chr}_plot.png`**: Visualization plots (if --plot specified)

## Zone Classification Logic

- **IZ (Initiation Zone)**: Region from TS → IS (replication initiation activity)
- **TZ (Termination Zone)**: Region from IS → TS (replication termination activity)
- **Transition zones**: Handle consecutive same-type features

## Example Analysis

```python
from rfd_analyzer import RFDAnalyzer

# Initialize analyzer
analyzer = RFDAnalyzer(min_peak_distance=10, min_peak_prominence=0.1)

# Load data
analyzer.load_bedgraph("your_data.bedgraph")

# Find peaks and valleys
analyzer.find_peaks_and_valleys()

# Classify zones
analyzer.classify_zones()

# Export results
analyzer.export_results("my_analysis")

# Plot specific chromosome
analyzer.plot_chromosome("chr1", "chr1_plot.png")
```

## Requirements

- Python 3.7+
- pandas
- numpy
- scipy
- matplotlib

## Installation

```bash
pip install pandas numpy scipy matplotlib
```

## Sample Data Generation

Use the included sample data generator for testing:

```bash
python generate_sample_data.py
```

This creates:
- `sample_rfd_chr1.bedgraph`: Single chromosome sample
- `sample_rfd_multi.bedgraph`: Multi-chromosome sample

## Notes

- RFD values should range from -1 to +1
- The tool assumes 10kb genomic bins
- Peak detection parameters may need adjustment for different datasets
- For noisy data, consider increasing `min_prominence` parameter
- For high-resolution data, consider decreasing `min_distance` parameter