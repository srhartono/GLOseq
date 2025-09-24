# Input Data Directory

Place your **real bedgraph files** here for analysis.

## File Format

Your bedgraph files should be tab-separated with 4 columns:
```
chr1    1000000    1010000    0.245
chr1    1010000    1020000    0.312
chr1    1020000    1030000    -0.156
...
```

Where:
- Column 1: Chromosome (chr1, chr2, etc.)
- Column 2: Start position 
- Column 3: End position
- Column 4: RFD value (-1 to +1)

## Usage

1. **Add your files**: Copy your `.bedgraph` files to this directory
2. **Run analysis**: Execute `python run_analysis.py` from the main directory
3. **Check results**: Find outputs in `../results/`

## For Testing

If you want to test with sample data first:
```bash
# From main directory
python setup_examples.py    # Copies example files here
python run_analysis.py      # Runs analysis
```

## File Naming

- Files should have `.bedgraph` extension
- Avoid spaces and special characters in filenames
- Example good names: `sample1.bedgraph`, `ctrl_rep1.bedgraph`

## Priority

The workflow will process files in this directory first. If this directory is empty, it will fall back to using files in `../examples/` for testing purposes.