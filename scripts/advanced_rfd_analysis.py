#!/usr/bin/env python3
"""
Advanced RFD analysis utilities and helper functions
"""

import pandas as pd
import numpy as np
from rfd_analyzer import RFDAnalyzer
import matplotlib.pyplot as plt

def validate_bedgraph_format(filepath: str) -> bool:
    """
    Validate that the input file is in proper bedgraph format (supports .bedgraph and .bedgraph.gz)
    
    Args:
        filepath: Path to bedgraph file (.bedgraph or .bedgraph.gz)
        
    Returns:
        True if format is valid, False otherwise
    """
    try:
        import gzip
        
        # Read first few lines, handling both regular and gzipped files
        if filepath.endswith('.gz'):
            with gzip.open(filepath, 'rt') as f:
                sample = pd.read_csv(f, sep='\t', nrows=10, header=None)
        else:
            sample = pd.read_csv(filepath, sep='\t', nrows=10, header=None)
        
        # Check number of columns
        if sample.shape[1] != 4:
            print(f"Error: Expected 4 columns, found {sample.shape[1]}")
            return False
        
        # Check if RFD values are in reasonable range
        rfd_values = sample.iloc[:, 3]
        if rfd_values.min() < -2 or rfd_values.max() > 2:
            print(f"Warning: RFD values may be out of expected range (-1 to +1)")
            print(f"Found range: {rfd_values.min():.3f} to {rfd_values.max():.3f}")
        
        # Check if positions are sorted
        positions = sample.iloc[:, 1]
        if not positions.is_monotonic_increasing:
            print("Warning: Positions may not be sorted")
        
        print("Bedgraph format validation passed")
        return True
        
    except Exception as e:
        print(f"Error validating file: {e}")
        return False

def filter_chromosomes(data: pd.DataFrame, include_chroms: list = None, 
                      exclude_chroms: list = None) -> pd.DataFrame:
    """
    Filter data to include/exclude specific chromosomes
    
    Args:
        data: Input DataFrame
        include_chroms: List of chromosomes to include
        exclude_chroms: List of chromosomes to exclude
        
    Returns:
        Filtered DataFrame
    """
    filtered_data = data.copy()
    
    if include_chroms:
        filtered_data = filtered_data[filtered_data['chr'].isin(include_chroms)]
        print(f"Filtered to include chromosomes: {include_chroms}")
    
    if exclude_chroms:
        filtered_data = filtered_data[~filtered_data['chr'].isin(exclude_chroms)]
        print(f"Excluded chromosomes: {exclude_chroms}")
    
    return filtered_data

def smooth_rfd_data(data: pd.DataFrame, window_size: int = 3) -> pd.DataFrame:
    """
    Apply smoothing to RFD data to reduce noise
    
    Args:
        data: Input DataFrame with RFD column
        window_size: Size of smoothing window
        
    Returns:
        DataFrame with smoothed RFD values
    """
    smoothed_data = data.copy()
    
    for chrom in data['chr'].unique():
        chrom_mask = data['chr'] == chrom
        rfd_values = data.loc[chrom_mask, 'RFD']
        
        # Apply rolling mean
        smoothed_values = rfd_values.rolling(
            window=window_size, 
            center=True, 
            min_periods=1
        ).mean()
        
        smoothed_data.loc[chrom_mask, 'RFD'] = smoothed_values
    
    print(f"Applied smoothing with window size {window_size}")
    return smoothed_data

def analyze_zone_statistics(analyzer: RFDAnalyzer) -> dict:
    """
    Calculate comprehensive statistics for zones
    
    Args:
        analyzer: RFDAnalyzer object with completed analysis
        
    Returns:
        Dictionary with statistics
    """
    if not analyzer.zones:
        raise ValueError("No zones found. Run analysis first.")
    
    stats = {}
    
    for chrom, chrom_data in analyzer.zones.items():
        zones = chrom_data['zones']
        
        # Separate by zone type
        iz_zones = [z for z in zones if z['zone_type'] == 'IZ']
        tz_zones = [z for z in zones if z['zone_type'] == 'TZ']
        
        # Calculate statistics
        iz_lengths = [z['length'] for z in iz_zones]
        tz_lengths = [z['length'] for z in tz_zones]
        
        chrom_stats = {
            'total_zones': len(zones),
            'iz_count': len(iz_zones),
            'tz_count': len(tz_zones),
            'iz_mean_length': np.mean(iz_lengths) if iz_lengths else 0,
            'iz_median_length': np.median(iz_lengths) if iz_lengths else 0,
            'tz_mean_length': np.mean(tz_lengths) if tz_lengths else 0,
            'tz_median_length': np.median(tz_lengths) if tz_lengths else 0,
            'total_iz_coverage': sum(iz_lengths),
            'total_tz_coverage': sum(tz_lengths)
        }
        
        stats[chrom] = chrom_stats
    
    return stats

def create_summary_report(analyzer: RFDAnalyzer, output_file: str = "analysis_summary.txt"):
    """
    Create a comprehensive summary report
    
    Args:
        analyzer: RFDAnalyzer object with completed analysis
        output_file: Output file for summary report
    """
    stats = analyze_zone_statistics(analyzer)
    
    with open(output_file, 'w') as f:
        f.write("RFD Analysis Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        # Overall statistics
        total_is = sum(len(chrom_data['peaks']['positions']) 
                      for chrom_data in analyzer.peaks.values())
        total_ts = sum(len(chrom_data['valleys']['positions']) 
                      for chrom_data in analyzer.peaks.values())
        
        f.write(f"Overall Statistics:\n")
        f.write(f"  Total Initiation Sites (IS): {total_is}\n")
        f.write(f"  Total Termination Sites (TS): {total_ts}\n")
        f.write(f"  Chromosomes analyzed: {len(analyzer.peaks)}\n\n")
        
        # Per-chromosome statistics
        for chrom, chrom_stats in stats.items():
            f.write(f"Chromosome {chrom}:\n")
            f.write(f"  Zones: {chrom_stats['total_zones']} total "
                   f"({chrom_stats['iz_count']} IZ, {chrom_stats['tz_count']} TZ)\n")
            f.write(f"  IZ lengths: mean={chrom_stats['iz_mean_length']:,.0f} bp, "
                   f"median={chrom_stats['iz_median_length']:,.0f} bp\n")
            f.write(f"  TZ lengths: mean={chrom_stats['tz_mean_length']:,.0f} bp, "
                   f"median={chrom_stats['tz_median_length']:,.0f} bp\n")
            f.write(f"  Coverage: {chrom_stats['total_iz_coverage']:,.0f} bp IZ, "
                   f"{chrom_stats['total_tz_coverage']:,.0f} bp TZ\n\n")
        
        # Analysis parameters
        f.write(f"Analysis Parameters:\n")
        f.write(f"  Minimum peak distance: {analyzer.min_peak_distance} bins\n")
        f.write(f"  Minimum peak prominence: {analyzer.min_peak_prominence}\n")
    
    print(f"Summary report saved to: {output_file}")

def batch_analyze_multiple_files(file_list: list, output_dir: str = "batch_results"):
    """
    Analyze multiple bedgraph files in batch
    
    Args:
        file_list: List of bedgraph file paths
        output_dir: Output directory for results
    """
    import os
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    results_summary = []
    
    for i, filepath in enumerate(file_list, 1):
        print(f"\nProcessing file {i}/{len(file_list)}: {filepath}")
        
        try:
            # Initialize analyzer
            analyzer = RFDAnalyzer()
            
            # Load and analyze
            analyzer.load_bedgraph(filepath)
            analyzer.find_peaks_and_valleys()
            analyzer.classify_zones()
            
            # Create output prefix
            filename = os.path.basename(filepath)
            # Remove both .bedgraph.gz and .bedgraph extensions
            if filename.endswith('.bedgraph.gz'):
                filename = filename.replace('.bedgraph.gz', '')
            elif filename.endswith('.bedgraph'):
                filename = filename.replace('.bedgraph', '')
            output_prefix = os.path.join(output_dir, filename)
            
            # Export results
            analyzer.export_results(output_prefix)
            create_summary_report(analyzer, f"{output_prefix}_summary.txt")
            
            # Collect summary statistics
            stats = analyze_zone_statistics(analyzer)
            total_is = sum(len(chrom_data['peaks']['positions']) 
                          for chrom_data in analyzer.peaks.values())
            total_ts = sum(len(chrom_data['valleys']['positions']) 
                          for chrom_data in analyzer.peaks.values())
            
            results_summary.append({
                'file': filepath,
                'total_is': total_is,
                'total_ts': total_ts,
                'chromosomes': len(analyzer.peaks),
                'status': 'Success'
            })
            
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            results_summary.append({
                'file': filepath,
                'total_is': 0,
                'total_ts': 0,
                'chromosomes': 0,
                'status': f'Error: {str(e)}'
            })
    
    # Save batch summary
    summary_df = pd.DataFrame(results_summary)
    summary_df.to_csv(os.path.join(output_dir, "batch_summary.csv"), index=False)
    print(f"\nBatch analysis complete. Results in {output_dir}/")

if __name__ == "__main__":
    # Example of advanced analysis workflow
    
    print("Advanced RFD Analysis Example\n")
    
    # 1. Validate file format
    if validate_bedgraph_format("sample_rfd_chr1.bedgraph.gz"):
        
        # 2. Load and preprocess data
        analyzer = RFDAnalyzer(min_peak_prominence=0.15)
        
        # Load data
        data = analyzer.load_bedgraph("sample_rfd_chr1.bedgraph.gz")
        
        # Optional: Apply smoothing for noisy data
        # smoothed_data = smooth_rfd_data(data, window_size=3)
        
        # 3. Run analysis
        analyzer.find_peaks_and_valleys()
        analyzer.classify_zones()
        
        # 4. Generate comprehensive outputs
        analyzer.export_results("advanced_analysis")
        create_summary_report(analyzer, "advanced_summary.txt")
        
        # 5. Display statistics
        stats = analyze_zone_statistics(analyzer)
        
        print("\n=== Zone Statistics ===")
        for chrom, chrom_stats in stats.items():
            print(f"{chrom}: {chrom_stats['iz_count']} IZ, {chrom_stats['tz_count']} TZ")
            print(f"  IZ avg length: {chrom_stats['iz_mean_length']:,.0f} bp")
            print(f"  TZ avg length: {chrom_stats['tz_mean_length']:,.0f} bp")
    
    print("\nAdvanced analysis complete!")