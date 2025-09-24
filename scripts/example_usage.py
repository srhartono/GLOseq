#!/usr/bin/env python3
"""
Example usage of the RFD analyzer for OK-seq data analysis
"""

from rfd_analyzer import RFDAnalyzer
import matplotlib.pyplot as plt

def example_analysis():
    """Run a complete example analysis"""
    
    print("=== RFD Analysis Example ===\n")
    
    # Initialize analyzer with custom parameters
    analyzer = RFDAnalyzer(
        min_peak_distance=10,      # Minimum 10 bins between peaks
        min_peak_prominence=0.1    # Minimum prominence for peak detection
    )
    
    # Load sample data
    print("1. Loading sample data...")
    data = analyzer.load_bedgraph("sample_rfd_chr1.bedgraph")
    print(f"   Loaded {len(data)} genomic bins")
    print(f"   RFD range: {data['RFD'].min():.3f} to {data['RFD'].max():.3f}\n")
    
    # Find peaks and valleys
    print("2. Finding initiation and termination sites...")
    peaks_data = analyzer.find_peaks_and_valleys()
    
    for chrom, chrom_data in peaks_data.items():
        n_is = len(chrom_data['peaks']['positions'])
        n_ts = len(chrom_data['valleys']['positions'])
        print(f"   {chrom}: {n_is} initiation sites, {n_ts} termination sites")
    print()
    
    # Classify zones
    print("3. Classifying replication zones...")
    zones_data = analyzer.classify_zones()
    
    for chrom, chrom_data in zones_data.items():
        zones = chrom_data['zones']
        iz_count = sum(1 for z in zones if z['zone_type'] == 'IZ')
        tz_count = sum(1 for z in zones if z['zone_type'] == 'TZ')
        print(f"   {chrom}: {iz_count} initiation zones, {tz_count} termination zones")
    print()
    
    # Export results
    print("4. Exporting results...")
    analyzer.export_results("example_output")
    print()
    
    # Generate plot
    print("5. Generating visualization...")
    analyzer.plot_chromosome("chr1", "example_plot.png")
    
    print("=== Analysis Complete ===")
    print("\nOutput files:")
    print("- example_output_critical_points.bed: Initiation and termination sites")
    print("- example_output_zones.bed: Classified zones")
    print("- example_output_zones_detailed.tsv: Detailed zone information")
    print("- example_plot.png: Visualization")

def analyze_specific_region():
    """Example of analyzing a specific genomic region"""
    
    print("\n=== Region-Specific Analysis ===\n")
    
    analyzer = RFDAnalyzer()
    analyzer.load_bedgraph("sample_rfd_chr1.bedgraph")
    analyzer.find_peaks_and_valleys()
    analyzer.classify_zones()
    
    # Get data for chr1
    chrom_data = analyzer.peaks['chr1']
    zones = analyzer.zones['chr1']['zones']
    
    # Find zones in a specific region (e.g., 2-4 Mb)
    region_start, region_end = 2000000, 4000000
    
    print(f"Analyzing region chr1:{region_start:,}-{region_end:,}")
    
    # Filter zones in region
    region_zones = [z for z in zones if z['start'] >= region_start and z['end'] <= region_end]
    
    print(f"Found {len(region_zones)} zones in this region:")
    for zone in region_zones[:5]:  # Show first 5
        print(f"  {zone['zone_type']}: {zone['start']:,} - {zone['end']:,} ({zone['length']:,} bp)")
    
    if len(region_zones) > 5:
        print(f"  ... and {len(region_zones) - 5} more zones")

def compare_parameters():
    """Example comparing different peak detection parameters"""
    
    print("\n=== Parameter Comparison ===\n")
    
    # Test different prominence thresholds
    prominence_values = [0.05, 0.1, 0.2]
    
    for prom in prominence_values:
        analyzer = RFDAnalyzer(min_peak_prominence=prom)
        analyzer.load_bedgraph("sample_rfd_chr1.bedgraph")
        peaks_data = analyzer.find_peaks_and_valleys()
        
        chrom_data = peaks_data['chr1']
        n_is = len(chrom_data['peaks']['positions'])
        n_ts = len(chrom_data['valleys']['positions'])
        
        print(f"Prominence {prom}: {n_is} IS, {n_ts} TS")

if __name__ == "__main__":
    # Run complete example
    example_analysis()
    
    # Additional examples
    analyze_specific_region()
    compare_parameters()