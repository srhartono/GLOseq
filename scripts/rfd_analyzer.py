#!/usr/bin/env python3
"""
RFD (Replication Fork Directionality) Analysis Tool for OK-seq Data

This script analyzes RFD data from bedgraph files to identify:
- Initiation Sites (IS): local maxima in RFD values
- Termination Sites (TS): local minima in RFD values  
- Initiation Zones (IZ): regions from minima to maxima
- Termination Zones (TZ): regions from maxima to minima

Input: Bedgraph file with columns: chr, start, end, RFD_value
Output: Annotated regions with zone classifications
"""

import pandas as pd
import numpy as np
from scipy.signal import find_peaks, argrelextrema
import matplotlib.pyplot as plt
import argparse
import sys
from typing import List, Tuple, Dict

class RFDAnalyzer:
    def __init__(self, min_peak_distance: int = 10, min_peak_prominence: float = 0.1):
        """
        Initialize the RFD analyzer.
        
        Args:
            min_peak_distance: Minimum distance between peaks (in bins)
            min_peak_prominence: Minimum prominence for peak detection
        """
        self.min_peak_distance = min_peak_distance
        self.min_peak_prominence = min_peak_prominence
        self.data = None
        self.peaks = {}
        self.zones = {}
    
    def load_bedgraph(self, filepath: str) -> pd.DataFrame:
        """
        Load RFD data from bedgraph file (supports both .bedgraph and .bedgraph.gz).
        
        Args:
            filepath: Path to bedgraph file (.bedgraph or .bedgraph.gz)
            
        Returns:
            DataFrame with RFD data
        """
        try:
            import gzip
            
            # Check if file is gzipped based on extension or by trying to read it
            if filepath.endswith('.gz'):
                # Read gzipped file
                with gzip.open(filepath, 'rt') as f:
                    self.data = pd.read_csv(f, sep='\t', 
                                          names=['chr', 'start', 'end', 'RFD'],
                                          comment='#')
            else:
                # Read regular file
                self.data = pd.read_csv(filepath, sep='\t', 
                                      names=['chr', 'start', 'end', 'RFD'],
                                      comment='#')
            
            # Sort by chromosome and position
            self.data = self.data.sort_values(['chr', 'start']).reset_index(drop=True)
            
            print(f"Loaded {len(self.data)} intervals from {filepath}")
            print(f"RFD values range: {self.data['RFD'].min():.3f} to {self.data['RFD'].max():.3f}")
            
            return self.data
            
        except Exception as e:
            print(f"Error loading file {filepath}: {e}")
            sys.exit(1)
    
    def find_peaks_and_valleys(self, chromosome: str = None) -> Dict:
        """
        Find local maxima (IS) and minima (TS) in RFD data.
        
        Args:
            chromosome: Specific chromosome to analyze (if None, analyze all)
            
        Returns:
            Dictionary with peak and valley information
        """
        if self.data is None:
            raise ValueError("No data loaded. Call load_bedgraph() first.")
        
        results = {}
        
        # Get unique chromosomes
        chromosomes = [chromosome] if chromosome else self.data['chr'].unique()
        
        for chrom in chromosomes:
            chrom_data = self.data[self.data['chr'] == chrom].copy()
            
            if len(chrom_data) < 3:
                print(f"Warning: Insufficient data for chromosome {chrom}")
                continue
            
            rfd_values = chrom_data['RFD'].values
            positions = chrom_data['start'].values
            
            # Find peaks (maxima) - Initiation Sites
            peaks, peak_properties = find_peaks(
                rfd_values,
                distance=self.min_peak_distance,
                prominence=self.min_peak_prominence
            )
            
            # Find valleys (minima) - Termination Sites
            valleys, valley_properties = find_peaks(
                -rfd_values,  # Invert for valley detection
                distance=self.min_peak_distance,
                prominence=self.min_peak_prominence
            )
            
            results[chrom] = {
                'data': chrom_data,
                'peaks': {
                    'indices': peaks,
                    'positions': positions[peaks] if len(peaks) > 0 else np.array([]),
                    'values': rfd_values[peaks] if len(peaks) > 0 else np.array([]),
                    'properties': peak_properties
                },
                'valleys': {
                    'indices': valleys,
                    'positions': positions[valleys] if len(valleys) > 0 else np.array([]),
                    'values': rfd_values[valleys] if len(valleys) > 0 else np.array([]),
                    'properties': valley_properties
                }
            }
            
            print(f"Chromosome {chrom}: Found {len(peaks)} initiation sites, {len(valleys)} termination sites")
        
        self.peaks = results
        return results
    
    def classify_zones(self) -> Dict:
        """
        Classify regions into Initiation Zones (IZ) and Termination Zones (TZ).
        
        IZ: from minima (TS) to maxima (IS)
        TZ: from maxima (IS) to minima (TS)
        
        Returns:
            Dictionary with zone classifications
        """
        if not self.peaks:
            raise ValueError("No peaks found. Call find_peaks_and_valleys() first.")
        
        zones = {}
        
        for chrom, chrom_data in self.peaks.items():
            peaks_pos = chrom_data['peaks']['positions']
            valleys_pos = chrom_data['valleys']['positions']
            
            # Combine and sort all critical points
            all_points = []
            
            # Add peaks (IS)
            for pos in peaks_pos:
                all_points.append((pos, 'IS'))
            
            # Add valleys (TS)
            for pos in valleys_pos:
                all_points.append((pos, 'TS'))
            
            # Sort by position
            all_points.sort(key=lambda x: x[0])
            
            # Classify zones between consecutive points
            zone_list = []
            
            for i in range(len(all_points) - 1):
                start_pos, start_type = all_points[i]
                end_pos, end_type = all_points[i + 1]
                
                # Determine zone type based on transition
                if start_type == 'TS' and end_type == 'IS':
                    zone_type = 'IZ'  # Initiation Zone: TS -> IS
                elif start_type == 'IS' and end_type == 'TS':
                    zone_type = 'TZ'  # Termination Zone: IS -> TS
                else:
                    # Handle consecutive same types (rare, but possible)
                    zone_type = f"Transition_{start_type}_to_{end_type}"
                
                zone_list.append({
                    'chr': chrom,
                    'start': start_pos,
                    'end': end_pos,
                    'zone_type': zone_type,
                    'start_feature': start_type,
                    'end_feature': end_type,
                    'length': end_pos - start_pos
                })
            
            zones[chrom] = {
                'critical_points': all_points,
                'zones': zone_list
            }
            
            print(f"Chromosome {chrom}: Classified {len(zone_list)} zones")
        
        self.zones = zones
        return zones
    
    def export_results(self, output_prefix: str):
        """
        Export results to various output files.
        
        Args:
            output_prefix: Prefix for output files
        """
        if not self.peaks or not self.zones:
            raise ValueError("No results to export. Run analysis first.")
        
        # Export critical points (IS and TS)
        critical_points = []
        for chrom, data in self.peaks.items():
            # Add IS (peaks) - Red color
            for i, pos in enumerate(data['peaks']['positions']):
                rfd_val = data['peaks']['values'][i]
                score = min(1000, max(0, int((abs(rfd_val) * 1000))))  # Scale RFD to 0-1000
                critical_points.append({
                    'chr': chrom,
                    'start': pos,
                    'end': pos + 10000,  # 10kb bins
                    'name': f'IS_{abs(rfd_val):.3f}',
                    'score': score,
                    'strand': '.',
                    'thickStart': pos,
                    'thickEnd': pos + 10000,
                    'itemRgb': '255,0,0'  # Red for IS
                })
            
            # Add TS (valleys) - Blue color
            for i, pos in enumerate(data['valleys']['positions']):
                rfd_val = data['valleys']['values'][i]
                score = min(1000, max(0, int((abs(rfd_val) * 1000))))  # Scale RFD to 0-1000
                critical_points.append({
                    'chr': chrom,
                    'start': pos,
                    'end': pos + 10000,  # 10kb bins
                    'name': f'TS_{abs(rfd_val):.3f}',
                    'score': score,
                    'strand': '.',
                    'thickStart': pos,
                    'thickEnd': pos + 10000,
                    'itemRgb': '0,0,255'  # Blue for TS
                })
        
        # Save critical points in BED9 format with colors
        critical_df = pd.DataFrame(critical_points)
        critical_df = critical_df.sort_values(['chr', 'start'])
        # Reorder columns for BED9 format
        bed9_columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']
        critical_df = critical_df[bed9_columns]
        critical_df.to_csv(f"{output_prefix}_critical_points.bed", 
                          sep='\t', index=False, header=False)
        
        # Export zones
        all_zones = []
        for chrom, data in self.zones.items():
            for zone in data['zones']:
                all_zones.append(zone)
        
        zones_df = pd.DataFrame(all_zones)
        zones_df = zones_df.sort_values(['chr', 'start'])
        
        # Save zones as BED9 file with colors
        zones_bed = zones_df.copy()
        # Add BED9 columns with colors
        zones_bed['name'] = zones_bed['zone_type'] + '_' + zones_bed['length'].astype(str)
        zones_bed['score'] = zones_bed['length'].apply(lambda x: min(1000, max(0, int(x/1000))))  # Scale length to score
        zones_bed['strand'] = '.'
        zones_bed['thickStart'] = zones_bed['start']
        zones_bed['thickEnd'] = zones_bed['end']
        # Set colors: IZ = green (0,128,0), TZ = orange (255,165,0)
        zones_bed['itemRgb'] = zones_bed['zone_type'].apply(
            lambda x: '255,100,0' if x == 'IZ' else
            '0,100,255' if x == "TZ" else
            '255,200,0' if x == "Transition_IS_to_IS" else
            '0,200,255' if x == "Transition_TS_to_TS" else '0,155,0'
        )
        
        # Select and reorder columns for BED9 format
        bed9_zones = zones_bed[['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']]
        bed9_zones.to_csv(f"{output_prefix}_zones.bed", 
                         sep='\t', index=False, header=False)
        
        # Save detailed zones with all info
        zones_df.to_csv(f"{output_prefix}_zones_detailed.tsv", 
                       sep='\t', index=False)
        
        print(f"Exported {len(critical_points)} critical points to {output_prefix}_critical_points.bed")
        print(f"Exported {len(all_zones)} zones to {output_prefix}_zones.bed")
        print(f"Exported detailed zone info to {output_prefix}_zones_detailed.tsv")
    
    def plot_chromosome(self, chromosome: str, output_file: str = None):
        """
        Plot RFD data for a specific chromosome with annotations.
        
        Args:
            chromosome: Chromosome to plot
            output_file: Optional output file for saving plot
        """
        if chromosome not in self.peaks:
            raise ValueError(f"No data found for chromosome {chromosome}")
        
        chrom_data = self.peaks[chromosome]
        data = chrom_data['data']
        
        plt.figure(figsize=(15, 6))
        
        # Plot RFD values
        plt.plot(data['start'], data['RFD'], 'b-', linewidth=1, alpha=0.7, label='RFD')
        
        # Plot peaks (IS)
        if len(chrom_data['peaks']['positions']) > 0:
            plt.scatter(chrom_data['peaks']['positions'], 
                       chrom_data['peaks']['values'],
                       color='red', s=50, marker='^', 
                       label=f"IS (n={len(chrom_data['peaks']['positions'])})", 
                       zorder=5)
        
        # Plot valleys (TS)
        if len(chrom_data['valleys']['positions']) > 0:
            plt.scatter(chrom_data['valleys']['positions'], 
                       chrom_data['valleys']['values'],
                       color='blue', s=50, marker='v', 
                       label=f"TS (n={len(chrom_data['valleys']['positions'])})", 
                       zorder=5)
        
        # Add zone coloring if available
        if chromosome in self.zones:
            for zone in self.zones[chromosome]['zones']:
                if zone['zone_type'] == 'IZ':
                    plt.axvspan(zone['start'], zone['end'], alpha=0.2, color='green', label='IZ')
                elif zone['zone_type'] == 'TZ':
                    plt.axvspan(zone['start'], zone['end'], alpha=0.2, color='orange', label='TZ')
        
        plt.xlabel('Geznomic Pzosition')
        plt.ylabel('RFD Value')
        plt.title(f'RFD Analysis - Chromosome {chromosome}')
        # plt.legend()
        plt.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        
        # plt.show()

def main():
    parser = argparse.ArgumentParser(description='Analyze RFD data from OK-seq')
    parser.add_argument('input_file', help='Input bedgraph file')
    parser.add_argument('-o', '--output', default='rfd_analysis', 
                       help='Output prefix for result files')
    parser.add_argument('-d', '--min-distance', type=int, default=1000000,
                       help='Minimum distance between peaks (bins)')
    parser.add_argument('-p', '--min-prominence', type=float, default=1,
                       help='Minimum peak prominence')
    parser.add_argument('-c', '--chromosome', 
                       help='Analyze specific chromosome only')
    parser.add_argument('--plot', action='store_true',
                       help='Generate plots for each chromosome')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = RFDAnalyzer(
        min_peak_distance=args.min_distance,
        min_peak_prominence=args.min_prominence
    )
    
    # Load data
    analyzer.load_bedgraph(args.input_file)
    
    # Find peaks and valleys
    analyzer.find_peaks_and_valleys(args.chromosome)
    
    # Classify zones
    analyzer.classify_zones()
    
    # Export results
    analyzer.export_results(args.output)
    
    # Generate plots if requested
    if args.plot:
        chromosomes = [args.chromosome] if args.chromosome else analyzer.peaks.keys()
        for chrom in chromosomes:
            analyzer.plot_chromosome(chrom, f"{args.output}_{chrom}_plot.png")

if __name__ == "__main__":
    main()