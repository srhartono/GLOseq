#!/usr/bin/env python3
"""
Enhanced TRC Peak Distribution Analysis

This script analyzes peaks and distributions of transcription-replication conflicts
and provides detailed statistics for co-directional, neutral, and head-on conflicts.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os
from scipy import stats
from scipy.signal import find_peaks
import logging
from typing import Dict, List, Tuple

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TRCPeakAnalyzer:
    def __init__(self):
        """Initialize the TRC Peak Distribution analyzer."""
        self.trc_results = None
        self.critical_points = None
        self.zones = None
        self.rfd_data = None
        
    def load_trc_results(self, trc_file: str) -> pd.DataFrame:
        """Load TRC analysis results."""
        logger.info(f"Loading TRC results from: {trc_file}")
        
        try:
            df = pd.read_csv(trc_file, sep='\t')
            logger.info(f"Loaded TRC results for {len(df)} genes")
            
            # Filter out genes without conflict data
            df = df[df['conflict_type'].isin(['head-on', 'co-directional', 'neutral'])]
            logger.info(f"Analyzing {len(df)} genes with valid conflict classifications")
            
            self.trc_results = df
            return df
            
        except Exception as e:
            logger.error(f"Error loading TRC results: {e}")
            return None
    
    def load_critical_points(self, critical_points_file: str) -> pd.DataFrame:
        """Load IS/TS critical points."""
        logger.info(f"Loading critical points from: {critical_points_file}")
        
        try:
            df = pd.read_csv(critical_points_file, sep='\t', header=None)
            df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end', 'color']
            
            # Extract site type (IS/TS) and RFD score from name
            df['site_type'] = df['name'].str.extract(r'(IS|TS)_')
            df['rfd_score'] = df['name'].str.extract(r'_([\d.-]+)').astype(float)
            df['midpoint'] = (df['start'] + df['end']) // 2
            
            logger.info(f"Loaded {len(df)} critical points")
            logger.info(f"Site type distribution: {df['site_type'].value_counts().to_dict()}")
            
            self.critical_points = df
            return df
            
        except Exception as e:
            logger.error(f"Error loading critical points: {e}")
            return None
    
    def analyze_conflict_distributions(self) -> Dict:
        """Analyze distributions of RFD values for each conflict type."""
        if self.trc_results is None:
            logger.error("No TRC results loaded")
            return None
        
        logger.info("Analyzing conflict type distributions...")
        
        distributions = {}
        
        for conflict_type in ['head-on', 'co-directional', 'neutral']:
            subset = self.trc_results[self.trc_results['conflict_type'] == conflict_type]
            
            if len(subset) == 0:
                logger.warning(f"No genes found for conflict type: {conflict_type}")
                continue
            
            rfd_values = subset['rfd_value'].dropna()
            
            # Basic statistics
            stats_dict = {
                'count': len(rfd_values),
                'mean': rfd_values.mean(),
                'median': rfd_values.median(),
                'std': rfd_values.std(),
                'min': rfd_values.min(),
                'max': rfd_values.max(),
                'q25': rfd_values.quantile(0.25),
                'q75': rfd_values.quantile(0.75),
                'iqr': rfd_values.quantile(0.75) - rfd_values.quantile(0.25)
            }
            
            # Distribution shape analysis
            skewness = stats.skew(rfd_values)
            kurtosis = stats.kurtosis(rfd_values)
            
            # Find peaks in RFD value distribution
            hist, bin_edges = np.histogram(rfd_values, bins=20)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Find peaks in the histogram
            peaks, properties = find_peaks(hist, height=1, distance=2)
            peak_positions = bin_centers[peaks] if len(peaks) > 0 else []
            peak_heights = hist[peaks] if len(peaks) > 0 else []
            
            distributions[conflict_type] = {
                'statistics': stats_dict,
                'shape': {
                    'skewness': skewness,
                    'kurtosis': kurtosis
                },
                'peaks': {
                    'positions': peak_positions.tolist(),
                    'heights': peak_heights.tolist(),
                    'count': len(peaks)
                },
                'raw_values': rfd_values.tolist()
            }
        
        return distributions
    
    def find_conflict_hotspots(self) -> Dict:
        """Find genomic regions with high concentrations of specific conflict types."""
        if self.trc_results is None:
            logger.error("No TRC results loaded")
            return None
        
        logger.info("Finding conflict hotspots...")
        
        hotspots = {}
        
        for conflict_type in ['head-on', 'co-directional']:
            subset = self.trc_results[self.trc_results['conflict_type'] == conflict_type]
            
            if len(subset) == 0:
                continue
            
            # Group by chromosome
            chr_hotspots = []
            
            for chromosome in subset['chromosome'].unique():
                chr_data = subset[subset['chromosome'] == chromosome].copy()
                
                if len(chr_data) < 3:  # Need at least 3 genes for hotspot analysis
                    continue
                
                # Sort by position
                chr_data = chr_data.sort_values('midpoint')
                
                # Find clusters of genes within 500kb
                window_size = 500000  # 500kb window
                hotspot_regions = []
                
                for i, gene in chr_data.iterrows():
                    # Count genes within window
                    window_start = gene['midpoint'] - window_size // 2
                    window_end = gene['midpoint'] + window_size // 2
                    
                    genes_in_window = chr_data[
                        (chr_data['midpoint'] >= window_start) & 
                        (chr_data['midpoint'] <= window_end)
                    ]
                    
                    if len(genes_in_window) >= 3:  # Hotspot threshold
                        hotspot_regions.append({
                            'chromosome': chromosome,
                            'center': gene['midpoint'],
                            'start': window_start,
                            'end': window_end,
                            'gene_count': len(genes_in_window),
                            'genes': genes_in_window['gene_name'].tolist(),
                            'mean_rfd': genes_in_window['rfd_value'].mean(),
                            'density': len(genes_in_window) / (window_size / 1000000)  # genes per Mb
                        })
                
                # Remove overlapping hotspots (keep the one with highest density)
                if hotspot_regions:
                    hotspot_df = pd.DataFrame(hotspot_regions)
                    hotspot_df = hotspot_df.sort_values('density', ascending=False)
                    
                    # Remove overlaps
                    filtered_hotspots = []
                    used_positions = set()
                    
                    for _, hotspot in hotspot_df.iterrows():
                        overlap = False
                        for used_pos in used_positions:
                            if abs(hotspot['center'] - used_pos) < window_size:
                                overlap = True
                                break
                        
                        if not overlap:
                            filtered_hotspots.append(hotspot.to_dict())
                            used_positions.add(hotspot['center'])
                    
                    chr_hotspots.extend(filtered_hotspots)
            
            hotspots[conflict_type] = chr_hotspots
        
        return hotspots
    
    def analyze_is_ts_association(self) -> Dict:
        """Analyze association between conflict types and IS/TS sites."""
        if self.trc_results is None or self.critical_points is None:
            logger.error("TRC results and/or critical points not loaded")
            return None
        
        logger.info("Analyzing IS/TS association with conflict types...")
        
        associations = {}
        
        # For each gene, find nearest IS/TS site
        enhanced_results = self.trc_results.copy()
        enhanced_results['nearest_is_distance'] = np.nan
        enhanced_results['nearest_ts_distance'] = np.nan
        enhanced_results['nearest_is_rfd'] = np.nan
        enhanced_results['nearest_ts_rfd'] = np.nan
        
        is_sites = self.critical_points[self.critical_points['site_type'] == 'IS']
        ts_sites = self.critical_points[self.critical_points['site_type'] == 'TS']
        
        for idx, gene in enhanced_results.iterrows():
            # Find nearest IS site
            chr_is = is_sites[is_sites['chr'] == gene['chromosome']]
            if len(chr_is) > 0:
                distances = np.abs(chr_is['midpoint'] - gene['midpoint'])
                nearest_idx = distances.idxmin()
                enhanced_results.loc[idx, 'nearest_is_distance'] = distances.min()
                enhanced_results.loc[idx, 'nearest_is_rfd'] = chr_is.loc[nearest_idx, 'rfd_score']
            
            # Find nearest TS site
            chr_ts = ts_sites[ts_sites['chr'] == gene['chromosome']]
            if len(chr_ts) > 0:
                distances = np.abs(chr_ts['midpoint'] - gene['midpoint'])
                nearest_idx = distances.idxmin()
                enhanced_results.loc[idx, 'nearest_ts_distance'] = distances.min()
                enhanced_results.loc[idx, 'nearest_ts_rfd'] = chr_ts.loc[nearest_idx, 'rfd_score']
        
        # Analyze associations
        for conflict_type in ['head-on', 'co-directional', 'neutral']:
            subset = enhanced_results[enhanced_results['conflict_type'] == conflict_type]
            
            if len(subset) == 0:
                continue
            
            associations[conflict_type] = {
                'is_distances': {
                    'mean': subset['nearest_is_distance'].mean(),
                    'median': subset['nearest_is_distance'].median(),
                    'std': subset['nearest_is_distance'].std()
                },
                'ts_distances': {
                    'mean': subset['nearest_ts_distance'].mean(),
                    'median': subset['nearest_ts_distance'].median(),
                    'std': subset['nearest_ts_distance'].std()
                },
                'is_rfd_correlation': subset['rfd_value'].corr(subset['nearest_is_rfd']),
                'ts_rfd_correlation': subset['rfd_value'].corr(subset['nearest_ts_rfd'])
            }
        
        return associations, enhanced_results
    
    def create_peak_distribution_plots(self, distributions: Dict, output_prefix: str):
        """Create comprehensive plots for peak distributions."""
        logger.info("Creating peak distribution plots...")
        
        # Set up the plot style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Transcription-Replication Conflict Analysis: Peak Distributions', fontsize=16, fontweight='bold')
        
        colors = {'head-on': 'red', 'co-directional': 'blue', 'neutral': 'gray'}
        
        # Plot 1: RFD value distributions (histograms)
        ax1 = axes[0, 0]
        for conflict_type, data in distributions.items():
            if 'raw_values' in data:
                ax1.hist(data['raw_values'], bins=20, alpha=0.6, 
                        label=f"{conflict_type} (n={data['statistics']['count']})",
                        color=colors.get(conflict_type, 'black'))
        ax1.set_xlabel('RFD Value')
        ax1.set_ylabel('Frequency')
        ax1.set_title('RFD Value Distributions by Conflict Type')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Box plots of RFD values
        ax2 = axes[0, 1]
        box_data = []
        box_labels = []
        box_colors = []
        for conflict_type, data in distributions.items():
            if 'raw_values' in data:
                box_data.append(data['raw_values'])
                box_labels.append(f"{conflict_type}\n(n={data['statistics']['count']})")
                box_colors.append(colors.get(conflict_type, 'black'))
        
        bp = ax2.boxplot(box_data, labels=box_labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        ax2.set_ylabel('RFD Value')
        ax2.set_title('RFD Value Distributions (Box Plots)')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Peak positions
        ax3 = axes[0, 2]
        x_offset = 0
        for conflict_type, data in distributions.items():
            if 'peaks' in data and len(data['peaks']['positions']) > 0:
                y_pos = [x_offset] * len(data['peaks']['positions'])
                ax3.scatter(data['peaks']['positions'], y_pos, 
                           s=np.array(data['peaks']['heights']) * 10,
                           color=colors.get(conflict_type, 'black'),
                           alpha=0.7, label=conflict_type)
                x_offset += 1
        ax3.set_xlabel('RFD Peak Position')
        ax3.set_ylabel('Conflict Type')
        ax3.set_title('Peak Positions in RFD Distributions')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Statistical summary bar chart
        ax4 = axes[1, 0]
        stats_data = []
        for conflict_type, data in distributions.items():
            stats = data['statistics']
            stats_data.append([stats['mean'], stats['median'], stats['std']])
        
        stats_df = pd.DataFrame(stats_data, 
                               columns=['Mean', 'Median', 'Std Dev'],
                               index=list(distributions.keys()))
        stats_df.plot(kind='bar', ax=ax4, color=['skyblue', 'lightgreen', 'salmon'])
        ax4.set_ylabel('RFD Value')
        ax4.set_title('Statistical Summary by Conflict Type')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45)
        
        # Plot 5: Distribution shape metrics
        ax5 = axes[1, 1]
        shape_data = []
        for conflict_type, data in distributions.items():
            shape = data['shape']
            shape_data.append([shape['skewness'], shape['kurtosis']])
        
        shape_df = pd.DataFrame(shape_data,
                               columns=['Skewness', 'Kurtosis'],
                               index=list(distributions.keys()))
        shape_df.plot(kind='bar', ax=ax5, color=['orange', 'purple'])
        ax5.set_ylabel('Value')
        ax5.set_title('Distribution Shape Metrics')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        plt.setp(ax5.xaxis.get_majorticklabels(), rotation=45)
        
        # Plot 6: Peak count summary
        ax6 = axes[1, 2]
        peak_counts = [data['peaks']['count'] for data in distributions.values()]
        conflict_types = list(distributions.keys())
        bars = ax6.bar(conflict_types, peak_counts, color=[colors.get(ct, 'black') for ct in conflict_types], alpha=0.7)
        ax6.set_ylabel('Number of Peaks')
        ax6.set_title('Peak Count by Conflict Type')
        ax6.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, count in zip(bars, peak_counts):
            ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        plot_file = f"{output_prefix}_peak_distributions.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved peak distribution plot: {plot_file}")
        plt.close()
    
    def create_comprehensive_report(self, distributions: Dict, hotspots: Dict, 
                                   associations: Dict, output_prefix: str):
        """Create a comprehensive text report."""
        report_file = f"{output_prefix}_peak_analysis_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("TRANSCRIPTION-REPLICATION CONFLICT PEAK ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            # Overview
            total_genes = sum([data['statistics']['count'] for data in distributions.values()])
            f.write(f"OVERVIEW:\n")
            f.write(f"Total genes analyzed: {total_genes}\n\n")
            
            # Distribution Analysis
            f.write("DISTRIBUTION ANALYSIS:\n")
            f.write("-" * 40 + "\n")
            
            for conflict_type, data in distributions.items():
                stats = data['statistics']
                shape = data['shape']
                peaks = data['peaks']
                
                f.write(f"\n{conflict_type.upper()} CONFLICTS (n={stats['count']}):\n")
                f.write(f"  RFD Statistics:\n")
                f.write(f"    Mean: {stats['mean']:.4f}\n")
                f.write(f"    Median: {stats['median']:.4f}\n")
                f.write(f"    Std Dev: {stats['std']:.4f}\n")
                f.write(f"    Range: {stats['min']:.4f} to {stats['max']:.4f}\n")
                f.write(f"    IQR: {stats['iqr']:.4f}\n")
                
                f.write(f"  Distribution Shape:\n")
                f.write(f"    Skewness: {shape['skewness']:.4f}")
                if shape['skewness'] > 0.5:
                    f.write(" (right-skewed)")
                elif shape['skewness'] < -0.5:
                    f.write(" (left-skewed)")
                else:
                    f.write(" (approximately symmetric)")
                f.write("\n")
                
                f.write(f"    Kurtosis: {shape['kurtosis']:.4f}")
                if shape['kurtosis'] > 0:
                    f.write(" (heavy-tailed)")
                else:
                    f.write(" (light-tailed)")
                f.write("\n")
                
                f.write(f"  Peak Analysis:\n")
                f.write(f"    Number of peaks: {peaks['count']}\n")
                if peaks['count'] > 0:
                    f.write(f"    Peak positions: {[f'{p:.3f}' for p in peaks['positions']]}\n")
                    f.write(f"    Peak heights: {peaks['heights']}\n")
                f.write("\n")
            
            # Hotspot Analysis
            if hotspots:
                f.write("HOTSPOT ANALYSIS:\n")
                f.write("-" * 40 + "\n")
                
                for conflict_type, hotspot_list in hotspots.items():
                    f.write(f"\n{conflict_type.upper()} HOTSPOTS:\n")
                    
                    if len(hotspot_list) == 0:
                        f.write("  No significant hotspots detected.\n")
                    else:
                        for i, hotspot in enumerate(hotspot_list[:5]):  # Top 5 hotspots
                            f.write(f"  Hotspot {i+1}:\n")
                            f.write(f"    Location: {hotspot['chromosome']}:{hotspot['start']}-{hotspot['end']}\n")
                            f.write(f"    Gene count: {hotspot['gene_count']}\n")
                            f.write(f"    Density: {hotspot['density']:.2f} genes/Mb\n")
                            f.write(f"    Mean RFD: {hotspot['mean_rfd']:.4f}\n")
                            f.write(f"    Genes: {', '.join(hotspot['genes'][:10])}")
                            if len(hotspot['genes']) > 10:
                                f.write(f" ... (+{len(hotspot['genes'])-10} more)")
                            f.write("\n\n")
            
            # IS/TS Association Analysis
            if associations:
                f.write("IS/TS SITE ASSOCIATION ANALYSIS:\n")
                f.write("-" * 40 + "\n")
                
                for conflict_type, assoc_data in associations.items():
                    f.write(f"\n{conflict_type.upper()} ASSOCIATIONS:\n")
                    f.write(f"  Distance to nearest IS site:\n")
                    f.write(f"    Mean: {assoc_data['is_distances']['mean']:.0f} bp\n")
                    f.write(f"    Median: {assoc_data['is_distances']['median']:.0f} bp\n")
                    
                    f.write(f"  Distance to nearest TS site:\n")
                    f.write(f"    Mean: {assoc_data['ts_distances']['mean']:.0f} bp\n")
                    f.write(f"    Median: {assoc_data['ts_distances']['median']:.0f} bp\n")
                    
                    f.write(f"  RFD Correlations:\n")
                    f.write(f"    With IS sites: {assoc_data['is_rfd_correlation']:.4f}\n")
                    f.write(f"    With TS sites: {assoc_data['ts_rfd_correlation']:.4f}\n\n")
            
            # Biological Interpretation
            f.write("BIOLOGICAL INTERPRETATION:\n")
            f.write("-" * 40 + "\n")
            
            # Head-on analysis
            if 'head-on' in distributions:
                ho_data = distributions['head-on']
                f.write("\nHEAD-ON CONFLICTS:\n")
                f.write("- Associated with increased replication stress and mutagenesis\n")
                f.write(f"- {ho_data['statistics']['count']} genes ({ho_data['statistics']['count']/total_genes*100:.1f}% of total)\n")
                f.write(f"- Mean RFD: {ho_data['statistics']['mean']:.3f}\n")
                
                if ho_data['peaks']['count'] > 1:
                    f.write(f"- Multiple peaks detected ({ho_data['peaks']['count']}), suggesting subgroups\n")
                elif ho_data['peaks']['count'] == 1:
                    f.write("- Single dominant peak, indicating consistent RFD pattern\n")
                else:
                    f.write("- No clear peaks, indicating dispersed RFD distribution\n")
            
            # Co-directional analysis  
            if 'co-directional' in distributions:
                cd_data = distributions['co-directional']
                f.write("\nCO-DIRECTIONAL ORIENTATION:\n")
                f.write("- Generally favorable for replication progression\n")
                f.write(f"- {cd_data['statistics']['count']} genes ({cd_data['statistics']['count']/total_genes*100:.1f}% of total)\n")
                f.write(f"- Mean RFD: {cd_data['statistics']['mean']:.3f}\n")
                
                if cd_data['peaks']['count'] > 1:
                    f.write(f"- Multiple peaks detected ({cd_data['peaks']['count']}), suggesting context-dependent effects\n")
            
            # Neutral analysis
            if 'neutral' in distributions:
                n_data = distributions['neutral']
                f.write("\nNEUTRAL ZONES:\n")
                f.write("- Regions with low RFD signal, potential replication timing boundaries\n")
                f.write(f"- {n_data['statistics']['count']} genes ({n_data['statistics']['count']/total_genes*100:.1f}% of total)\n")
                f.write(f"- Mean RFD: {n_data['statistics']['mean']:.3f}\n")
            
            f.write("\n" + "=" * 80 + "\n")
        
        logger.info(f"Comprehensive report saved: {report_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze peaks and distributions in transcription-replication conflicts",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--trc-results', required=True,
                       help='TRC analysis results TSV file')
    parser.add_argument('--critical-points', required=True,
                       help='Critical points BED file (IS/TS sites)')
    parser.add_argument('--output', required=True,
                       help='Output prefix for results')
    parser.add_argument('--create-plots', action='store_true',
                       help='Create distribution plots (requires matplotlib)')
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = TRCPeakAnalyzer()
    
    # Load data
    trc_results = analyzer.load_trc_results(args.trc_results)
    if trc_results is None:
        sys.exit(1)
    
    critical_points = analyzer.load_critical_points(args.critical_points)
    if critical_points is None:
        logger.warning("Could not load critical points - IS/TS analysis will be skipped")
    
    # Analyze distributions
    distributions = analyzer.analyze_conflict_distributions()
    if distributions is None:
        sys.exit(1)
    
    # Find hotspots
    hotspots = analyzer.find_conflict_hotspots()
    
    # Analyze IS/TS associations
    associations = None
    enhanced_results = None
    if critical_points is not None:
        associations, enhanced_results = analyzer.analyze_is_ts_association()
    
    # Create plots if requested
    if args.create_plots:
        try:
            analyzer.create_peak_distribution_plots(distributions, args.output)
        except Exception as e:
            logger.error(f"Error creating plots: {e}")
    
    # Create comprehensive report
    analyzer.create_comprehensive_report(distributions, hotspots, associations, args.output)
    
    # Export enhanced results with IS/TS distances
    if enhanced_results is not None:
        enhanced_file = f"{args.output}_enhanced_trc_results.tsv"
        enhanced_results.to_csv(enhanced_file, sep='\t', index=False)
        logger.info(f"Enhanced results with IS/TS distances saved: {enhanced_file}")
    
    logger.info("Peak distribution analysis completed successfully!")

if __name__ == "__main__":
    main()