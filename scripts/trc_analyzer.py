#!/usr/bin/env python3
"""
Transcription-Replication Conflict (TRC) Analysis Tool

This script analyzes transcription-replication conflicts by comparing:
1. Gene coordinates and transcription directions
2. Replication fork directionality (RFD) from OK-seq data
3. Identifies head-on and co-directional conflicts

Conflict Types:
- Head-on: Gene transcription direction opposite to replication fork direction
- Co-directional: Gene transcription direction same as replication fork direction

Input: 
- Gene coordinate file (BED/GTF/GFF format)
- RFD bedgraph file or RFD analysis results
Output: Conflict classification and statistics
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
import gzip
from typing import List, Tuple, Dict, Optional
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TRCAnalyzer:
    def __init__(self):
        """Initialize the Transcription-Replication Conflict analyzer."""
        self.genes = None
        self.rfd_data = None
        self.zones = None
        self.conflicts = None
        
    def load_gene_coordinates(self, filepath: str, file_format: str = "auto") -> pd.DataFrame:
        """
        Load gene coordinates from various formats.
        
        Args:
            filepath: Path to gene coordinate file
            file_format: Format of input file ("bed", "gtf", "gff", "auto")
            
        Returns:
            DataFrame with standardized gene coordinates
        """
        logger.info(f"Loading gene coordinates from: {filepath}")
        
        # Auto-detect format from file extension
        if file_format == "auto":
            ext = filepath.lower().split('.')[-1]
            if ext in ['bed']:
                file_format = "bed"
            elif ext in ['gtf', 'gff', 'gff3']:
                file_format = "gtf"
            else:
                logger.warning(f"Unknown format for {filepath}, assuming BED format")
                file_format = "bed"
        
        try:
            if file_format == "bed":
                # BED format: chr, start, end, name, score, strand, ...
                df = pd.read_csv(filepath, sep='\t', header=None, comment='#')
                
                # Ensure we have at least 6 columns for strand information
                if df.shape[1] < 6:
                    logger.error("BED file must have at least 6 columns (chr, start, end, name, score, strand)")
                    return None
                
                df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand'] + [f'col_{i}' for i in range(6, df.shape[1])]
                
            elif file_format in ["gtf", "gff"]:
                # GTF/GFF format: chr, source, feature, start, end, score, strand, frame, attributes
                df = pd.read_csv(filepath, sep='\t', header=None, comment='#')
                df.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
                
                # Filter for genes only (not exons, CDS, etc.)
                df = df[df['feature'].isin(['gene', 'transcript', 'mRNA'])]
                
                # Extract gene name from attributes
                def extract_gene_name(attrs):
                    if pd.isna(attrs):
                        return "unknown"
                    # Try different attribute patterns
                    for pattern in ['gene_id "([^"]+)"', 'gene_name "([^"]+)"', 'Name=([^;]+)', 'ID=([^;]+)']:
                        import re
                        match = re.search(pattern, str(attrs))
                        if match:
                            return match.group(1)
                    return "unknown"
                
                df['name'] = df['attributes'].apply(extract_gene_name)
            
            # Standardize chromosome names (remove 'chr' prefix if present for consistency)
            df['chr'] = df['chr'].astype(str)
            
            # Filter out invalid strands and convert to numeric
            df = df[df['strand'].isin(['+', '-', '1', '-1'])]
            df['strand_numeric'] = df['strand'].map({'+': 1, '-': -1, '1': 1, '-1': -1})
            
            # Ensure coordinates are integers
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            
            # Calculate gene length and midpoint
            df['length'] = df['end'] - df['start']
            df['midpoint'] = (df['start'] + df['end']) // 2
            
            logger.info(f"Loaded {len(df)} genes from {len(df['chr'].unique())} chromosomes")
            logger.info(f"Strand distribution: {df['strand'].value_counts().to_dict()}")
            
            self.genes = df
            return df
            
        except Exception as e:
            logger.error(f"Error loading gene coordinates: {e}")
            return None
    
    def load_rfd_data(self, filepath: str) -> pd.DataFrame:
        """
        Load RFD data from bedgraph file.
        
        Args:
            filepath: Path to RFD bedgraph file (.bedgraph or .bedgraph.gz)
            
        Returns:
            DataFrame with RFD data
        """
        logger.info(f"Loading RFD data from: {filepath}")
        
        try:
            # Handle gzipped files
            if filepath.endswith('.gz'):
                with gzip.open(filepath, 'rt') as f:
                    df = pd.read_csv(f, sep='\t', header=None, comment='#')
            else:
                df = pd.read_csv(filepath, sep='\t', header=None, comment='#')
            
            # Standard bedgraph format: chr, start, end, value
            df.columns = ['chr', 'start', 'end', 'rfd_value']
            
            # Standardize chromosome names
            df['chr'] = df['chr'].astype(str)
            
            # Calculate bin midpoint and width
            df['midpoint'] = (df['start'] + df['end']) // 2
            df['bin_width'] = df['end'] - df['start']
            
            logger.info(f"Loaded {len(df)} RFD intervals from {len(df['chr'].unique())} chromosomes")
            logger.info(f"RFD value range: {df['rfd_value'].min():.3f} to {df['rfd_value'].max():.3f}")
            
            self.rfd_data = df
            return df
            
        except Exception as e:
            logger.error(f"Error loading RFD data: {e}")
            return None
    
    def load_rfd_zones(self, critical_points_file: str, zones_file: Optional[str] = None):
        """
        Load RFD analysis results (initiation/termination sites and zones).
        
        Args:
            critical_points_file: Path to critical points BED file
            zones_file: Optional path to zones BED file
        """
        logger.info("Loading RFD analysis results...")
        
        zones_data = {}
        
        # Load critical points (IS/TS)
        try:
            critical_df = pd.read_csv(critical_points_file, sep='\t', header=None)
            critical_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end', 'color']
            
            # Extract site type from name (IS or TS)
            critical_df['site_type'] = critical_df['name'].str.extract(r'(IS|TS)_')
            critical_df['midpoint'] = (critical_df['start'] + critical_df['end']) // 2
            
            zones_data['critical_points'] = critical_df
            logger.info(f"Loaded {len(critical_df)} critical points")
            
        except Exception as e:
            logger.error(f"Error loading critical points: {e}")
            return None
        
        # Load zones if provided
        if zones_file and os.path.exists(zones_file):
            try:
                zones_df = pd.read_csv(zones_file, sep='\t', header=None)
                zones_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thick_start', 'thick_end', 'color']
                
                # Extract zone type from name (IZ or TZ)
                zones_df['zone_type'] = zones_df['name'].str.extract(r'(IZ|TZ)_')
                zones_df['midpoint'] = (zones_df['start'] + zones_df['end']) // 2
                zones_df['length'] = zones_df['end'] - zones_df['start']
                
                zones_data['zones'] = zones_df
                logger.info(f"Loaded {len(zones_df)} replication zones")
                
            except Exception as e:
                logger.warning(f"Could not load zones file: {e}")
        
        self.zones = zones_data
        return zones_data
    
    def interpolate_rfd_at_position(self, chromosome: str, position: int) -> Optional[float]:
        """
        Interpolate RFD value at a specific genomic position.
        
        Args:
            chromosome: Chromosome name
            position: Genomic position
            
        Returns:
            Interpolated RFD value or None if not available
        """
        if self.rfd_data is None:
            return None
        
        # Filter RFD data for the specific chromosome
        chr_data = self.rfd_data[self.rfd_data['chr'] == chromosome].copy()
        
        if len(chr_data) == 0:
            return None
        
        # Sort by position
        chr_data = chr_data.sort_values('start')
        
        # Find overlapping interval
        overlapping = chr_data[(chr_data['start'] <= position) & (chr_data['end'] > position)]
        
        if len(overlapping) > 0:
            # Direct overlap - return the RFD value
            return overlapping.iloc[0]['rfd_value']
        
        # No direct overlap - interpolate between nearest intervals
        before = chr_data[chr_data['end'] <= position]
        after = chr_data[chr_data['start'] > position]
        
        if len(before) > 0 and len(after) > 0:
            # Linear interpolation between nearest points
            before_interval = before.iloc[-1]  # Last interval before position
            after_interval = after.iloc[0]     # First interval after position
            
            # Use midpoints for interpolation
            x1, y1 = before_interval['midpoint'], before_interval['rfd_value']
            x2, y2 = after_interval['midpoint'], after_interval['rfd_value']
            
            # Linear interpolation
            if x2 != x1:
                interpolated_rfd = y1 + (y2 - y1) * (position - x1) / (x2 - x1)
                return interpolated_rfd
        
        # Fallback: return nearest value
        if len(before) > 0:
            return before.iloc[-1]['rfd_value']
        elif len(after) > 0:
            return after.iloc[0]['rfd_value']
        
        return None
    
    def determine_replication_direction(self, rfd_value: float) -> int:
        """
        Determine replication fork direction from RFD value.
        
        Args:
            rfd_value: RFD value
            
        Returns:
            1 for rightward replication, -1 for leftward replication, 0 for neutral
        """
        # RFD > 0: rightward replication fork movement
        # RFD < 0: leftward replication fork movement
        # RFD ≈ 0: neutral or conflict zone
        
        threshold = 0.05  # Neutral zone threshold
        
        if rfd_value > threshold:
            return 1   # Rightward replication
        elif rfd_value < -threshold:
            return -1  # Leftward replication
        else:
            return 0   # Neutral/conflict zone
    
    def classify_transcription_replication_conflict(self, gene_strand: int, replication_direction: int) -> str:
        """
        Classify transcription-replication conflict type.
        
        Args:
            gene_strand: Gene transcription direction (1 for +, -1 for -)
            replication_direction: Replication direction (1 for rightward, -1 for leftward)
            
        Returns:
            Conflict type: "head-on", "co-directional", "neutral", or "unknown"
        """
        if replication_direction == 0:
            return "neutral"
        
        if gene_strand == replication_direction:
            return "co-directional"
        elif gene_strand == -replication_direction:
            return "head-on"
        else:
            return "unknown"
    
    def analyze_conflicts(self) -> pd.DataFrame:
        """
        Analyze transcription-replication conflicts for all genes.
        
        Returns:
            DataFrame with conflict analysis results
        """
        if self.genes is None:
            logger.error("No gene data loaded. Use load_gene_coordinates() first.")
            return None
        
        if self.rfd_data is None:
            logger.error("No RFD data loaded. Use load_rfd_data() first.")
            return None
        
        logger.info("Analyzing transcription-replication conflicts...")
        
        results = []
        
        for idx, gene in self.genes.iterrows():
            # Get RFD value at gene midpoint
            rfd_value = self.interpolate_rfd_at_position(gene['chr'], gene['midpoint'])
            
            if rfd_value is not None:
                # Determine replication direction
                replication_direction = self.determine_replication_direction(rfd_value)
                
                # Classify conflict type
                conflict_type = self.classify_transcription_replication_conflict(
                    gene['strand_numeric'], replication_direction
                )
                
                # Additional analysis: check if gene overlaps with replication zones
                zone_overlap = self.check_zone_overlap(gene['chr'], gene['start'], gene['end'])
                
                results.append({
                    'gene_name': gene.get('name', f"gene_{idx}"),
                    'chromosome': gene['chr'],
                    'start': gene['start'],
                    'end': gene['end'],
                    'midpoint': gene['midpoint'],
                    'length': gene['length'],
                    'gene_strand': gene['strand'],
                    'gene_strand_numeric': gene['strand_numeric'],
                    'rfd_value': rfd_value,
                    'replication_direction': replication_direction,
                    'conflict_type': conflict_type,
                    'zone_type': zone_overlap.get('zone_type', 'none'),
                    'zone_distance': zone_overlap.get('distance', np.nan)
                })
            else:
                # No RFD data available for this gene
                results.append({
                    'gene_name': gene.get('name', f"gene_{idx}"),
                    'chromosome': gene['chr'],
                    'start': gene['start'],
                    'end': gene['end'],
                    'midpoint': gene['midpoint'],
                    'length': gene['length'],
                    'gene_strand': gene['strand'],
                    'gene_strand_numeric': gene['strand_numeric'],
                    'rfd_value': np.nan,
                    'replication_direction': 0,
                    'conflict_type': 'no_data',
                    'zone_type': 'none',
                    'zone_distance': np.nan
                })
        
        conflicts_df = pd.DataFrame(results)
        
        # Generate summary statistics
        self.generate_conflict_summary(conflicts_df)
        
        self.conflicts = conflicts_df
        return conflicts_df
    
    def check_zone_overlap(self, chromosome: str, start: int, end: int) -> Dict:
        """
        Check if a genomic region overlaps with replication zones.
        
        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
            
        Returns:
            Dictionary with overlap information
        """
        if self.zones is None or 'zones' not in self.zones:
            return {'zone_type': 'none', 'distance': np.nan}
        
        zones_df = self.zones['zones']
        chr_zones = zones_df[zones_df['chr'] == chromosome]
        
        if len(chr_zones) == 0:
            return {'zone_type': 'none', 'distance': np.nan}
        
        # Check for direct overlap
        overlapping = chr_zones[
            (chr_zones['start'] <= end) & (chr_zones['end'] >= start)
        ]
        
        if len(overlapping) > 0:
            # Return the zone type with maximum overlap
            overlaps = []
            for _, zone in overlapping.iterrows():
                overlap_start = max(start, zone['start'])
                overlap_end = min(end, zone['end'])
                overlap_length = overlap_end - overlap_start
                overlaps.append((zone['zone_type'], overlap_length))
            
            # Sort by overlap length and return the largest
            overlaps.sort(key=lambda x: x[1], reverse=True)
            return {'zone_type': overlaps[0][0], 'distance': 0}
        
        # No direct overlap - find nearest zone
        gene_midpoint = (start + end) // 2
        chr_zones['distance'] = np.abs(chr_zones['midpoint'] - gene_midpoint)
        nearest = chr_zones.loc[chr_zones['distance'].idxmin()]
        
        return {
            'zone_type': f"near_{nearest['zone_type']}", 
            'distance': nearest['distance']
        }
    
    def generate_conflict_summary(self, conflicts_df: pd.DataFrame):
        """
        Generate and display summary statistics for TRC analysis.
        
        Args:
            conflicts_df: DataFrame with conflict analysis results
        """
        logger.info("=== TRANSCRIPTION-REPLICATION CONFLICT ANALYSIS SUMMARY ===")
        
        total_genes = len(conflicts_df)
        logger.info(f"Total genes analyzed: {total_genes}")
        
        # Conflict type distribution
        conflict_counts = conflicts_df['conflict_type'].value_counts()
        logger.info("\nConflict Type Distribution:")
        for conflict_type, count in conflict_counts.items():
            percentage = (count / total_genes) * 100
            logger.info(f"  {conflict_type}: {count} ({percentage:.1f}%)")
        
        # Chromosome distribution
        chr_distribution = conflicts_df.groupby(['chromosome', 'conflict_type']).size().unstack(fill_value=0)
        logger.info(f"\nConflicts by chromosome:")
        for chromosome in chr_distribution.index:
            head_on = chr_distribution.loc[chromosome].get('head-on', 0)
            co_dir = chr_distribution.loc[chromosome].get('co-directional', 0)
            total_chr = head_on + co_dir
            if total_chr > 0:
                ho_pct = (head_on / total_chr) * 100
                logger.info(f"  {chromosome}: {head_on} head-on ({ho_pct:.1f}%), {co_dir} co-directional")
        
        # Zone overlap analysis
        if 'zone_type' in conflicts_df.columns:
            zone_counts = conflicts_df['zone_type'].value_counts()
            logger.info("\nReplication Zone Overlap:")
            for zone_type, count in zone_counts.items():
                percentage = (count / total_genes) * 100
                logger.info(f"  {zone_type}: {count} ({percentage:.1f}%)")
        
        # RFD value statistics by conflict type
        rfd_stats = conflicts_df.groupby('conflict_type')['rfd_value'].describe()
        logger.info("\nRFD Value Statistics by Conflict Type:")
        logger.info(rfd_stats)
    
    def export_results(self, output_prefix: str, output_formats: List[str] = ['tsv', 'bed']):
        """
        Export TRC analysis results to various formats.
        
        Args:
            output_prefix: Output file prefix
            output_formats: List of output formats ('tsv', 'bed', 'gff')
        """
        if self.conflicts is None:
            logger.error("No conflict analysis results to export. Run analyze_conflicts() first.")
            return
        
        logger.info(f"Exporting results with prefix: {output_prefix}")
        
        if 'tsv' in output_formats:
            # Export detailed TSV file
            tsv_file = f"{output_prefix}_trc_analysis.tsv"
            self.conflicts.to_csv(tsv_file, sep='\t', index=False)
            logger.info(f"Exported detailed results to: {tsv_file}")
        
        if 'bed' in output_formats:
            # Export BED files for each conflict type
            for conflict_type in self.conflicts['conflict_type'].unique():
                if conflict_type in ['head-on', 'co-directional']:
                    subset = self.conflicts[self.conflicts['conflict_type'] == conflict_type]
                    
                    bed_file = f"{output_prefix}_{conflict_type.replace('-', '_')}_genes.bed"
                    
                    # Create BED format with color coding
                    bed_data = subset[['chromosome', 'start', 'end', 'gene_name', 'rfd_value', 'gene_strand']].copy()
                    bed_data['score'] = 1000  # Max BED score
                    bed_data['thick_start'] = bed_data['start']
                    bed_data['thick_end'] = bed_data['end']
                    
                    # Color coding: + strand = red, - strand = blue
                    bed_data['color'] = bed_data['gene_strand'].apply(lambda x: '255,0,0' if x == '+' else '0,0,255')
                    
                    # Reorder columns for BED format
                    bed_data = bed_data[['chromosome', 'start', 'end', 'gene_name', 'score', 'gene_strand', 
                                       'thick_start', 'thick_end', 'color']]
                    
                    bed_data.to_csv(bed_file, sep='\t', index=False, header=False)
                    logger.info(f"Exported {conflict_type} genes to: {bed_file}")
        
        if 'gff' in output_formats:
            # Export GFF format
            gff_file = f"{output_prefix}_trc_analysis.gff"
            
            with open(gff_file, 'w') as f:
                f.write("##gff-version 3\n")
                
                for _, gene in self.conflicts.iterrows():
                    if gene['conflict_type'] in ['head-on', 'co-directional']:
                        # GFF format: seqname, source, feature, start, end, score, strand, frame, attribute
                        f.write(f"{gene['chromosome']}\t")
                        f.write(f"TRC_analysis\t")
                        f.write(f"gene\t")
                        f.write(f"{gene['start']}\t")
                        f.write(f"{gene['end']}\t")
                        f.write(f"{gene['rfd_value']:.3f}\t")
                        f.write(f"{gene['gene_strand']}\t")
                        f.write(f".\t")
                        f.write(f"ID={gene['gene_name']};conflict_type={gene['conflict_type']};")
                        f.write(f"rfd_value={gene['rfd_value']:.3f}\n")
            
            logger.info(f"Exported GFF results to: {gff_file}")

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Analyze transcription-replication conflicts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with BED gene file and RFD bedgraph
  python trc_analyzer.py -g genes.bed -r rfd_data.bedgraph -o output_prefix
  
  # Analysis with GTF gene file and RFD analysis results
  python trc_analyzer.py -g genes.gtf -c critical_points.bed -z zones.bed -o output_prefix
  
  # Include additional output formats
  python trc_analyzer.py -g genes.bed -r rfd_data.bedgraph -o output_prefix --formats tsv bed gff
        """
    )
    
    # Input files
    parser.add_argument('-g', '--genes', required=True,
                       help='Gene coordinates file (BED, GTF, or GFF format)')
    parser.add_argument('-r', '--rfd-data', 
                       help='RFD bedgraph file (.bedgraph or .bedgraph.gz)')
    parser.add_argument('-c', '--critical-points',
                       help='RFD critical points BED file (IS/TS sites)')
    parser.add_argument('-z', '--zones',
                       help='RFD zones BED file (IZ/TZ regions)')
    
    # Output options
    parser.add_argument('-o', '--output', required=True,
                       help='Output file prefix')
    parser.add_argument('--formats', nargs='+', default=['tsv', 'bed'],
                       choices=['tsv', 'bed', 'gff'],
                       help='Output formats (default: tsv bed)')
    
    # Analysis options
    parser.add_argument('--gene-format', choices=['auto', 'bed', 'gtf', 'gff'], default='auto',
                       help='Gene file format (default: auto-detect)')
    parser.add_argument('--rfd-threshold', type=float, default=0.05,
                       help='RFD threshold for neutral zone (default: 0.05)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not args.rfd_data and not args.critical_points:
        logger.error("Either --rfd-data or --critical-points must be provided")
        sys.exit(1)
    
    # Initialize analyzer
    analyzer = TRCAnalyzer()
    
    # Load gene coordinates
    genes_df = analyzer.load_gene_coordinates(args.genes, args.gene_format)
    if genes_df is None:
        logger.error("Failed to load gene coordinates")
        sys.exit(1)
    
    # Load RFD data
    if args.rfd_data:
        rfd_df = analyzer.load_rfd_data(args.rfd_data)
        if rfd_df is None:
            logger.error("Failed to load RFD data")
            sys.exit(1)
    
    # Load RFD analysis results
    if args.critical_points:
        zones = analyzer.load_rfd_zones(args.critical_points, args.zones)
        if zones is None:
            logger.error("Failed to load RFD analysis results")
            sys.exit(1)
    
    # Analyze conflicts
    conflicts_df = analyzer.analyze_conflicts()
    if conflicts_df is None:
        logger.error("Failed to analyze conflicts")
        sys.exit(1)
    
    # Export results
    analyzer.export_results(args.output, args.formats)
    
    logger.info("TRC analysis completed successfully!")

if __name__ == "__main__":
    main()