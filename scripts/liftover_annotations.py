#!/usr/bin/env python3
"""
Liftover script for GLOseq critical points and zones.

This script lifts over genomic coordinates from one assembly to another using 
pyliftover and chain files. Supports conversions between hg19, hg38, and hs1.

Usage:
    python liftover_annotations.py input.bed --from hg19 --to hg38 --output output.bed
    python liftover_annotations.py critical_points.bed --from hg38 --to hs1 --chain-file misc/hg38ToHs1.over.chain.gz

Author: GLOseq Pipeline
"""

import argparse
import pandas as pd
import os
import sys
from pathlib import Path
import logging
import gzip
from typing import Optional, Tuple, Dict, List

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try:
    from pyliftover import LiftOver
    PYLIFTOVER_AVAILABLE = True
except ImportError:
    PYLIFTOVER_AVAILABLE = False
    logger.warning("pyliftover not installed. Install with: pip install pyliftover")

class GenomeLiftOver:
    """Class for lifting over genomic coordinates between assemblies."""
    
    # Default chain file locations
    CHAIN_FILES = {
        ('hg19', 'hg38'): 'misc/hg19ToHg38.over.chain.gz',
        ('hg38', 'hg19'): 'misc/hg38ToHg19.over.chain.gz',
        ('hg38', 'hs1'): 'misc/hg38ToHs1.over.chain.gz',
        ('hs1', 'hg38'): 'misc/hs1ToHg38.over.chain.gz',
        ('hg19', 'hs1'): 'misc/hg19ToHs1.over.chain.gz',
        ('hs1', 'hg19'): 'misc/hs1ToHg19.over.chain.gz',
    }
    
    def __init__(self, from_assembly: str, to_assembly: str, chain_file: str = None):
        """
        Initialize the liftover object.
        
        Args:
            from_assembly: Source genome assembly (e.g., 'hg19', 'hg38', 'hs1')
            to_assembly: Target genome assembly (e.g., 'hg19', 'hg38', 'hs1')
            chain_file: Path to chain file (optional, will use default if available)
        """
        self.from_assembly = from_assembly.lower()
        self.to_assembly = to_assembly.lower()
        
        if not PYLIFTOVER_AVAILABLE:
            raise ImportError("pyliftover is required. Install with: pip install pyliftover")
        
        # Determine chain file
        if chain_file:
            self.chain_file = chain_file
        else:
            key = (self.from_assembly, self.to_assembly)
            if key in self.CHAIN_FILES:
                self.chain_file = self.CHAIN_FILES[key]
            else:
                raise ValueError(f"No default chain file for {from_assembly} -> {to_assembly}")
        
        # Check if chain file exists
        if not os.path.exists(self.chain_file):
            raise FileNotFoundError(f"Chain file not found: {self.chain_file}")
        
        # Initialize LiftOver
        logger.info(f"Loading chain file: {self.chain_file}")
        self.liftover = LiftOver(self.chain_file)
        logger.info(f"Initialized liftover: {from_assembly} -> {to_assembly}")
    
    def liftover_position(self, chromosome: str, position: int) -> Optional[Tuple[str, int]]:
        """
        Lift over a single genomic position.
        
        Args:
            chromosome: Chromosome name (e.g., 'chr1', '1')
            position: Genomic position (0-based or 1-based)
            
        Returns:
            Tuple of (new_chromosome, new_position) or None if liftover fails
        """
        # Ensure chromosome has 'chr' prefix for liftover
        if not chromosome.startswith('chr'):
            chromosome = f'chr{chromosome}'
        
        try:
            result = self.liftover.convert_coordinate(chromosome, position)
            if result:
                return result[0]  # Return first result (chr, pos)
            else:
                return None
        except Exception as e:
            logger.warning(f"Liftover failed for {chromosome}:{position} - {e}")
            return None
    
    def liftover_interval(self, chromosome: str, start: int, end: int) -> Optional[Tuple[str, int, int]]:
        """
        Lift over a genomic interval.
        
        Args:
            chromosome: Chromosome name
            start: Start position
            end: End position
            
        Returns:
            Tuple of (new_chromosome, new_start, new_end) or None if liftover fails
        """
        start_result = self.liftover_position(chromosome, start)
        end_result = self.liftover_position(chromosome, end)
        
        if start_result and end_result:
            # Check if both coordinates map to the same chromosome
            if start_result[0] == end_result[0]:
                new_chr = start_result[0]
                new_start = min(start_result[1], end_result[1])
                new_end = max(start_result[1], end_result[1])
                return (new_chr, new_start, new_end)
            else:
                logger.warning(f"Start and end coordinates map to different chromosomes: "
                             f"{start_result[0]} vs {end_result[0]}")
                return None
        else:
            return None

def detect_bed_format(filepath: str) -> str:
    """
    Detect the format of a BED file (BED3, BED4, BED9, etc.)
    
    Args:
        filepath: Path to BED file
        
    Returns:
        Format description string
    """
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        if not first_line or first_line.startswith('#'):
            # Skip header lines
            for line in f:
                if line.strip() and not line.startswith('#'):
                    first_line = line.strip()
                    break
    
    columns = first_line.split('\t')
    num_cols = len(columns)
    
    if num_cols >= 9:
        return "BED9"
    elif num_cols >= 6:
        return "BED6"
    elif num_cols >= 4:
        return "BED4"
    elif num_cols >= 3:
        return "BED3"
    else:
        return "Unknown"

def liftover_bed_file(input_file: str, output_file: str, liftover_obj: GenomeLiftOver, 
                     keep_failed: bool = False) -> Dict[str, int]:
    
    """
    Lift over coordinates in a BED file.
    
    Args:
        input_file: Input BED file path
        output_file: Output BED file path
        liftover_obj: GenomeLiftOver object
        keep_failed: If True, keep intervals that failed liftover with original coordinates
        
    Returns:
        Dictionary with statistics
    """
    logger.info(f"Processing BED file: {input_file}")
    bed_format = detect_bed_format(input_file)
    logger.info(f"Detected format: {bed_format}")
    
    # Read input file
    df = pd.read_csv(input_file, sep='\t', header=None, comment='#')
    
    # Ensure we have at least 3 columns
    if df.shape[1] < 3:
        raise ValueError("BED file must have at least 3 columns (chr, start, end)")
    
    # Column names based on format
    if df.shape[1] >= 9:
        column_names = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']
        df.columns = column_names[:df.shape[1]]
    else:
        basic_names = ['chr', 'start', 'end', 'name', 'score', 'strand']
        df.columns = basic_names[:df.shape[1]]
    
    # Track statistics
    stats = {'total': len(df), 'successful': 0, 'failed': 0}
    
    # Lists to store results
    successful_rows = []
    failed_rows = []
    
    logger.info(f"Lifting over {len(df)} intervals...")
    output_dir = os.path.dirname(output_file)
    logger.info(f"Output directory: {output_dir}")
    if output_dir and not os.path.exists(output_dir):
        logger.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    for idx, row in df.iterrows():
        chr_name = str(row['chr'])
        start_pos = int(row['start'])
        end_pos = int(row['end'])
        
        # Perform liftover
        result = liftover_obj.liftover_interval(chr_name, start_pos, end_pos)
        
        if result:
            new_chr, new_start, new_end = result
            # Update coordinates
            new_row = row.copy()
            new_row['chr'] = new_chr
            new_row['start'] = new_start
            new_row['end'] = new_end
            
            # Update thickStart and thickEnd if they exist
            if 'thickStart' in new_row:
                new_row['thickStart'] = new_start
            if 'thickEnd' in new_row:
                new_row['thickEnd'] = new_end
                
            successful_rows.append(new_row)
            stats['successful'] += 1
        else:
            if keep_failed:
                failed_rows.append(row)
            stats['failed'] += 1
    
    # Combine results
    if successful_rows:
        result_df = pd.DataFrame(successful_rows)
        if keep_failed and failed_rows:
            failed_df = pd.DataFrame(failed_rows)
            result_df = pd.concat([result_df, failed_df], ignore_index=True)
        
        # Sort by chromosome and position
        result_df = result_df.sort_values(['chr', 'start'])
        
        # Save output
        
        result_df.to_csv(output_file, sep='\t', index=False, header=False)
        logger.info(f"Saved {len(result_df)} intervals to {output_file}")
    else:
        logger.warning("No intervals were successfully lifted over!")
        if keep_failed and failed_rows:
            failed_df = pd.DataFrame(failed_rows)
            failed_df.to_csv(output_file, sep='\t', index=False, header=False)
            logger.info(f"Saved {len(failed_df)} original intervals (liftover failed) to {output_file}")
    
    return stats

def main():
    parser = argparse.ArgumentParser(
        description="Liftover genomic coordinates in BED files between genome assemblies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Liftover critical points from hg19 to hg38
  python liftover_annotations.py critical_points_hg19.bed --from hg19 --to hg38 --output critical_points_hg38.bed
  
  # Liftover zones from hg38 to hs1 with custom chain file
  python liftover_annotations.py zones_hg38.bed --from hg38 --to hs1 --chain-file misc/hg38ToHs1.over.chain.gz --output zones_hs1.bed
  
  # Batch process multiple files
  python liftover_annotations.py results/*.bed --from hg19 --to hg38 --output-dir results_hg38/

Supported assemblies: hg19, hg38, hs1
Chain files should be placed in the misc/ directory.
        """
    )
    
    parser.add_argument('input_files', nargs='+', 
                       help='Input BED file(s) to liftover')
    parser.add_argument('--from', dest='from_assembly', required=True,
                       choices=['hg19', 'hg38', 'hs1'],
                       help='Source genome assembly')
    parser.add_argument('--to', dest='to_assembly', required=True,
                       choices=['hg19', 'hg38', 'hs1'],
                       help='Target genome assembly')
    parser.add_argument('--output', '-o',
                       help='Output file (for single input file)')
    parser.add_argument('--output-dir',
                       help='Output directory (for multiple input files)')
    parser.add_argument('--chain-file',
                       help='Path to chain file (optional, uses default if available)')
    parser.add_argument('--keep-failed', action='store_true',
                       help='Keep intervals that failed liftover with original coordinates')
    parser.add_argument('--suffix', default='_lifted',
                       help='Suffix for output files when using --output-dir (default: _lifted)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if len(args.input_files) == 1 and args.output:
        # Single file mode
        output_mode = 'single'
    elif len(args.input_files) >= 1 and args.output_dir:
        # Multiple files mode
        output_mode = 'multiple'
        os.makedirs(args.output_dir, exist_ok=True)
    else:
        logger.error("Either specify --output for single file or --output-dir for multiple files")
        sys.exit(1)
    
    if args.from_assembly == args.to_assembly:
        logger.error("Source and target assemblies cannot be the same")
        sys.exit(1)
    
    # Check if pyliftover is available
    if not PYLIFTOVER_AVAILABLE:
        logger.error("pyliftover is not installed. Install with: pip install pyliftover")
        sys.exit(1)
    
    try:
        # Initialize liftover
        liftover_obj = GenomeLiftOver(args.from_assembly, args.to_assembly, args.chain_file)
        
        total_stats = {'total': 0, 'successful': 0, 'failed': 0}
        
        # Process each input file
        for input_file in args.input_files:
            if not os.path.exists(input_file):
                logger.warning(f"Input file not found: {input_file}")
                continue
            
            # Determine output file
            if output_mode == 'single':
                output_file = args.output
            else:
                base_name = os.path.splitext(os.path.basename(input_file))[0]
                output_file = os.path.join(args.output_dir, f"{base_name}{args.suffix}.bed")
            
            # Perform liftover
            try:
                stats = liftover_bed_file(input_file, output_file, liftover_obj, args.keep_failed)
                
                # Update total stats
                for key in total_stats:
                    total_stats[key] += stats[key]
                
                logger.info(f"File {input_file}: {stats['successful']}/{stats['total']} intervals lifted successfully")
                
            except Exception as e:
                logger.error(f"Failed to process {input_file}: {e}")
                continue
        
        # Print summary
        success_rate = (total_stats['successful'] / total_stats['total']) * 100 if total_stats['total'] > 0 else 0
        logger.info(f"\nSummary:")
        logger.info(f"Total intervals processed: {total_stats['total']}")
        logger.info(f"Successfully lifted: {total_stats['successful']} ({success_rate:.1f}%)")
        logger.info(f"Failed liftover: {total_stats['failed']}")
        
    except Exception as e:
        logger.error(f"Liftover failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()