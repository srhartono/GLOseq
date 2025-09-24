#!/usr/bin/env python3
"""
Generate sample RFD data for testing the RFD analyzer
"""

import numpy as np
import pandas as pd

def generate_sample_rfd_data(output_file: str = "sample_rfd.bedgraph", 
                            chromosome: str = "chr1",
                            start_pos: int = 1000000,
                            num_bins: int = 1000,
                            bin_size: int = 10000):
    """
    Generate synthetic RFD data that mimics real OK-seq data patterns.
    
    Args:
        output_file: Output bedgraph filename
        chromosome: Chromosome name
        start_pos: Starting genomic position
        num_bins: Number of 10kb bins to generate
        bin_size: Size of each bin (default 10kb)
    """
    
    # Create position arrays
    positions = np.arange(start_pos, start_pos + num_bins * bin_size, bin_size)
    
    # Generate synthetic RFD signal with realistic patterns
    # Base sinusoidal pattern with multiple frequencies
    x = np.linspace(0, 8 * np.pi, num_bins)
    
    # Primary replication pattern (low frequency)
    primary_signal = 0.6 * np.sin(x)
    
    # Secondary oscillations (higher frequency)
    secondary_signal = 0.3 * np.sin(3 * x) 
    
    # Tertiary fine structure
    tertiary_signal = 0.15 * np.sin(8 * x)
    
    # Add some noise
    noise = np.random.normal(0, 0.1, num_bins)
    
    # Combine signals
    rfd_values = primary_signal + secondary_signal + tertiary_signal + noise
    
    # Ensure values are in realistic range (-1 to +1)
    rfd_values = np.clip(rfd_values, -1.0, 1.0)
    
    # Create DataFrame
    data = pd.DataFrame({
        'chr': chromosome,
        'start': positions,
        'end': positions + bin_size,
        'RFD': rfd_values
    })
    
    # Save to bedgraph format
    data.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Generated {len(data)} bins of sample RFD data")
    print(f"Genomic range: {chromosome}:{start_pos:,}-{positions[-1] + bin_size:,}")
    print(f"RFD range: {rfd_values.min():.3f} to {rfd_values.max():.3f}")
    print(f"Saved to: {output_file}")
    
    return data

def generate_multi_chromosome_data(output_file: str = "multi_chr_rfd.bedgraph"):
    """Generate sample data for multiple chromosomes"""
    
    all_data = []
    
    # Generate data for chromosomes 1, 2, and 3
    for chr_num in [1, 2, 3]:
        chr_name = f"chr{chr_num}"
        
        # Different parameters for each chromosome
        if chr_num == 1:
            start_pos, num_bins = 1000000, 800
        elif chr_num == 2:
            start_pos, num_bins = 2000000, 600
        else:  # chr3
            start_pos, num_bins = 500000, 400
        
        # Generate data for this chromosome
        chr_data = generate_sample_rfd_data(
            output_file=f"temp_{chr_name}.bedgraph",
            chromosome=chr_name,
            start_pos=start_pos,
            num_bins=num_bins
        )
        
        all_data.append(chr_data)
    
    # Combine all chromosomes
    combined_data = pd.concat(all_data, ignore_index=True)
    combined_data = combined_data.sort_values(['chr', 'start'])
    
    # Save combined file
    combined_data.to_csv(output_file, sep='\t', index=False, header=False)
    
    # Clean up temporary files
    import os
    for chr_num in [1, 2, 3]:
        temp_file = f"temp_chr{chr_num}.bedgraph"
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    print(f"\nCombined multi-chromosome data saved to: {output_file}")
    print(f"Total bins: {len(combined_data)}")
    
    return combined_data

if __name__ == "__main__":
    # Generate single chromosome sample
    generate_sample_rfd_data("sample_rfd_chr1.bedgraph")
    
    # Generate multi-chromosome sample
    generate_multi_chromosome_data("sample_rfd_multi.bedgraph")