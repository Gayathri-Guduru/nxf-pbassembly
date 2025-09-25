#!/usr/bin/env python3
# extract_16S.py - Extract 16S ribosomal RNA sequences from Prokka output

import pandas as pd
from Bio import SeqIO
import os
import sys
import argparse

def extract_16s_sequences(tsv_file, ffn_file, output_dir, sample_id):
    """Extract 16S ribosomal RNA sequences from Prokka output."""
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Output file paths
    csv_output = os.path.join(output_dir, f"{sample_id}_16S.csv")
    fasta_output = os.path.join(output_dir, f"{sample_id}_16S.fasta")
    
    try:
        # Read the TSV file
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Filter for 16S rRNA entries
        filtered_df = df[(df['ftype'] == 'rRNA') & (df['product'].str.contains('16S ribosomal RNA'))]
        
        # Save filtered entries to CSV
        filtered_df.to_csv(csv_output, index=False)
        
        if filtered_df.empty:
            # Create an empty fasta file
            with open(fasta_output, 'w') as f:
                pass
        else:
            # Get the locus tags
            locus_tags = filtered_df['locus_tag'].tolist()
            
            # Read the FFN file and extract sequences for the locus tags
            sequences = []
            for record in SeqIO.parse(ffn_file, "fasta"):
                if record.id in locus_tags:
                    sequences.append(record)
            
            # Write the extracted sequences to a new FASTA file
            if sequences:
                SeqIO.write(sequences, fasta_output, "fasta")
            else:
                # Create an empty fasta file
                with open(fasta_output, 'w') as f:
                    pass
    except Exception as e:
        # Create empty output files in case of error
        with open(csv_output, 'w') as f:
            f.write(f"Error: {str(e)}\n")
        with open(fasta_output, 'w') as f:
            pass
        return 1  # Error exit code
    
    return 0  # Success exit code

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract 16S rRNA sequences from Prokka output.')
    parser.add_argument('tsv_file', help='Path to the TSV annotation file')
    parser.add_argument('ffn_file', help='Path to the FFN nucleotide sequences file')
    parser.add_argument('output_dir', help='Directory to save output files')
    parser.add_argument('sample_id', help='Sample identifier')
    
    args = parser.parse_args()
    
    sys.exit(extract_16s_sequences(args.tsv_file, args.ffn_file, args.output_dir, args.sample_id))
