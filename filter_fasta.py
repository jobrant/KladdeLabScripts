#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
FASTA Alignment Filter

This script filters FASTA alignment files for use with methylmapper and R clustering scripts.
It filters aligned reads based on a minimum length requirement relative to the reference sequence.

Usage:
    filter_fasta.py <input_fasta> <percent_length> [-o OUTPUT]

Arguments:
    input_fasta     Input FASTA alignment file (output from bsDraw.py)
    percent_length  Minimum percent length relative to reference sequence (e.g., 90)
    -o, --output    Optional output filename (default: input_name.filtered.fa)

Example:
    filter_fasta.py input.fa 90 -o filtered_output.fa

Created on: Mon Jun 17 14:46:26 2019
Author: jasonbrant
'''

from Bio import SeqIO
import argparse

def read_fasta(fname: str) -> tuple[str, str]:
    '''
    Reads in FASTA file and yields tuples of sequence IDs and sequences.
    
    Args:
        fname (str): Path to input FASTA file
        
    Yields:
        tuple[str, str]: Tuple containing the read ID and sequence
        
    Raises:
        FileNotFoundError: If input file doesn't exist
    '''

    with open(fname) as fp:
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            line = line.strip('-')
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

def filter_seqs(filename: str, perc_len: float) -> set[str]:
    '''
    This function determines the length of the reference sequence (first entry in FASTA) 
    and filters out alignments that don't meet user supplied parameters based on the 
    reference sequence length.
    
    Args:
        filename (str): Path to the FASTA file
        perc_len (float): Minimum percentage of reference length to keep

    Returns:
        set[str]: Set of sequence IDs that meet the length criteria
    '''

    # Get reference length from first entry
    first_entry = next(read_fasta(filename))  # Get first entry only
    ref_len = len(first_entry[1])
    
    # Calculate minimum required length
    min_required_length = ref_len * perc_len/100
    
    alignment_lengths = []
    to_keep = set()

    # Process remaining sequences
    for name, seq in read_fasta(filename):
        if name.startswith(">m"):  # Only process mapped sequences
            seq_len = len(seq)
            alignment_lengths.append(seq_len)
            if seq_len >= min_required_length:
                to_keep.add(name.strip('>'))

    print("\nReference Length\t{}\n\nMax alignment length\t{}\nMin alignment length\t{}\n"
          "Minimum required length ({}% of ref)\t{}\n"
          .format(ref_len, max(alignment_lengths), min(alignment_lengths), 
                 perc_len, min_required_length))
    
    print("Kept {} out of {} sequences".format(len(to_keep), len(alignment_lengths)))
    return to_keep

def output_set(filename: str, nameset: set[str], output_name: str) -> None:
    '''
    Writes out the reference sequence and filtered alignments to a new FASTA file.
    
    Args:
        filename (str): Path to input FASTA file
        nameset (set[str]): Set of sequence IDs to keep
        output_name (str): Path for output FASTA file
        
    Returns:
        None
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        PermissionError: If unable to write to output location
    '''

    seqiter = SeqIO.parse(open(filename), 'fasta')
    with open(output_name, "w+") as handle:
        SeqIO.write((seq for seq in seqiter if seq.id in nameset), handle, 'fasta')

parser = argparse.ArgumentParser(
    description="Filter FASTA alignments based on minimum length requirement relative to reference sequence."
)

parser.add_argument(
    "filename", 
    help="Input FASTA alignment file (output from bsDraw.py)"
)
parser.add_argument(
    "percent", 
    type=float,
    help="Minimum percent length relative to reference sequence (e.g., 90)"
)
parser.add_argument(
    "-o", "--output",
    help="Optional output filename (default: input_name.filtered.fa)"
)
        
if __name__ == "__main__":
    args = parser.parse_args()
    if args.output is None:
        output_name = args.filename[:-3] + ".filtered.fa"
    else:
        output_name = args.output
    names = filter_seqs(args.filename, args.percent)
    output_set(args.filename, names, output_name)
    


