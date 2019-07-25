#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon Jun 17 14:46:26 2019

@author: jasonbrant
Run filter_fasta.py -h for help. 
Use this to filter a fasta alignment file for use with methylmapper and Rhonda's R clustering scripts. 
This will read an alignment fasta, determine the reference length, and the max length of the alignments. 
This will then write out a new alignment fasta of those aligned reads that meet user defined length paramters.
Command line argument 1 is the fasta file, argument 2 is percentage of max aligned read length. 
So 90 would mean each read must be at least 90% the length of the max aligned read length. 
By default, the output fasta file will be named like the input file, but ending in filtered.fa. 
Optionally, you can specify an output file name using the -o option.

"""

from Bio import SeqIO
import argparse


# reads in fasta file and returns a tuple containing the read ID and
# the sequence
def read_fasta(fname):
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

def filter_seqs(filename, perc_len):
    ref_len = 0
    alignment_lengths = []
    to_keep = set()
    
    for name, seq in read_fasta(filename):
        if name.startswith(">"):
            if not name.startswith(">m"):
                ref_len = len(seq)
        if name.startswith(">m"):
            alignment_lengths.append(len(seq))
    for name, seq in read_fasta(filename):
        if len(seq) >= (max(alignment_lengths) * perc_len/100):
            to_keep.add(name.strip('>'))

    print("\nReference Length\t{}\n\nMax alignment length\t{}\nMin alignment length\t{}\n".format(ref_len, max(alignment_lengths), min(alignment_lengths)))
    return to_keep

def output_set(filename, nameset, output_name):
    seqiter = SeqIO.parse(open(filename), 'fasta')
    with open(output_name, "w+") as handle:
        SeqIO.write((seq for seq in seqiter if seq.id in nameset), handle, 'fasta')

parser = argparse.ArgumentParser(description = "Filter an alignment fasta (produced by bsDraw) by percentage of max aligned read length")
parser.add_argument("filename", help = "Input fasta file to be filtered")
parser.add_argument("percent", help = "Percent", type = float)
parser.add_argument("-o", "--output", help = "Specify an optional output file name")
        
if __name__ == "__main__":
    args = parser.parse_args()
    if args.output is None:
        output_name = args.filename[:-3] + ".filtered.fa"
    else:
        output_name = args.output
    names = filter_seqs(args.filename, args.percent)
    output_set(args.filename, names, output_name)
    


