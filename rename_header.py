#!/usr/bin/env python
# 2019 Brant and Riva
"""
Renames fasta header from fasta files generated by bsDraw.
This is neccessary when using these files with rrbs.py.

"""

import os
import os.path
import sys
import argparse

def readNames(filename):
    dict = {}
    with open(filename, "r") as f:
        for line in f:
            #print(line)
            tok = line.split('\t')
            dict[tok[3]] = tok[5]
    return dict

def processFiles(dir, dict):
    for filename in os.listdir(dir):
        if filename.endswith(".fa"):
            # Generate output file name by stripping out .fa extension
            outfile = os.path.splitext(filename)[0] + ".renamed.fa"
            processFile(filename, outfile, d)

def processFile(filename, outfile, dict):
    # Open output file
    with open(outfile, "w") as out:
        # Open original file
        with open(filename, "r") as f:
            # Read header line
            hdr = f.readline().rstrip("\r\n")
            # Extract current name
            frag = hdr[1:]
            # Replace it with new one
            #out.write(">" + dict[frag] + "\n")
            out.write(">" + dict[frag])
            # Then simply copy all other lines
            for line in f:
                out.write(line)

parser = argparse.ArgumentParser(description = "Run rename_header.py -h for help. "
                                               "Renames fasta header from fasta files generated by bsDraw. "
                                               "This is necessary when using these files with rrbs.py. "
                                               "Provide a file, in bed-like format with the coordinates of fragment, "
                                               "and the old and new header names. This usually made by modifying the "
                                               "report file from fastools digest. (tsv format with chr, start, end, "
                                               "fragment_ID, length, new_ID ")

parser.add_argument("names_file", help = "bed-like file containing the reference fragment coords plus old and new names. ")
parser.add_argument("dir_name", help = "directory to read files from. ")

if __name__ == "__main__":
    args = parser.parse_args()
    d = readNames(args.names_file)
    processFiles(args.dir_name, d)
