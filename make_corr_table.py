#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:17:53 2019

@author: jasonbrant
"""

import argparse
import sys
import os

def processFiles(dir1_name, dir2_name, site):
    for filename in os.listdir(dir1_name):
        if filename.startswith(site) and filename.endswith(".freqs.csv")\
           and os.path.exists(dir2_name + "/" + filename):
            # Generate output file name by stripping out .fa extension
            outfile = os.path.splitext(filename)[0] + ".correlation.txt"
            outfile = open(outfile, "w+")
            outfile.write("Pos\trep1\trep2\n")
            list1 = get_freqs(dir1_name + "/" + filename)
            list2 = get_freqs(dir2_name + "/" + filename)
            for pair in zip(list1, list2):
                (pos1, freq1), (pos2, freq2) = pair
                
                outfile.write("{}\t{}\t{}\n".format(pos1, freq1, freq2))
            outfile.close()

def get_freqs(filename):
    freqs_list = []
    with open(filename, 'r') as rep1:
        for line in rep1:
            if line.startswith("#"):
                continue
            if line.strip() == "":
                break
            if line.startswith("Pos"):
                continue
            line = line.split()
            freqs_list.append((line[0], line[2]))
    return freqs_list

parser = argparse.ArgumentParser(description = "This script will make a table that can be used to perform correlation "
                                               "analysis with the output of methylmapper (.freqs.csv files). "
                                               "Typically used for finding correlation of NGC replicates. "
                                               "When calling this script, provide the path to each replicates directory "
                                               "and the site you are interested in (e.g. HCG). "
                                               "Returns a three column table with site position and site frequency for "
                                               "each replicate.")
parser.add_argument("dir1_name", help = "Directory containing the first replicate")
parser.add_argument("dir2_name", help = "Directory containing the second replicate")
parser.add_argument("site", help = "Specifies what frequency table site to read. Should be either HCG or GCH.")
        
if __name__ == "__main__":
    args = parser.parse_args()
    processFiles(args.dir1_name, args.dir2_name, args.site)
