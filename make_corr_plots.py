#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:17:53 2019

@author: jasonbrant
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib
import matplotlib.pyplot as plt

def processFiles(dir1_name, dir2_name, site):
    total_array = []
    for filename in os.listdir(dir1_name):
        array = []
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
                array.append((freq1,freq2))
                total_array.append((freq1,freq2))
            df = pd.DataFrame(array)
            df = df.astype(float)
            corr1 = df[0].corr(df[1])
            print("{}\t{}\n".format(filename, corr1))
            outfile.close()
            b, m = polyfit(df[0], df[1], 1)
            plt.scatter(df[0], df[1])
            plt.plot(df[0], b + m * df[0], '-')
            plt.suptitle(os.path.splitext(filename)[0].split('.')[0], y = 0.97, fontsize = 12)
            subtitle = "$r^2$ = {:.3f}"
            plt.title(subtitle.format(corr1), fontsize = 10)
            plt.xlabel('rep1')
            plt.ylabel('rep2')
            plt.savefig(os.path.splitext(filename)[0] + ".png")
            plt.close()

    total_df = pd.DataFrame(total_array)
    total_df = total_df.astype(float)
    corr2 = total_df[0].corr(total_df[1])
    print("Total correlation number is:\t{}".format(corr2))
    b, m = polyfit(total_df[0], total_df[1], 1)
    plt.scatter(total_df[0], total_df[1])
    plt.plot(total_df[0], b + m * total_df[0], '-')
    plt.suptitle(dir1_name.split('/')[0], y = 0.97, fontsize = 12)
    total_suptitle = "$r^2$ = {:.3f}"
    plt.title(total_suptitle.format(corr2), fontsize = 10)
    plt.savefig(dir1_name.split('/')[0] + ".png")
    plt.close()


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

parser = argparse.ArgumentParser(description = "This script will generate correlation plots for each gene/fragment "
                                               "between replicates using the output of methylmapper (.freqs.csv files). "
                                               "Typically used for finding correlation of NGC replicates. "
                                               "When calling this script, provide the path to each replicate's directory "
                                               "and the site you are interested in (e.g. HCG). "
                                               "This will generate a scatter plot for each gene/fragment with a "
                                               "regression line and r-squared and a three column table with site "
                                               "position and site frequency. It also generates a plot and table for the "
                                               "combined data for every gene/fragment in the directories. ")

parser.add_argument("dir1_name", help = "Directory containing the first replicate")
parser.add_argument("dir2_name", help = "Directory containing the second replicate")
parser.add_argument("site", help = "Specifies what frequency table site to read. Should be either HCG or GCH.")
        
if __name__ == "__main__":
    args = parser.parse_args()
    processFiles(args.dir1_name, args.dir2_name, args.site)

