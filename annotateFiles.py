#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon July 31 2019

@author: jobrant
"""

#from __future__ import print_function
#import subprocess
import sys
import argparse
import os
import pandas as pd
from subprocess import run

# # check python version
# if sys.version_info[0] == 3:
#     print("Python version 3 detected")
#     run = subprocess.run
# else:
#     print("Python version 2 detected")
#     run = subprocess.call

def processFiles(dir_name, file_prefix):
    outfiles = []
    for filename in os.listdir(dir_name):
        if filename.startswith(file_prefix) and filename.endswith("xlsx"):
            xlsx = pd.ExcelFile(filename)
            sheet_names = xlsx.sheet_names
            temp1 = xlsx.parse(sheet_names[0])
            temp2 = xlsx.parse(sheet_names[1])
            temp1['Strand'] = 0
            temp2['Strand'] = 0
            #temp1['Unique_ID'] = temp1[['Chrom', 'Start']].apply(lambda x: ".".join(x), axis = 1)
            #temp2['Unique_ID'] = temp2[['Chrom', 'Start']].apply(lambda x: ".".join(x), axis = 1)
            temp1['Unique_ID'] = temp1['Chrom'] + '.' + temp1['Start'].astype(str)
            temp2['Unique_ID'] = temp2['Chrom'] + '.' + temp2['Start'].astype(str)
            header = ["Unique_ID", "Chrom", "Start", "End", "Strand", "Gene", "ID", "Accession", "Class", "ExNum",
                      "log2(FC)"]
            outfile1 = os.path.splitext(filename)[0] + sheet_names[0] + ".tsv"
            outfile2 = os.path.splitext(filename)[0] + sheet_names[1] + ".tsv"
            outfiles.append(outfile1)
            outfiles.append(outfile2)
            outfile1 = open(outfile1, "w+")
            outfile2 = open(outfile2, "w+")
            temp1.to_csv(outfile1, columns = header, sep = '\t', index = False)
            temp2.to_csv(outfile2, columns = header, sep = '\t', index = False)
            outfile1.close()
            outfile2.close()
    return(outfiles)

HOMER_COMMAND="annotatePeaks.pl {file} /apps/homer/4.10/share/homer/4.10-0/.//data/genomes/hg38/ -organism human -annStats {file2}"

def submitHomer(outfiles):
    my_env = os.environ.copy()
    annoFiles = []
    for f in outfiles:
        fout = os.path.splitext(f)[0] + ".annotated.tsv"
        fout_base_name = '.'.join(fout.split('.')[0:3])
        annoFiles.append(fout_base_name)
        fout = open(fout, "w+")
        ann_stats = os.path.splitext(f)[0] + ".annStats.txt"
        ann_stats_out = open(ann_stats, "w+")
        run(HOMER_COMMAND.format(file = f, file2 = ann_stats), shell = True, stdout = fout, stderr = ann_stats_out, env = my_env)
        fout.close()
        ann_stats_out.close()
    annoFiles = set(annoFiles)
    return(annoFiles)

def mergeFiles(annoFiles):
    for f in annoFiles:
        baseNames = f.split('.')
        annotated1 = pd.read_csv(str(f) + ".annotGenes" + "-" + str(baseNames[0]) + ".annotated.tsv", delimiter='\t', encoding='utf-8')
        annotated2 = pd.read_csv(str(f) + ".annotGenes" + "-" + str(baseNames[2]) + ".annotated.tsv", delimiter='\t', encoding='utf-8')
        annotated1.rename(columns={annotated1.columns[0]: "Unique_ID"}, inplace=True)
        annotated2.rename(columns={annotated2.columns[0]: "Unique_ID"}, inplace=True)
        original1 = pd.read_csv(str(f) + ".annotGenes" + "-" + str(baseNames[0]) + ".tsv", delimiter = '\t', encoding = 'utf-8')
        original2 = pd.read_csv(str(f) + ".annotGenes" + "-" + str(baseNames[2]) + ".tsv", delimiter='\t', encoding='utf-8')
        combined1 = pd.merge(original1, annotated1, on = 'Unique_ID')
        combined2 = pd.merge(original2, annotated2, on = 'Unique_ID')
        print(annotated1.columns)
        print(annotated2.columns)
        print(original1.columns)
        print(original2.columns)
        print(combined1.columns)
        print(combined2.columns)
        fout = str(f) + ".annotGenes" + "-" + str(baseNames[0]) + ".annotated.xlsx"
        writer = pd.ExcelWriter(fout, engine = 'xlsxwriter')
        fout = open(fout, "w+")
        combined1.to_excel(writer, sheet_name = "Genes-" + str(baseNames[0]), index=False)
        combined2.to_excel(writer, sheet_name="Genes-" + str(baseNames[2]), index=False)
        fout.close()

parser = argparse.ArgumentParser(description = "Annotation testing")
parser.add_argument("dir_name", help = "Directory containing the files")
parser.add_argument("file_prefix", help = "Filename prefix")

if __name__ == "__main__":
    args = parser.parse_args()
    outfiles = processFiles(args.dir_name, args.file_prefix)
    annoFiles = submitHomer(outfiles)
    mergeFiles(annoFiles)

