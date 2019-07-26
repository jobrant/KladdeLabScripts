#!/usr/bin/env python

"""
This program determines which NFRs are potentially cut by BsaWI sites.
This reads the records file from rrbs.py and gives as output a file with
all cut sites and the patch lengths, positions (5', 3' or both) and coordinates,
as well as the fragment length and IDs of fragment and read.
Also prints out to the screen the total number of NFRs analyzed 
and the number of cuts at each location and the percentage of cut.


"""

import argparse
import sys

class Patch(object):
    start = 0
    end = 0
    length = 0
    direction = "f"

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.length = abs(end - start + 1)

class PatchFinder(object):
    state = "start"
    patch_start = 0

    def reset(self):
        self.state = "start"
        self.patch_start = 0

    def consume(self, ch, pos):
        if self.state == "start":
            if ch == ' ':
                pass
            elif ch == '*':
                self.state = 'BeginPatch'
                self.patch_start = pos
            elif ch == '#':
                return None
        elif self.state == 'BeginPatch':
            if ch == '+':
                self.state = "InPatch"
                pass
            else:
                return None
        elif self.state == 'InPatch':
            if ch in "*+":
                pass
            else:
                return Patch(self.patch_start, pos - 1)

        return False

    def search(self, s, direction):
        self.reset()
        if direction == 1:
            indexes = range(len(s))
        else:
            indexes = range(len(s)-1, -1, -1) 

        for i in indexes:
            a = self.consume(s[i], i)
            if a is None:
                return None
            if a.__class__.__name__ == 'Patch':
                if direction == -1:
                    s = a.end + 2
                    a.end = a.start
                    a.start = s
                return a
        return None

class Main(object):
    filename = ""
    outfile = ""
    minlength = 0
    PF = PatchFinder()
    nlines = 0
    cut5 = 0
    cut3 = 0
    cut2 = 0
    cut1 = 0


    def run(self):
        with open(self.outfile, "w") as out:
            out.write("Cut\tPatch1_Start\tPatch1_End\tPatch1_Length\tPatch2_Start\tPatch2_End\tPatch2_Length\tFragment_Length\tFragment_ID\tRead_ID\n")
            with open(self.filename, "r") as f:
                for line in f:
                    self.nlines += 1
                    line = line.rstrip("\n").split("\t")
                    item = line[0].split(':')
                    item = item[1].split('-')
                    FragLen = int(item[1]) - int(item[0])

                    map = line[3]
                    p1 = self.PF.search(map, 1)
                    p2 = self.PF.search(map, -1)
                    
                    x = (1 if p1 else 0) + (2 if p2 else 0)

                    if x == 0:
                        pass
                    elif x == 1:
                        if p1.length > int(self.minlength):
                            out.write("5\t{}\t{}\t{}\tNA\tNA\tNA\t{}\t{}\t{}\n".format(p1.start, p1.end, p1.length, FragLen, line[0], line[1]))
                            self.cut5 += 1
                    elif x == 2:
                        if p2.length > int(self.minlength):
                            out.write("3\tNA\tNA\tNA\t{}\t{}\t{}\t{}\t{}\t{}\n".format(p2.start, p2.end, p2.length, FragLen, line[0], line[1]))
                            self.cut3 += 1
                    elif x == 3:
                        if p1.length > int(self.minlength):
                            if p1.start == p2.start and p1.end == p2.end:
                                out.write("1\t{}\t{}\t{}\tNA\tNA\tNA\t{}\t{}\t{}\n".format(p1.start, p1.end, p1.length, FragLen, line[0], line[1]))
                                self.cut1 += 1
                            else:
                                out.write("2\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(p1.start, p1.end, p1.length, p2.start, p2.end, p2.length, FragLen, line[0], line[1]))
                                self.cut2 += 1
        Perc5 = 100 * (float(self.cut5)) / (float(self.nlines))
        Perc3 = 100 * (float(self.cut3)) / (float(self.nlines))
        Perc2 = 100 * (float(self.cut2)) / (float(self.nlines))
        Perc1 = 100 * (float(self.cut1)) / (float(self.nlines))

        print("Total number of NFRs analyzed = {}\nNumber cut on 5' end = {}\nPercent cut on 5' end = {:0.2f}%\nNumber cut on 3' end = {}\nPercent cut on 3' end = {:0.2f}%\nNumber with cuts on both ends = {}\nPercent with two cut NFRs = {:0.2f}%\nNumber with one NFR cut on both ends = {}\nPercent with one NFR cut on both ends = {:0.2f}%".format(self.nlines, self.cut5, Perc5, self.cut3, Perc3, self.cut2, Perc2, self.cut1, Perc1))

parser = argparse.ArgumentParser(description = "This program determines which NFRs are potentially cut by BsaWI sites. "
                                               "This reads the records file from rrbs.py and gives as output a file "
                                               "with all cut sites and the patch lengths, positions (5', 3' or both) "
                                               "and coordinates, as well as the fragment length and IDs of fragment "
                                               "and read. Also prints out to the screen the total number of NFRs "
                                               "analyzed and the number of cuts at each location and the percentage of "
                                               "cut")

parser.add_argument("filename", help = "Name of the records.tsv file to be analyzed")
parser.add_argument("outfile", help = "Name of output file")
parser.add_argument("minlength", help = "Minimum length of NFR to be counted")

if __name__ == "__main__":
    args = parser.parse_args()
    M = Main()
    M.filename = args.filename
    M.outfile = args.outfile
    M.minlength = args.minlength
    M.run()


