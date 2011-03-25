#!/usr/bin/env python
"""Convert CNV wig file to a BED file for overlapping.

This takes the CNV wiggle file from:

http://intermine.modencode.org/release-21/objectDetails.do?id=935000034

and converts it into a BED file of regions which are likely to have copy
number variations. The config variable defines the score required
for a region to have a CNV, along with the number of contiguous segments
that must have that score. Segment sizes are 1kb.
"""
import os
import sys
import csv
import itertools
import subprocess

config = dict(threshold=100.0, group_size=5)

def main(wig_file):
    simple_file = convert_to_simple(wig_file)
    write_bed(simple_file, config)

def write_bed(in_file, config):
    out_file = "%s.bed" % os.path.splitext(in_file)[0]
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect="excel-tab")
        for chrom, start, end in simple_pass(in_file, config):
            writer.writerow([chrom, start, end])

def simple_pass(in_file, config):
    info = read_simple(in_file)
    for items in itertools.izip(*[read_simple(in_file)]*config["group_size"]):
        chrs = set(x[0] for x in items)
        positions = [x[1] for x in items]
        scores = [x[2] for x in items]
        if len(chrs) == 1 and min(scores) > config["threshold"]:
            dist = (positions[-1] - positions[0]) + (positions[1] - positions[0])
            yield chrs.pop(), positions[0], positions[0] + dist

def read_simple(in_file):
    with open(in_file) as in_handle:
        for line in in_handle:
            chr, pos, score = line.split()
            yield "chr%s" % chr, int(pos), float(score)

def convert_to_simple(wig_file):
    out_file = "%s.txt" % os.path.splitext(wig_file)[0]
    if not os.path.exists(out_file):
        cl = ["wiggle_to_simple.py", wig_file]
        with open(out_file, "w") as out_handle:
            subprocess.check_call(cl, stdout=out_handle)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
