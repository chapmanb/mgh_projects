#!/usr/bin/env python
"""Examine distribution of reads in Chromosome Entry Sites (CES) versus genome.

Plot the read distribution of BAM files in chromosome entry sites versus
a random background on chromosome X.

Usage:
  ces_distribution.py <config file>
"""
import os
import sys
import random
import collections

import yaml
import pysam
import rpy2.robjects as rpy

from bcbio.picard import PicardRunner
from bcbio.bam.counts import NormalizedBam, random_regions

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = PicardRunner(config["programs"]["picard"])
    bams = [NormalizedBam(align["name"], align["file"], picard)
            for align in config["alignments"]]

    for rconfig in config["regions"]:
        with open(rconfig["file"]) as in_handle:
            assert rconfig["name"] == "CES"
            regions = read_ces_regions(in_handle, rconfig["size"])
        rregions = random_regions(regions, rconfig["random"], rconfig["size"])
        counts = region_counts(bams, regions, rconfig["name"])
        rcounts = region_counts(bams, rregions, "background")
        out_file = "%s-distribution%s.pdf" % (rconfig["name"],
                                              "-smoothed" if rconfig["smoothed"] else "")
        plot_count_distribution(out_file, counts + rcounts, rconfig["smoothed"])

def plot_count_distribution(out_file, all_counts, smoothed):
    tbl = {'reads.per.million': rpy.FloatVector([x[0] for x in all_counts]),
           'exp': rpy.StrVector([x[1] for x in all_counts]),
           'region': rpy.StrVector([x[2] for x in all_counts])}
    rpy.r.assign("exp.data", rpy.r['data.frame'](**tbl))
    rpy.r.assign("out.file", out_file)
    rpy.r.assign("smoothed", smoothed)
    rpy.r('''
    library(ggplot2)
    p <- ggplot(exp.data, aes(x=reads.per.million, colour=region))
    p <- p + if (smoothed) geom_density(aes(y = ..scaled..), fill=NA) else
                           stat_bin(aes(y = ..ncount..), binwidth=0.1,
                                    geom="line", position="identity")
    p <- p +
         scale_x_log10(breaks=c(0.1, 1, 10, 100)) +
         facet_wrap(~ exp, ncol=1) +
         xlab("Reads per million") +
         ylab("Normalized counts") +
         opts(panel.background = theme_blank(),
              panel.grid.minor = theme_blank(),
              panel.grid.major = theme_blank(),
              axis.ticks = theme_blank(),
              axis.text.y = theme_blank())
    ggsave(out.file, p, width=8, height=11)
    ''')

def region_counts(bams, regions, region_name):
    counts = []
    for bam in bams:
        for space, start, end in regions:
            counts.append([bam.read_count(space, start, end),
                           bam.name, region_name])
    return counts

def read_ces_regions(in_handle, size):
    """Read chrX target regions defined in an input file.
    """
    spread = size // 2
    space = "chrX"
    regions = []
    for line in in_handle:
        parts = line.split()
        try:
            peak = int(parts[-1])
        except (ValueError, IndexError):
            continue
        regions.append([space, peak-spread, peak+spread])
    return regions

if __name__ == "__main__":
    main(*sys.argv[1:])
