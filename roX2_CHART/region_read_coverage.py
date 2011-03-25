#!/usr/bin/env python
"""Plot read coverage over defined regions in multiple BAM files.

Usage:
  region_read_coverage.py <config file>
"""
import os
import sys

import yaml
import rpy2.robjects as rpy

from bcbio.picard import PicardRunner
from bcbio.bam.counts import NormalizedBam

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = PicardRunner(config["programs"]["picard"])
    bams = [NormalizedBam(align["name"], align["file"], picard)
            for align in config["alignments"]]
    for pconfig in config["detail_plots"]:
        out_file = "%s-coverage.pdf" % pconfig["name"]
        coverage = region_coverage(bams, pconfig["space"], pconfig["start"],
                                   pconfig["end"])
        plot_coverage(coverage, pconfig["name"], out_file)

def plot_coverage(coverage, title, out_file):
    tbl = {'position': rpy.IntVector([x[0] for x in coverage]),
           'coverage': rpy.FloatVector([x[1] for x in coverage]),
           'exp': rpy.StrVector([x[2] for x in coverage])}
    rpy.r.assign("exp.data", rpy.r['data.frame'](**tbl))
    rpy.r.assign("out.file", out_file)
    rpy.r.assign("title", title)
    rpy.r('''
    library(ggplot2)
    p <- ggplot(exp.data, aes(x=position, y=coverage)) +
         geom_area() +
         facet_wrap(~ exp, ncol=1) +
         xlab(paste("chrX:", title, sep=" ")) +
         ylab("Coverage (reads per million)") +
         opts(panel.background = theme_blank(),
              panel.grid.minor = theme_blank(),
              panel.grid.major = theme_blank(),
              axis.ticks = theme_blank())
    ggsave(out.file, p, width=8, height=11)
    ''')

def region_coverage(bams, space, start, end):
    coverage = []
    for bam in bams:
        for pos, cov in bam.coverage_pileup(space, start, end):
            coverage.append((pos, cov, bam.name))
    return coverage

if __name__ == "__main__":
    main(*sys.argv[1:])
