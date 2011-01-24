#!/usr/bin/env python
"""Create scatter plot comparing BAM files for reproducibility.

Randomly chooses chromosomal regions, plotting scatter of read counts in
each file.

Usage:
  compare_bam_files.py <config file>
"""
import sys
import os

import yaml
import rpy2.robjects as rpy

from bcbio.bam import counts
from bcbio.picard import PicardRunner

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    region_size = int(config["algorithm"]["region_size"])
    num_regions = int(config["algorithm"]["num_regions"])
    bam_files = prepare_bam_files(config)
    regions = counts.random_regions(bam_files[0].all_regions(),
                                    num_regions, region_size)
    values = []
    for i, bam1 in enumerate(bam_files):
        for bam2 in bam_files[i+1:]:
            values.extend(compare_bam_files(bam1, bam2, regions))
    scatter_plot(values, bam1.name, bam2.name, region_size,
                 config["files"]["coverage_pdf"])

def scatter_plot(counts, name1, name2, region_size, out_file):
    tbl = {'name': rpy.StrVector([x[0] for x in counts]),
           'a': rpy.FloatVector([x[1] for x in counts]),
           'b': rpy.FloatVector([x[2] for x in counts])}
    rpy.r.assign("exp.data", rpy.r['data.frame'](**tbl))
    rpy.r.assign("out.file", out_file)
    rpy.r.assign("x.label", name1)
    rpy.r.assign("y.label", name2)
    rpy.r.assign("rsize", region_size)
    rpy.r('''
    library(ggplot2)
    p <- ggplot(exp.data, aes(a, b)) + geom_point(position="jitter") +
         scale_x_log10() + scale_y_log10() +
         facet_wrap(~ name, ncol=1) +
         xlab('') + ylab('') +
         opts(title = paste("Coverage of ", rsize, "bp regions"))
    ggsave(out.file, p, width=8, height=11)
    ''')

def compare_bam_files(bam1, bam2, regions):
    counts1 = []
    counts2 = []
    for space, start, end in regions:
        counts1.append(bam1.read_count(space, start, end))
        counts2.append(bam2.read_count(space, start, end))
    name = "%s v %s" % (bam1.name, bam2.name)
    return [(name, a, b) for (a, b) in zip(counts1, counts2) if a > 0 or b > 0]

def prepare_bam_files(config):
    picard = PicardRunner(config["programs"]["picard"])
    ready_bams = []
    for bconfig in config["bams"]:
        ready_bams.append(
            counts.NormalizedBam(bconfig["name"], bconfig["file"], picard, False))
    return ready_bams

if __name__ == "__main__":
    main(*sys.argv[1:])
