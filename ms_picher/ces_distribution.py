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
        out_file = "%s-distribution.pdf" % rconfig["name"]
        plot_count_distribution(out_file, counts + rcounts)

def plot_count_distribution(out_file, all_counts):
    tbl = {'reads.per.million': rpy.FloatVector([x[0] for x in all_counts]),
           'exp': rpy.StrVector([x[1] for x in all_counts]),
           'region': rpy.StrVector([x[2] for x in all_counts])}
    rpy.r.assign("exp.data", rpy.r['data.frame'](**tbl))
    rpy.r.assign("out.file", out_file)
    rpy.r('''
    library(ggplot2)
    p <- ggplot(exp.data, aes(x=reads.per.million, colour=region)) +
         stat_bin(aes(y = ..ncount..), binwidth=0.1,
                  geom="line", position="identity") +
         scale_x_log10() + facet_wrap(~ exp, ncol=1) +
         xlab("Reads (per million)") +
         ylab("Counts (normalized)")
    ggsave(out.file, p, width=8, height=11)
    ''')

def region_counts(bams, regions, region_name):
    counts = []
    for bam in bams:
        for space, start, end in regions:
            counts.append([bam.read_count(space, start, end),
                           bam.name, region_name])
    return counts

class NormalizedBam:
    """Prepare and query an alignment BAM file for normalized read counts.
    """
    def __init__(self, name, fname, picard):
        self.name = name
        self._bam = pysam.Samfile(fname, "rb")
        picard.run_fn("picard_index", fname)
        #self._total = 1e6
        self._total = sum(1 for _ in self._bam.fetch())
        print name, self._total

    def read_count(self, space, start, end):
        """Retrieve the normalized read count in the provided region.
        """
        read_counts = 0
        for read in self._bam.fetch(space, start, end):
            read_counts += 1
        return self._normalize(read_counts, self._total)

    def _normalize(self, count, total):
        """Normalize to reads per million.
        """
        return float(count) / float(total) * 1e6

def random_regions(base, n, size):
    """Generate n random sized regions in the provided base spread.
    """
    base_info = collections.defaultdict(list)
    for space, start, end in base:
        base_info[space].append(start)
        base_info[space].append(end)
    regions = []
    print "Generating", n, "regions of size", size
    for _ in range(n):
        space = random.choice(base_info.keys())
        pos = random.randint(min(base_info[space]), max(base_info[space]))
        regions.append([space, pos-size, pos+size])
    return regions

def read_ces_regions(in_handle, size):
    """Read chrX target regions defined in an input file.
    """
    space = "chrX"
    regions = []
    for line in in_handle:
        parts = line.split()
        try:
            peak = int(parts[-1])
        except (ValueError, IndexError):
            continue
        regions.append([space, peak-size, peak+size])
    return regions

if __name__ == "__main__":
    main(*sys.argv[1:])
