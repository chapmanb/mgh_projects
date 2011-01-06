#!/usr/bin/env python
"""Examine distribution of reads in Chromosome Entry Sites (CES) versus genome.

Plot the read distribution of BAM files in chromosome entry sites versus
a random background on chromosome X.

Usage:
  ces_distribution.py <config file>
"""
import os
import sys

import yaml

from bcbio import picard

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard_run = picard.PicardRunner(config["programs"]["picard"])
    for align in config["alignments"]:
        picard.run.picard_index(picard_run, align["file"])

if __name__ == "__main__":
    main(*sys.argv[1:])
