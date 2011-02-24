#!/usr/bin/env python
"""Convert BAM input file into fasta formatted output.

Usage:
    bam_align_to_fasta.py <Picard dir> <BAM dir>

Requires:
    Picard (http://picard.sourceforge.net/)
    fastx toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
"""
import glob
import os
import sys
import subprocess

from Mgh.Picard import PicardRunner

def main(picard_dir, bam_dir):
    picard = PicardRunner(picard_dir)
    for fname in glob.glob(os.path.join(bam_dir, "*-sort-filter.bam")):
        fastq_out = bam_to_fastq(fname, picard)
        fastq_to_fasta(fastq_out)

def bam_to_fastq(in_file, picard):
    out_file = "%s.fastq" % os.path.splitext(os.path.basename(in_file))[0]
    if not os.path.exists(out_file):
        opts = [("INPUT", in_file), ("FASTQ", out_file)]
        picard.run("SamToFastq", opts)
    return out_file

def fastq_to_fasta(fastq_file):
    out_file = "%s.fa" % os.path.splitext(fastq_file)[0]
    cl = ["fastq_to_fasta", "-n", "-i", fastq_file, "-o", out_file]
    subprocess.check_call(cl)

if __name__ == "__main__":
    main(*sys.argv[1:])
