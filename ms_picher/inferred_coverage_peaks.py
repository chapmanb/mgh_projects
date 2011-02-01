#!/usr/bin/env python
"""Prepare peak calls for PIChER reads using inferred coverage.

Driver script that for each experiment:
- Does bowtie alignments
- Prepares inferred coverage BED files from paired end alignments
- Uses BED files to predict peaks
- Provides a BigWig inferred coverage display file

Usage:
  inferred_coverage_peaks.py <YAML configuration file>
"""
import os
import sys
import subprocess

import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    align_dir = os.getcwd()
    for cur_fastq in config["fastq"]:
        for ref in config["ref"]:
            fastq_one, fastq_two = _get_align_files(cur_fastq, config)
            base = "%s-%s" % (cur_fastq["name"], os.path.basename(ref))
            align_file = bowtie_to_sam(fastq_one, fastq_two, ref, base,
                                       align_dir, config)
            bam_file = sam_to_bam(align_file, fastq_one, fastq_two, ref,
                                  config_file)
            wig_file = bam_to_wig(bam_file)

def bam_to_wig(bam_file):
    cl = ["bam_to_wiggle.py", bam_file]
    subprocess.check_call(cl)

def sam_to_bam(align_file, fastq_one, fastq_two, ref, config_file):
    dir, fname = os.path.split(ref)
    sam_ref = os.path.join(dir, os.pardir, "seq", "%s.fa" % fname)
    cl = ["picard_sam_to_bam.py", config_file, align_file, sam_ref, fastq_one]
    if fastq_two:
        cl.append(fastq_two)
    subprocess.check_call(cl)
    return "%s-sort.bam" % os.path.splitext(align_file)[0]

def _get_align_files(cur_fastq, config):
    fnames = [os.path.join(config["fastq_dir"], f) for f in cur_fastq["files"]]
    assert len(fnames) in [1, 2]
    if len(fnames) == 1:
        fnames.append(None)
    return tuple(fnames)

def bowtie_to_sam(fastq_file, pair_file, ref_file, out_base, align_dir, config):
    """Before a standard or paired end alignment with bowtie.
    """
    out_file = os.path.join(align_dir, "%s.sam" % out_base)
    if not os.path.exists(out_file):
        cl = [config["program"]["bowtie"]]
        cl += ["-q", "--solexa1.3-quals",
               "-v", config["algorithm"]["max_errors"],
               "-X", 1000, # matches bwa sampe default size
               "--best",
               "--strata",
               "--sam",
               ref_file]
        cl += config["algorithm"]["align_criteria"].split()
        if pair_file:
            cl += ["-1", fastq_file, "-2", pair_file]
        else:
            cl += [fastq_file]
        cl += [out_file]
        cl = [str(i) for i in cl]
        print " ".join(cl)
        child = subprocess.check_call(cl)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
