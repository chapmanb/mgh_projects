#!/usr/bin/env python
"""Initial descriptive analysis of Millie's flowcells.
"""
import sys
import os
import subprocess
import glob

import yaml
import rpy2.robjects as robjects

from bcbio.picard import utils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    peak_beds = []
    for exp in config["exps"]:
        peak_file = macs_peaks(exp["align"], config)
        peak_beds.append(peak_file)
    unique_peaks(peak_beds, config)
    analyze_regions(config)

def analyze_regions(config):
    peak_files = sorted(glob.glob(os.path.join(config["dirs"]["peaks"],
        "*_peaks.bed")))
    feature_files = []
    for peak_file in peak_files:
        print os.path.basename(peak_file)
        feature_file = add_peak_features(peak_file)
        feature_files.append(feature_file)
    feature_summary(feature_files)

def feature_summary(csv_files):
    for csv_file in csv_files:
        print csv_file
        robjects.r.assign("feat.file", csv_file)
        robjects.r('''
          feature.tbl <- read.csv(feat.file)
          feature.sum <- summary(feature.tbl$ftype)
          feature.sum['total'] = sum(feature.sum)
          print(feature.sum)
        ''')

def add_peak_features(in_file):
    out_file = "%s-features.csv" % os.path.splitext(in_file)[0]
    if not os.path.exists(out_file):
        cl = ["Rscript", "peaks_to_features.R", in_file, out_file]
        subprocess.check_call(cl)
    return out_file

def unique_peaks(peak_beds, config):
    for i, one_file in enumerate(peak_beds):
        for two_file in peak_beds[i+1:]:
            out_file = _get_compare_out(one_file, two_file, "not")
            _intersect_files(one_file, two_file, out_file, config)
            out_file = _get_compare_out(two_file, one_file, "not")
            _intersect_files(two_file, one_file, out_file, config)
            out_file = _get_compare_out(one_file, two_file, "and")
            _intersect_files(one_file, two_file, out_file, config, True)
            out_file = _get_compare_out(two_file, one_file, "and")
            _intersect_files(two_file, one_file, out_file, config, True)

def _get_compare_out(one_file, two_file, combine_str):
    one_unique, one_shared = os.path.basename(one_file).split("_", 1)
    two_unique, two_shared = os.path.basename(two_file).split("_", 1)
    assert one_shared == two_shared
    out_file = os.path.join(os.path.dirname(one_file),
        "%s_%s_%s_%s" % (one_unique, combine_str, two_unique, one_shared))
    return out_file

def _intersect_files(one_file, two_file, out_file, config, union=False):
    if not os.path.exists(out_file):
        cl = [config["program"]["bed_intersect"], "-m", "10"]
        if not union:
            cl.append("-v")
        cl += [one_file, two_file]
        with open(out_file, "w") as out_handle:
            proc = subprocess.Popen(cl, stdout=out_handle)
            proc.wait()

def macs_peaks(in_align, config):
    bed_file = _convert_to_bed(in_align)
    peak_dir = config["dirs"]["peaks"]
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    sample_name = os.path.join(peak_dir,
            os.path.splitext(os.path.basename(in_align))[0])
    peak_file = "%s_peaks.bed" % sample_name
    if not os.path.exists(peak_file):
        cl = ["macs14", "-t", bed_file, "--name=%s" % sample_name,
              "--format=BED"]
        subprocess.check_call(cl)
    return peak_file

def _convert_to_bed(bam_file):
    bed_file = "%s.bed" % os.path.splitext(bam_file)[0]
    if not os.path.exists(bed_file):
        cl = ["bamToBed", "-i", bam_file]
        with open(bed_file, "w") as out_handle:
            proc = subprocess.Popen(cl, stdout=out_handle)
            proc.wait()
    return bed_file

if __name__ == "__main__":
    main(*sys.argv[1:])
