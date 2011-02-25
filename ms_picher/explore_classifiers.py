#!/usr/bin/env python
"""Explore classifiers for separating positives from a set of peaks.

This uses a list of statistics along with a group of positive regions
to build a classifier tree differentiating the positive and negative results.
The output is a list of classifiers to use.

Usage:
  explore_classifiers.py <config.yaml>
"""
import os
import sys
import csv
import collections

import yaml
import tablib
import rpy2.robjects as rpy
from bx.intervals.intersection import IntervalTree

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    cconfig = config["classify"]
    with open(cconfig["positive"]) as in_handle:
        pos_itree = _peaks_to_intervaltree(in_handle)
    for exp in cconfig["exps"]:
        print exp
        with open(exp) as in_handle:
            data = tablib.Dataset()
            tablib.formats.tsv.import_set(data, in_handle.read())
        data.append(col= _get_groups(data, pos_itree), header="group")
        build_classifier(_dataset_to_rpy_dict(data), _group_prior(data))

def _group_prior(data):
    counts = collections.defaultdict(int)
    for item in data["group"]:
        counts[item] += 1
    return float(counts["no"]) / (float(counts["yes"]) + float(counts["no"]))

def build_classifier(tbl, yes_prior):
    rpy.r.assign("exp.data", rpy.r["data.frame"](**tbl))
    rpy.r.assign("yes.prior", yes_prior)
    rpy.r('''
      library(rpart)
      #library(RWeka)
      #library(randomForest)
      print(head(exp.data))
      #fit <- DecisionStump(group ~ readcount + width + blastoligo,
      #                     data=exp.data)
      #print(fit)
      #fit <- randomForest(group ~ readcount + width + blastoligo,
      #                    data=exp.data, sampsize=c(yes=100, no=100))
      fit <- rpart(group ~ readcount + width + blastoligo,
                   data=exp.data, method="class",
                   parms=list(prior=c(yes=yes.prior, no=1-yes.prior)))
      print(summary(fit))
    ''')

def _dataset_to_rpy_dict(data):
    out = {}
    for i, header in enumerate(data.headers):
        try:
            float(data[0][i])
            rtype = rpy.FloatVector
        except ValueError:
            rtype = rpy.StrVector
        out[header] = rtype([item[i] for item in data])
    return out

def _get_groups(data, pos_itree):
    groups = []
    for (space, start, end) in (l[:3] for l in data):
        try:
            ols = pos_itree[space].find(int(start), int(end))
        except KeyError:
            ols = []
        groups.append("yes" if len(ols) > 0 else "no")
    return groups

def _read_bed(in_handle):
    reader = csv.reader(in_handle, dialect="excel-tab")
    for chrom, start, end in ((l[0], int(l[1]), int(l[2])) for l in reader):
        yield dict(chrom=chrom, start=max(start, 1), end=end)

def _peaks_to_intervaltree(peak_file):
    itree = collections.defaultdict(IntervalTree)
    for peak in _read_bed(peak_file):
        itree[peak["chrom"]].insert(peak["start"], peak["end"], peak)
    return itree

if __name__ == "__main__":
    main(*sys.argv[1:])
