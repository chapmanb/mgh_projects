#!/usr/bin/env python
"""Intersect peaks to generate Venn diagram like stats of overlaps.

Usage:
    intersect_classes.py <config file>
"""
import sys
import subprocess
import collections

import yaml
from bx.intervals.intersection import IntervalTree

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    for space in [None] + config["classify"]["chrs"]:
        display_overlaps(config, space)

def display_overlaps(config, cur_space):
    print cur_space or "all"
    info = []
    for cur_bed in config["classify"]["final"] + \
                   [config["classify"]["positive"]]:
        info.append((cur_bed["name"].replace(" ", "_"),) +
                    parse_bed(cur_bed["file"], cur_space))
    for (name1, bed1, _) in info:
        print name1, len(bed1)
    for i, (name1, bed1, _) in enumerate(info):
        for i2, (name2, _, itree2) in enumerate(info[i+1:]):
            print name1, name2, len(intersect(bed1, itree2))
            for (name3, _, itree3) in info[i+i2+2:]:
                print name1, name2, name3, len(intersect(bed1, itree2, itree3))

    [(name1, bed1, itree1), (name2, bed2, itree2), (_, _, itree3)] = info
    print_unique(name1, bed1, itree2, itree3)
    print_unique(name2, bed2, itree1, itree3)

def print_unique(name, bed, itree_o, itree_c):
    all_ol = intersect(bed, itree_o, itree_c)
    cnl_ol = intersect(bed, itree_c)
    unique = sorted(list(set(cnl_ol) - set(all_ol)))
    print name
    for space, start, end in unique:
        print space, start, end

def intersect(one, two, three=None, union=True):
    final = []
    for space, start, end in one:
        two_matches = two[space].find(start, end)
        if _passes(two_matches, union):
            if three:
                three_matches = three[space].find(start, end)
                if _passes(three_matches, union):
                    final.append((space, start, end))
            else:
                final.append((space, start, end))
    return final

def _passes(matches, union=True):
    return (union and len(matches) > 0 or
            not union and len(matches) == 0)

def parse_bed(in_file, required_space=None):
    itree = collections.defaultdict(IntervalTree)
    regions = []
    with open(in_file) as in_handle:
        for space, start, end in (l.split()[:3] for l in in_handle):
            if required_space is None or space == required_space:
                start = int(start)
                end = int(end)
                regions.append((space, start, end))
                itree[space].insert(start, end)
    return regions, itree

if __name__ == "__main__":
    main(*sys.argv[1:])
