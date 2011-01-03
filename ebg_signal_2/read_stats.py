#!/usr/bin/env python
"""Provide graphs of read sizes and nucleotide distribution across reads.

This provides high level summaries of information on aligning reads.
"""
import sys
import os
import itertools
import collections

import yaml
import pysam
import numpy
import rpy2.robjects as robjects

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    exp_sizes = dict()
    exp_bases = dict()
    for exp in config['experiments']:
        print exp['name']
        bam_file = os.path.join(config['directories']['align'],
                "%s-sort-filter.bam" % exp['name'])
        cur_size, cur_bases = size_base_stats(bam_file)
        exp_sizes[exp['name']] = cur_size
        exp_bases[exp['name']] = cur_bases
    plot_bases(exp_bases)
    plot_sizes(exp_sizes)

def plot_bases(bases):
    """Plot the distribution of bases at each position.
    """
    out_file = "base_distribution.pdf"
    plot_data = {"base": [], "position" : [], "count" : [], "experiment" : []}
    for exp_name, pos_info in bases.iteritems():
        for pos, base_info in pos_info.iteritems():
            for base, count in base_info.iteritems():
                plot_data['base'].append(base)
                plot_data['position'].append(pos)
                plot_data['experiment'].append(exp_name)
                plot_data['count'].append(count)
    for name, convert_to in [('base', robjects.StrVector),
                             ('count', robjects.IntVector),
                             ('position', robjects.IntVector),
                             ('experiment', robjects.StrVector)]:
        plot_data[name] = convert_to(plot_data[name])
    robjects.r.assign('exp.data', robjects.r['data.frame'](**plot_data))
    robjects.r.assign('save.file', out_file)
    robjects.r('''
      library(ggplot2)
      print(head(exp.data))
      p <- ggplot(exp.data, aes(x=position, y=count, fill=base))
      p <- p + geom_bar(position="fill", stat="identity")
      p <- p + facet_wrap(~experiment)
      ggsave(save.file, p, width=11, height=8)
    ''')

def plot_sizes(sizes):
    """Plot the distribution of RNA sizes per experiment.
    """
    out_file = "size_distribution.pdf"
    plot_data = {'size' : [], 'experiment': [], 'count' : []}
    exp_totals = dict()
    for exp_name, exp_sizes in sizes.iteritems():
        total = 0
        for _, count in exp_sizes.iteritems():
            total += count
        exp_totals[exp_name] = float(total)

    for exp_name, exp_sizes in sizes.iteritems():
        for size, count in exp_sizes.iteritems():
            plot_data['size'].append(size)
            plot_data['experiment'].append(exp_name)
            plot_data['count'].append(count / exp_totals[exp_name])
    for name, convert_to in [('size', robjects.IntVector),
                             ('count', robjects.FloatVector),
                             ('experiment', robjects.StrVector)]:
        plot_data[name] = convert_to(plot_data[name])
    robjects.r.assign('exp.data', robjects.r['data.frame'](**plot_data))
    robjects.r.assign('save.file', out_file)
    robjects.r('''
      library(ggplot2)
      print(head(exp.data))
      p <- ggplot(exp.data, aes(x=size, y=count))
      p <- p + geom_bar(stat="identity")
      p <- p + facet_wrap(~experiment)
      ggsave(save.file, p, width=11, height=8)
    ''')

def size_base_stats(bam_file):
    """Count read size stats and per base GATC distribution.
    """
    sam_reader = pysam.Samfile(bam_file, "rb")
    sizes = collections.defaultdict(int)
    bases = collections.defaultdict(lambda : collections.defaultdict(int))
    #for read in itertools.islice(sam_reader, 100000):
    for read in sam_reader:
        if not read.is_unmapped:
            sizes[len(read.seq)] += 1
            for pos, base in enumerate(read.seq):
                bases[pos+1][base] += 1
    return dict(sizes), dict(bases)

if __name__ == "__main__":
    main(*sys.argv[1:])
