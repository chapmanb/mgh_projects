#!/usr/bin/env python
"""Explore statistics of read clustering between experiments.

Usage:
    cluster_stats.py <config>
"""
import sys
import os
import collections
import itertools
import csv
import math
import multiprocessing

import WPSRC

import yaml
import pysam
import numpy
from bx.intervals.cluster import ClusterTree
import rpy2.robjects as robjects

def main(config_file):
    do_export = True
    do_phasing = False
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    make_dirs = ['cluster']
    for confname in make_dirs:
        dirname = config['directories'][confname]
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    counts = []
    sizes = []
    phasing_info = []
    for exp in config['experiments']:
        bam_file = os.path.join(config['directories']['align'],
                "%s-sort-filter.bam" % exp['name'])
        print "Building clusters for %s" % exp["name"]
        clusters, start_map = build_clusters(bam_file,
                config['analysis']['cluster_distance'],
                config['analysis']['min_cluster_items'])
        if do_export:
            export_clusters(clusters,
                            config['analysis']['min_cluster_export'],
                            config['analysis']['max_cluster_export'],
                            exp['name'],
                            config['directories']['cluster'])
        if do_phasing:
            print "Phasing for %s" % exp["name"]
            phasing = analyze_phasing(clusters,
                            config['analysis']['min_cluster_export'],
                            config['analysis']['max_cluster_export'],
                            exp['name'], start_map)
            phasing_info.append((exp['name'], phasing))
        print "Calculating statistics", exp["name"]
        psizes, pcounts = cluster_stats(clusters)
        counts.append((exp['name'], pcounts))
        sizes.append((exp['name'], psizes))
    plot_percents(counts, "count")
    plot_percents(sizes, "size")
    if do_phasing:
        plot_phasing(phasing_info, config['analysis']['min_cluster_export'],
                     config['analysis']['max_cluster_export'])

def analyze_phasing(clusters, min_size, max_size, name, start_map):
    """Provide a report on phasing information in clusters.

    Utitilizes the WPSRC algorithm of psrnac to do the counting.

    http://code.google.com/p/psrnac/
    """
    pool = multiprocessing.Pool(processes=5)
    cluster_gen = _in_size_clusters(clusters, min_size, max_size, start_map)
    phasings = pool.map(_calc_phasing, cluster_gen)
    print name, len(phasings), numpy.mean(phasings), numpy.median(phasings), \
            numpy.std(phasings)
    return phasings

def _in_size_clusters(clusters, min_size, max_size, start_map):
    tcount = 0
    dcount = 0
    for chrom, cluster in clusters.iteritems():
        for (start, end, ids) in cluster.getregions():
            if (end - start) > min_size and (end - start) <= max_size:
                tcount += 1
                if tcount % 75 == 0 and dcount < 60:
                    dcount += 1
                    starts = sorted(list(set([start_map[i] for i in ids])))
                    yield (start, end, starts)

def _calc_phasing(args):
    start, end, starts = args
    most_phased = int(math.floor((end - start) / 21.0))
    phasing = WPSRC.basic(starts)
    max_phased = max(count for (pos, count) in phasing)
    return float(max_phased) / float(most_phased)

def plot_phasing(phasing, min_size, max_size):
    """Provide a boxplot of phasing summary statistics.
    """
    title = "Relative phasing of %sbp to %sbp regions" % (min_size, max_size)
    out_file = "phasing_by_exp.pdf"
    plot_data = {'exp' : [], 'relative.phasing' : []}
    for exp_name, vals in phasing:
        for val in vals:
            plot_data['exp'].append(exp_name)
            plot_data['relative.phasing'].append(val)
    plot_data['exp'] = robjects.StrVector(plot_data['exp'])
    plot_data['relative.phasing'] = robjects.FloatVector(plot_data['relative.phasing'])
    robjects.r.assign('exp.data', robjects.r['data.frame'](**plot_data))
    robjects.r.assign('save.file', out_file)
    robjects.r.assign('title', title)
    robjects.r('''
      library(ggplot2)
      print(head(exp.data))
      p <- ggplot(exp.data, aes(x=exp, y=relative.phasing))
      p <- p + geom_boxplot()
      p <- p + opts(title=title)
      ggsave(save.file, p, width=11, height=8)
    ''')

def export_clusters(clusters, min_size, max_size, name, out_dir):
    """Export a BED file with cluster locations in a specific size range.
    """
    out_file = os.path.join(out_dir,
            "%s-cluster-%s-%s.bed" % (name, min_size, max_size))
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle, dialect='excel-tab')
        for chrom, cluster in clusters.iteritems():
            for (start, end, _) in cluster.getregions():
                if (end - start) > min_size and (end - start) <= max_size:
                    writer.writerow([chrom, start, end])

def plot_percents(named_vals, ptype):
    out_file = "cluster-%s-summary.pdf" % ptype
    plot_data = {ptype : [], 'percent': [], 'exp' : []}
    for exp_name, vals in named_vals:
        for index, percent in vals.iteritems():
            plot_data[ptype].append(index)
            plot_data['percent'].append(percent)
            plot_data['exp'].append(exp_name)
    plot_data[ptype] = robjects.IntVector(plot_data[ptype])
    plot_data['percent'] = robjects.FloatVector(plot_data['percent'])
    plot_data['exp'] = robjects.StrVector(plot_data['exp'])
    robjects.r.assign('exp.data', robjects.r['data.frame'](**plot_data))
    robjects.r.assign('save.file', out_file)
    robjects.r('''
      library(ggplot2)
      print(head(exp.data))
      p <- ggplot(exp.data, aes(x=%s, y=percent, group=exp))
      p <- p + geom_point(aes(colour=exp)) + stat_smooth(se=FALSE, aes(colour=exp))
      p <- p + xlim(0, 400)
      ggsave(save.file, p, width=11, height=8)
    ''' % ptype)
      #p <- p + facet_grid(exp ~ .)
      #+ scale_x_log10()

def cluster_stats(clusters):
    counts = collections.defaultdict(int)
    sizes = collections.defaultdict(int)
    for chrom, cluster in clusters.iteritems():
        for start, end, all_ids in cluster.getregions():
            counts[len(all_ids)] += 1
            sizes[end - start] += 1
    return _normalize_stats(sizes), _normalize_stats(counts)

def _normalize_stats(info):
    thresh_max = 5
    total = sum(v for v in info.values() if v > thresh_max)
    final = dict()
    for key, val in info.iteritems():
        if val > thresh_max:
            final[key] = float(val) / float(total)
    return final

def build_clusters(bam_file, distance, min_num):
    clusters = collections.defaultdict(lambda: ClusterTree(distance, min_num))
    sam_reader = pysam.Samfile(bam_file, "rb")
    i = 0
    start_map = dict()
    #for read in itertools.islice(sam_reader, 100000):
    for read in sam_reader:
        if not read.is_unmapped:
            chrom = sam_reader.getrname(read.rname)
            if read.is_reverse:
                end = read.pos
                start = max(0, end - read.rlen)
            else:
                start = read.pos
                end = start + read.rlen
            clusters[chrom].insert(start, end, i)
            start_map[i] = start
            i+= 1
    return dict(clusters), start_map

if __name__ == "__main__":
    main(*sys.argv[1:])
