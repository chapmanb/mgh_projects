#!/usr/bin/env python
"""Trim, filter, align and calculate differential expression of short RNA inputs.

Usage:
    shortrna_analysis.py <config>

Where the config file is a YAML configuration file specifying experiments and
directories for analysis. This is the base reusable script for performing the
analysis.
"""
import sys
import os
import csv
import math
import subprocess
import collections
import contextlib

import yaml
import pysam
from bx.intervals.intersection import IntervalTree

from bcbio.pipeline.alignment import sam_to_sort_bam

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    trim_dir = config['directories']['trim']
    align_dir = config['directories']['align']
    count_dir = config['directories']['counts']
    diff_dir = config['directories']['diff']
    all_sizes = range(config["analysis"]["min_filter_size"],
                      config["analysis"]["max_filter_size"] + 1)
    for dirname in [trim_dir, align_dir, count_dir, diff_dir]:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    exp_counts = dict()
    for approach in config['analysis']['diff_approaches']:
        exp_counts[approach] = dict()
        for exp in config['experiments']:
            fasta_file = run_trim(exp, config, trim_dir)
            align_bam = align_and_sort(fasta_file, exp, config, align_dir)
            filter_bam = filter_rna(align_bam, config["annotations"]["filter"],
                                    all_sizes)
            count_file = generate_exp_counts(approach, exp, config, all_sizes,
                                             filter_bam, count_dir)
            exp_counts[approach][exp['name']] = count_file
    for diff in config['differential']:
        differential_analysis(diff, exp_counts, config, diff_dir)

def differential_analysis(diff, exp_counts, config, diff_dir):
    """Prepare count files and do differential analysis.
    """
    if diff["algorithm"] == "edgeR":
        _diff_analysis_edgeR(diff, exp_counts, config, diff_dir)
    elif diff["algorithm"] == "fold_change":
        out_file = os.path.join(diff_dir, "%s-%s-diffs.csv" %
                                (diff['name'], diff['approach']))
        if not os.path.exists(out_file):
            _diff_fold_change(diff, exp_counts, config, out_file)
    else:
        raise NotImplementedError(diff["algorithm"])

def _diff_fold_change(diff, exp_counts, config, out_file):
    """Identify differences based on a fold change difference.
    """
    assert len(diff["control"]) == 1
    assert len(diff["experimental"]) == 1
    thresh = float(config["analysis"]["fold_change"])
    cname = diff["control"][0]
    ename = diff["experimental"][0]
    c_counts = _read_count_file(exp_counts[diff["approach"]][cname])
    e_counts = _read_count_file(exp_counts[diff["approach"]][ename])
    all_regions = sorted(list(set(c_counts.keys() + e_counts.keys())))
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["chrom", "start", "end", "seq", cname, ename, "fold-change"])
        for region in all_regions:
            c_score = c_counts.get(region, 1.0)
            e_score = e_counts.get(region, 1.0)
            diff = c_score / e_score
            if diff > thresh or diff < (1.0 / thresh):
                writer.writerow(list(region) + [c_score, e_score, "%.1f" % diff])

def _read_count_file(in_file):
    counts = {}
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        reader.next() # header
        for chrom, start, end, seq, count in reader:
            counts[(chrom, start, end, seq)] = float(count)
    return counts

def _diff_analysis_edgeR(diff, exp_counts, config, diff_dir):
    """Calculate differential expression counts using edgeR.
    """
    count_file = os.path.join(diff_dir, "%s-%s.csv" % (diff['name'],
                                                       diff['approach']))
    if not os.path.exists(count_file):
        _write_count_file(diff, exp_counts[diff['approach']], count_file)
    diff_file = "%s-diffs%s" % os.path.splitext(count_file)
    if not os.path.exists(diff_file):
        cl = config["programs"]["diffexp"].split()
        cl += [count_file]
        subprocess.check_call(cl)

def _write_count_file(diff, exp_counts, count_file):
    with open(count_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(['name'] + diff['control'] + diff['experimental'])
        writer.writerow([''] + ['1' for _ in diff['control']] +
                               ['2' for _ in diff['experimental']])
        readers = [csv.reader(open(exp_counts[e])) for
                    e in diff['control'] + diff['experimental']]
        [r.next() for r in readers] # header
        for parts in readers[0]:
            out_info = [parts[0]]
            for line in [parts] + [r.next() for r in readers[1:]]:
                out_info.append(sum(int(i) for i in line[1:]))
            writer.writerow([str(i) for i in out_info])

def generate_exp_counts(approach, exp, config, all_sizes, in_bam, count_dir):
    """Generate count file for an experiment using the specific approach.
    """
    out_file = os.path.join(count_dir, "%s-%s-counts.csv" %
                            (exp['name'], approach))
    if not os.path.exists("%s.bai" % in_bam):
        pysam.index(in_bam)
    if not os.path.exists(out_file):
        if approach == "known":
            counts_of_interest(config, in_bam, config["annotations"]["find"],
                               all_sizes, out_file)
        elif approach == "read":
            raw_read_counts(config, in_bam, all_sizes, out_file)
        else:
            raise NotImplementedError(approach)
    return out_file

def raw_read_counts(config, in_bam, all_sizes, out_file):
    """Prepare counts for individual non-combined reads.

    Reads above the mean threshold are returned, normalized
    by reads per million.
    """
    sam_rdr = pysam.Samfile(in_bam, "rb")
    region_counts = collections.defaultdict(int)
    for i, read in enumerate(sam_rdr):
        if read.rlen in all_sizes:
            loc = (sam_rdr.getrname(read.rname), read.pos, read.pos + read.rlen,
                   read.seq)
            region_counts[loc] += 1
    total = float(sum(_totals_by_size(in_bam, all_sizes)))
    with open(out_file, "w") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["space", "start", "end", "seq", "count"])
        regions = region_counts.keys()
        regions.sort()
        for (space, start, end, seq) in regions:
            val = region_counts[(space, start, end, seq)]
            norm_val = "%.1f" % (val / total * 1e6)
            writer.writerow([space, start, end, seq, norm_val])

def counts_of_interest(config, in_bam, find_bed, all_sizes, out_file):
    """Provide count information in defined regions of interest from a BED file.
    """
    with contextlib.nested(open(out_file, "w"), open(find_bed)) as \
                          (out_handle, in_handle):
        writer = csv.writer(out_handle)
        writer.writerow(["name"] + ["%sbp" % s for s in all_sizes])
        writer.writerow(["total"] +
                        [str(s) for s in _totals_by_size(in_bam, all_sizes)])
        sam_reader = pysam.Samfile(in_bam, "rb")
        for c, s, e, name in csv.reader(in_handle, dialect="excel-tab"):
            sizes = collections.defaultdict(int)
            for read in sam_reader.fetch(c, int(s), int(e)):
                sizes[read.rlen] += 1
            writer.writerow([name] + [str(sizes[s]) for s in all_sizes])
        sam_reader.close()

def _totals_by_size(in_bam, all_sizes):
    """Get total read counts for a BAM file organized by size.
    """
    sam_reader = pysam.Samfile(in_bam, "rb")
    sizes = collections.defaultdict(int)
    names_seen = dict()
    for read in sam_reader:
        if not names_seen.has_key(read.qname):
            sizes[read.rlen] += 1
            names_seen[read.qname] = ""
    sam_reader.close()
    return [sizes[s] for s in all_sizes]

def filter_rna(in_bam, filter_bed, filter_sizes=None):
    """Produce a new BAM file with only reads aligning to non-filtered regions.
    """
    out_bam = "%s-filter%s" % os.path.splitext(in_bam)
    if not os.path.exists(out_bam):
        to_filter = dict()
        intervals = _filter_intervals(filter_bed)
        in_sam = pysam.Samfile(in_bam, "rb")
        for read in in_sam:
            if not read.is_unmapped:
                chrom = in_sam.getrname(read.rname)
                start, end = (read.pos, read.pos + read.rlen)
                if len(intervals[chrom].find(start, end)) > 0:
                    to_filter[read.qname] = ""
        in_sam.close()
        in_sam = pysam.Samfile(in_bam, "rb")
        out_sam = pysam.Samfile(out_bam, "wb", template=in_sam)
        remain = dict()
        for read in in_sam:
            if not read.is_unmapped:
                if (filter_sizes is None or read.rlen in filter_sizes):
                    if not to_filter.has_key(read.qname):
                        remain[read.qname] = ""
                        out_sam.write(read)
        in_sam.close()
        out_sam.close()
        print out_bam, "Filtered", len(to_filter), "Remain", len(remain)
    return out_bam

def _filter_intervals(in_bed):
    trees = collections.defaultdict(lambda: IntervalTree())
    with open(in_bed) as in_handle:
        for parts in csv.reader(in_handle, dialect="excel-tab"):
            chrom, start, end, name = parts[:4]
            trees[chrom].insert(int(start), int(end), name)
    return dict(trees)

def align_and_sort(in_file, exp, config, align_dir):
    out_sam = os.path.join(align_dir, "%s.sam" % (os.path.splitext(
        os.path.basename(in_file))[0]))
    out_bam = "%s.bam" % os.path.splitext(out_sam)[0]
    sort_bam = "%s-sort.bam" % os.path.splitext(out_sam)[0]
    if not os.path.exists(sort_bam):
        if not os.path.exists(out_sam):
            cl = [config["programs"]["bowtie"],
                  "-S", "-f", "--all", "--best", "--strata",
                  "-v", str(config["analysis"]["align_errors"]),
                  config["analysis"]["bowtie_genome"],
                  in_file, out_sam]
            print cl
            subprocess.check_call(cl)
        sort_bam = sam_to_sort_bam(out_sam, config["analysis"]["seq_genome"],
                                   in_file, None, exp["name"], "", exp["name"],
                                   config)
        for to_remove in [out_sam, out_bam]:
            if os.path.exists(to_remove):
                os.remove(to_remove)
    return sort_bam

def run_trim(exp, config, out_dir):
    out_file = os.path.join(out_dir, "%s.fa" % exp['name'])
    if not os.path.exists(out_file):
        cl = config['programs']['adaptor_trim'].split()
        cl += [exp['fastq'], out_file, exp['adaptor'],
               str(config['analysis']['trim_errors']),
               "--min_size=%s" % config['analysis']['min_filter_size'],
               "--max_size=%s" % config['analysis']['max_filter_size']]
        print cl
        subprocess.check_call(cl)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
