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
import subprocess
import collections
import contextlib

import yaml
import pysam
from bx.intervals.intersection import IntervalTree

from Mgh.Picard import sam

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    trim_dir = config['directories']['trim']
    align_dir = config['directories']['align']
    count_dir = config['directories']['counts']
    diff_dir = config['directories']['diff']
    for dirname in [trim_dir, align_dir, count_dir, diff_dir]:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    exp_counts = dict()
    for exp in config['experiments']:
        fasta_file = run_trim(exp, config, trim_dir)
        align_bam = align_and_sort(fasta_file, exp, config, align_dir)
        filter_bam = filter_rna(align_bam, config["annotations"]["filter"])
        count_file = counts_of_interest(exp, config, filter_bam,
                config["annotations"]["find"], count_dir)
        exp_counts[exp['name']] = count_file
    for diff in config['differential']:
        differential_analysis(diff, exp_counts, config, diff_dir)

def differential_analysis(diff, exp_counts, config, diff_dir):
    """Prepare count files and do differential analysis.
    """
    count_file = os.path.join(diff_dir, "%s.csv" % diff['name'])
    if not os.path.exists(count_file):
        _write_count_file(diff, exp_counts, count_file)
    diff_file = os.path.join(diff_dir, "%s-diffs.csv" % diff['name'])
    if not os.path.exists(diff_file):
        cl = config["programs"]["diffexp"].split()
        cl += [count_file]
        subprocess.check_call(cl)
    return diff_file

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

def counts_of_interest(exp, config, in_bam, find_bed, count_dir):
    all_sizes = range(config["analysis"]["min_filter_size"],
                      config["analysis"]["max_filter_size"] + 1)
    out_file = os.path.join(count_dir, "%s-counts.csv" % exp['name'])
    if not os.path.exists(out_file):
        pysam.index(in_bam)
        with contextlib.nested(open(out_file, "w"), open(find_bed)) as \
                              (out_handle, in_handle):
            writer = csv.writer(out_handle)
            writer.writerow(["name"] + ["%sbp" % s for s in all_sizes])

            sam_reader = pysam.Samfile(in_bam, "rb")
            sizes = collections.defaultdict(int)
            names_seen = dict()
            for read in sam_reader:
                if not names_seen.has_key(read.qname):
                    sizes[read.rlen] += 1
                    names_seen[read.qname] = ""
            writer.writerow(["total"] + [str(sizes[s]) for s in all_sizes])
            sam_reader.close()

            sam_reader = pysam.Samfile(in_bam, "rb")
            for c, s, e, name in csv.reader(in_handle, dialect="excel-tab"):
                sizes = collections.defaultdict(int)
                for read in sam_reader.fetch(c, int(s), int(e)):
                    sizes[read.rlen] += 1
                writer.writerow([name] + [str(sizes[s]) for s in all_sizes])
            sam_reader.close()
    return out_file

def filter_rna(in_bam, filter_bed):
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
                if read.is_reverse:
                    end = read.pos
                    start = end - read.rlen
                else:
                    start = read.pos
                    end = start + read.rlen
                if len(intervals[chrom].find(start, end)) > 0:
                    to_filter[read.qname] = ""
        in_sam.close()
        in_sam = pysam.Samfile(in_bam, "rb")
        out_sam = pysam.Samfile(out_bam, "wb", template=in_sam)
        remain = dict()
        for read in in_sam:
            if not read.is_unmapped:
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
        sort_bam = sam.samtool_sam_to_sortbam(out_sam,
                config["analysis"]["seq_genome"], config["programs"]["samtools"])
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
