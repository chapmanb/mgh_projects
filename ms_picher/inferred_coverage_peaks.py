#!/usr/bin/env python
"""Prepare peak calls for PIChER reads using inferred coverage.

Driver script that for each experiment:
- Does bowtie alignments
- Prepares inferred coverage BED files from paired end alignments
- Uses BED files to predict peaks with MACS
- Provides a BigWig inferred coverage display file

Usage:
  inferred_coverage_peaks.py <YAML configuration file>
"""
import os
import sys
import csv
import subprocess

import yaml
import pysam

from bcbio.picard import PicardRunner

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    align_dir = config["align_dir"]
    if not os.path.exists(align_dir):
        os.makedirs(align_dir)
    picard = PicardRunner(config["program"]["picard"])
    for ref in config["ref"]:
        bed_files = {}
        for cur_fastq in config["fastq"]:
            fastq_one, fastq_two = _get_align_files(cur_fastq, config)
            base = "%s-%s" % (cur_fastq["name"], os.path.basename(ref))
            align_file = bowtie_to_sam(fastq_one, fastq_two, ref, base,
                                       align_dir, config)
            bam_file = sam_to_bam(align_file, fastq_one, fastq_two, ref,
                                  config_file)
            bed_files[cur_fastq["name"]] = bam_to_bed(bam_file)
            inferred_bed = generate_inferred_coverage(bam_file)
            generate_bigwig(inferred_bed, ref, picard)
        peaks = []
        for peak in config["peaks"]:
            peaks.append(call_macs_peaks(bed_files[peak["condition"]],
                                         bed_files[peak["background"]],
                                         peak["condition"], ref,
                                         config["peak_dir"]))
        print peaks

def call_macs_peaks(exp_bed, back_bed, out_base, ref, peak_dir):
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    genome_size = _size_from_fai(_get_ref_fai(ref))
    out_name = os.path.join(peak_dir,
                            "%s-%s-macs" % (out_base, os.path.basename(ref)))
    cl = ["macs14", "-t", exp_bed, "-c", back_bed, "--name=%s" % out_name,
          "--format=BED", "--gsize=%s" % genome_size]
    subprocess.check_call(cl)
    return out_name

def generate_bigwig(bed_file, ref, picard):
    genome_ref = _get_ref_fai(ref)
    bam_file = "%s.bam" % os.path.splitext(bed_file)[0]
    if not os.path.exists(bam_file):
        cl = ["bedToBam", "-i", bed_file, "-g", genome_ref]
        with open(bam_file, "wb") as out_handle:
            subprocess.check_call(cl, stdout=out_handle)
    sort_bam = picard.run_fn("picard_sort", bam_file)
    cl = ["bam_to_wiggle.py", sort_bam]
    subprocess.check_call(cl)

def _size_from_fai(index_file):
    total = 0
    with open(index_file) as in_handle:
        for line in in_handle:
            parts = line.split()
            total += int(parts[1])
    return total

def _get_ref_fai(ref):
    base_dir, name = os.path.split(ref)
    return os.path.join(base_dir, os.pardir, "seq", "%s.fa.fai" % name)

def generate_inferred_coverage(bam_file):
    """Create an inferred coverage BED file from paired end regions.
    """
    bed_file = "%s-inferred.bed" % os.path.splitext(bam_file)[0]
    if not os.path.exists(bed_file):
        with open(bed_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            samfile = pysam.Samfile(bam_file, "rb")
            for read in samfile.fetch():
                if (read.is_paired and read.is_read1 and not read.is_unmapped
                    and not read.mate_is_unmapped and read.rname == read.mrnm):
                    chrom = samfile.getrname(read.rname)
                    pos = [read.pos, read.pos + read.rlen, read.mpos,
                           read.mpos + read.rlen]
                    writer.writerow([chrom, min(pos), max(pos), read.qname])
    return bed_file

def bam_to_bed(bam_file):
    bed_file = "%s.bed" % os.path.splitext(bam_file)[0]
    if not os.path.exists(bed_file):
        with open(bed_file, "w") as out_handle:
            cl = ["bamToBed", "-i", bam_file]
            subprocess.check_call(cl, stdout=out_handle)
    return bed_file

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
