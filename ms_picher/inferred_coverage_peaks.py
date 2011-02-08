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
import collections

import yaml
import pysam
from bx.intervals.intersection import IntervalTree
from bx.seq.twobit import TwoBitFile
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from bcbio.picard import PicardRunner
from bcbio.picard.utils import chdir

# ## High level functions to drive process

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    if config.get("align_dir", None) is not None:
        align_and_call(config, config_file)
    else:
        call_from_alignments(config)

def align_and_call(config, config_file):
    """Do peak calls and alignments, starting with fastq files.
    """
    picard = PicardRunner(config["program"]["picard"])
    align_dir = config["align_dir"]
    if not os.path.exists(align_dir):
        os.makedirs(align_dir)
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
        call_peaks(bed_files, ref, config)

def call_from_alignments(config):
    """Do peak calls, starting with alignment files.
    """
    bed_files = {}
    for align in config["alignments"]:
        bed_files[align["name"]] = bam_to_bed(align["file"])
    for ref in config["ref"]:
        call_peaks(bed_files, ref, config)

def call_peaks(bed_files, ref, config):
    """Top level peak calling functionality.
    """
    tmp_dir = config["tmp_dir"]
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    xreact_bases = config["algorithm"]["cross_react_bases"]
    min_width = config["algorithm"].get("min_width", 0)
    peaks = []
    oligos = []
    for peak in config["peaks"]:
        peaks.append(call_macs_peaks(bed_files[peak["condition"]],
                                     bed_files[peak["background"]],
                                     peak["condition"], ref,
                                     config["peak_dir"]))
        oligos.append(peak["oligos"])
    peak_screen = (PeakCrossReactScreener(ref, xreact_bases, min_width, tmp_dir)
                   if xreact_bases > 0 else None)
    ol_percent = config["algorithm"].get("overlap_percent", None)
    if ol_percent:
        output_combine_peaks(peaks, oligos, peak_screen, ol_percent)
    else:
        for i, peak in enumerate(peaks):
            filter_peaks(peak, oligos[i], peak_screen)

# ## Combine peaks from all experiments
# Read peaks from a MACS style BED file for all experiments and
# generate an output file with only peaks that overlap from all experiments.
# This results in a smaller set of final peaks for easier interpretation.

class PeakOverlapper:
    """Keep tally of peaks in base intervals that overlap with comparisons.
    """
    def __init__(self, compares, chrom, overlap_percent, peak_screener,
                 oligos):
        """Initialize with a list of IntervalTrees for comparison.
        """
        self._compares = compares
        self._pct = overlap_percent
        self._chrom = chrom
        self._screener = peak_screener
        self._oligos = oligos
        self.final = []

    def check_region(self, node):
        cmp_bases = set(range(node.start, node.end))
        count = 0
        for cmptree in self._compares:
            overlaps = cmptree.find(node.start, node.end)
            for ol in overlaps:
                if self._current_overlap(ol, cmp_bases) >= self._pct:
                    count += 1
                    break
        if count == len(self._compares):
            if (self._screener is None or
                self._screener.check(self._chrom, node.start, node.end, self._oligos)):
                self.final.append((node.start, node.end))

    def _current_overlap(self, ol, cmp_bases):
        ol_bases = set(range(ol['start'], ol['end']))
        return len(ol_bases & cmp_bases) / float(len(cmp_bases))

def _combine_peaks(peak_files, oligos, peak_screener, overlap_percent):
    """Retrieve only peaks that overlap in all three experiments.
    """
    itrees = [_bed_to_intervaltree(f) for f in peak_files]
    for chrom, base_itree in itrees[0].iteritems():
        overlapper = PeakOverlapper([t[chrom] for t in itrees[1:]], chrom,
                                    overlap_percent, peak_screener, oligos[0])
        base_itree.traverse(overlapper.check_region)
        for start, end in overlapper.final:
            yield chrom, start, end

def output_combine_peaks(peak_files, oligos, peak_screener, overlap_percent):
    out_file = "%s-combined%s" % os.path.splitext(peak_files[0])
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for chrom, start, end in _combine_peaks(peak_files, oligos,
                                                    peak_screener, overlap_percent):
                writer.writerow([chrom, start, end])

def filter_peaks(peak_file, oligos, peak_screener):
    out_file = "%s-filter%s" % os.path.splitext(peak_file)
    if not os.path.exists(out_file):
        with open(peak_file) as in_handle:
            with open(out_file, "w") as out_handle:
                writer = csv.writer(out_handle, dialect="excel-tab")
                for chrom, start, end in _read_bed(in_handle):
                    if peak_screener.check(chrom, start, end, oligos):
                        writer.writerow([chrom, start, end])

def _read_bed(in_handle):
    reader = csv.reader(in_handle, dialect="excel-tab")
    for chrom, start, end in ((l[0], int(l[1]), int(l[2])) for l in reader):
        yield chrom, start, end

def _bed_to_intervaltree(bed_file):
    itree = collections.defaultdict(IntervalTree)
    with open(bed_file) as in_handle:
        reader = csv.reader(in_handle, dialect="excel-tab")
        for chrom, start, end in ((l[0], int(l[1]), int(l[2])) for l in reader):
            itree[chrom].insert(start, end, dict(start=start, end=end))
    return itree

# ## Remove peak regions where DNA sequence matches binding oligos
# A potential source of pull down noise are DNA segments that overlap with
# the oligos used in the pull down experiment. This detects and removes these
# regions to yield a set of clean peaks for overlap experiments.

class PeakCrossReactScreener:
    def __init__(self, ref, xreact_bases, min_width, tmp_dir):
        self._2bit = self._get_twobit(ref)
        self._xreact = xreact_bases
        self._min_width = min_width
        self._tmpdir = tmp_dir

    def _get_twobit(self, ref):
        """Retrieve bx-python 2bit class from /reference hierarchy.
        """
        fname = os.path.join(os.path.dirname(ref), os.pardir, "ucsc",
                             "%s.2bit" % os.path.basename(ref))
        return TwoBitFile(open(fname))

    def check(self, chrom, start, end, oligos):
        # check for too small regions
        if self._min_width and end-start < self._min_width:
            return False
        # check for regions containing oligo matches
        seq = self._2bit[chrom].get(start, end)
        blast_db = self._make_blastdb(seq)
        input_file = self._make_input_file(oligos)
        blast_out = os.path.join(self._tmpdir, "out.blast")
        cl = NcbiblastnCommandline(query=input_file, db=blast_db,
                                   out=blast_out, outfmt=5, num_alignments=1,
                                   num_descriptions=0, task="blastn-short")
        subprocess.check_call(str(cl).split())
        with open(blast_out) as blast_handle:
            for rec in NCBIXML.parse(blast_handle):
                if self._has_oligo_match(rec):
                    return False
        return True

    def _has_oligo_match(self, rec):
        """Determine if the BLAST record indicates matches to an input oligo.
        """
        for align in rec.alignments:
            for hsp in align.hsps:
                if hsp.identities >= self._xreact:
                    return True
        return False

    def _make_input_file(self, oligos):
        out_file = os.path.join(self._tmpdir, "in.fa")
        with open(out_file, "w") as out_handle:
            for i, oligo in enumerate(oligos):
                out_handle.write(">o%s\n%s\n" % (i, oligo))
        return out_file

    def _make_blastdb(self, seq):
        out_file = "ref.fa"
        with chdir(self._tmpdir):
            with open(out_file, "w") as out_handle:
                out_handle.write(">ref\n%s\n" % seq.upper())
            cl = ["makeblastdb", "-in", out_file, "-dbtype", "nucl",
                  "-out", out_file]
            with open("/dev/null", "w") as stdout:
                subprocess.check_call(cl, stdout=stdout)
        return os.path.join(self._tmpdir, out_file)

# ## Generate inferred coverage BED file

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

# ## Run external programs
# These functions drive various external programs (bowtie, MACS) and handle
# format conversions.

def call_macs_peaks(exp_bed, back_bed, out_base, ref, peak_dir):
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    genome_size = _size_from_fai(_get_ref_fai(ref))
    out_name = os.path.join(peak_dir, "%s-%s-macs" %
                            (out_base.replace(" ", "_"), os.path.basename(ref)))
    peak_file = "%s_peaks.bed" % out_name
    if not os.path.exists(peak_file):
        cl = ["macs14", "-t", exp_bed, "-c", back_bed, "--name=%s" % out_name,
              "--format=BED", "--gsize=%s" % genome_size]
        subprocess.check_call(cl)
    return peak_file

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
