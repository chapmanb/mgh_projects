#!/usr/bin/env python
"""Prepare peak calls for PIChER reads, along with inferred coverage displays.

Driver script that for each experiment:
- Does unique Bowtie alignments.
- Prepares inferred coverage BigWig files from paired end alignments.
- Uses alignments to predict peaks with MACS and BayesPeak.
- Overlaps peak predictions and classify using peak characteristics.

Usage:
  inferred_coverage_peaks.py <YAML configuration file>

Requires:
  Bowtie http://bowtie-bio.sourceforge.net/
  MACS http://liulab.dfci.harvard.edu/MACS/
  R with BayesPeak http://www.bioconductor.org/packages/release/bioc/html/BayesPeak.html
  Picard http://picard.sourceforge.net/
  pysam http://code.google.com/p/pysam/
  rpy2 http://rpy.sourceforge.net/rpy2.html
  bx-python https://bitbucket.org/james_taylor/bx-python
  Biopython http://biopython.org
"""
import os
import sys
import csv
import glob
import copy
import subprocess
import collections
import contextlib
import multiprocessing

import yaml
import pysam
import rpy2.robjects as rpy
from bx.intervals.intersection import IntervalTree
from bx.seq.twobit import TwoBitFile
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from bcbio.picard import PicardRunner
from bcbio.picard.utils import chdir
from bcbio.bam.counts import NormalizedBam

# ## High level functions to drive process

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    if config.get("align_dir", None) is not None:
        align_and_call(config, config_file)
    else:
        call_from_alignments(config)
    for cconfig in config.get("custom_peak_combine", []):
        post_combine_peaks(config, cconfig)

def align_and_call(config, config_file):
    """Do peak calls and alignments, starting with fastq files.
    """
    for ref in config["ref"]:
        map_fn = _multi_map_fn(config)
        bed_info = map_fn(_align_and_bigwig,
                          [(ref, cur_fastq, config, config_file) for cur_fastq in config["fastq"]])
        bed_files = {}
        for name, bed in bed_info:
            bed_files[name] = bed
        call_multi_peaks(bed_files, ref, config)

def _multi_map_fn(config):
    """Get a mapping function, potentially for multiple cores.
    """
    cores = int(config["algorithm"].get("num_cores", 1))
    if cores > 1:
        pool = multiprocessing.Pool(cores)
        return pool.map
    else:
        return map

def _align_and_bigwig(args):
    def do_work(ref, cur_fastq, config, config_file):
        picard = PicardRunner(config["program"]["picard"])
        align_dir = config["align_dir"]
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        fastq_one, fastq_two = _get_align_files(cur_fastq, config)
        base = "%s-%s" % (cur_fastq["name"], os.path.basename(ref))
        align_file = bowtie_to_sam(fastq_one, fastq_two, ref, base,
                                   align_dir, config)
        bam_file = sam_to_bam(align_file, fastq_one, fastq_two, ref,
                              config_file)
        picard.run_fn("picard_index", bam_file)
        bed_file = bam_to_bed(bam_file)
        if fastq_two is not None:
            inferred_bed = generate_inferred_coverage(bam_file)
            generate_bigwig(inferred_bed, ref, picard)
        else:
            generate_bigwig(bed_file, ref, picard)
        return cur_fastq["name"], bed_file
    return do_work(*args)

def call_from_alignments(config):
    """Do peak calls, starting with alignment files.
    """
    bed_files = {}
    for align in config["alignments"]:
        bed_files[align["name"]] = bam_to_bed(align["file"])
    for ref in config["ref"]:
        call_multi_peaks(bed_files, ref, config)

def call_multi_peaks(bed_files, ref, config):
    peak_algorithms = config["algorithm"]["peaks"]
    if isinstance(peak_algorithms, dict):
        peak_algorithms = [peak_algorithms]
    for peak_algorithm in peak_algorithms:
        cur_config = copy.deepcopy(config)
        cur_config["algorithm"] = peak_algorithm
        call_peaks(bed_files, ref, cur_config)

def call_peaks(bed_files, ref, config):
    """Top level peak calling functionality.
    """
    tmp_dir = config["tmp_dir"]
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    picard = PicardRunner(config["program"]["picard"])
    xreact_bases = config["algorithm"].get("cross_react_bases", 0)
    min_width = config["algorithm"].get("min_width", 0)
    min_height = config["algorithm"].get("min_height", 0)
    min_count = config["algorithm"].get("min_count", 0)
    min_score = config["algorithm"].get("min_score", 0)
    min_reads = config["algorithm"].get("min_reads", 0)
    output = config["algorithm"].get("output", "bed")
    peak_caller = config["algorithm"].get("peak_caller", "macs")
    ol_percent = config["algorithm"].get("overlap_percent", None)
    peak_call_fn = _peak_caller_fn(peak_caller)
    peaks = []
    oligos = []
    bam_files = []
    for peak in config["peaks"]:
        peaks.append(peak_call_fn(bed_files[peak["condition"]],
                                  bed_files[peak["background"]],
                                  peak["condition"], ref,
                                  config))
        oligos.append(peak.get("oligos", []))
        bam_file = "%s.bam" % os.path.splitext(bed_files[peak["condition"]])[0]
        if not os.path.exists(bam_file):
            bam_file = None
        else:
            bam_file = NormalizedBam(bam_file, bam_file, picard)
        bam_files.append(bam_file)
    peak_screen = PeakCrossReactScreener(ref, xreact_bases, min_width, min_height,
                                         min_count, min_score, min_reads, output=="bed",
                                         tmp_dir)
    if ol_percent:
        output_combine_peaks(peaks, oligos, bam_files, peak_screen, ol_percent, config)
    else:
        for i, peak in enumerate(peaks):
            filter_peaks(peak, oligos[i], bam_files[i], peak_screen, config)

def _peak_caller_fn(name):
    fns = {"pics" : call_pics_peaks,
           "macs" : call_macs_peaks,
           "bayespeak" : call_bayespeak_peaks}
    return fns[name]

def post_combine_peaks(config, cconfig):
    """Combine peaks from multiple peak call methods.
    """
    peak_dir = config["peak_dir"]
    tmp_dir = config["tmp_dir"]
    final_file = os.path.join(peak_dir, cconfig["name"])
    overlaps = [os.path.join(peak_dir, f) for f in cconfig["overlap"]]
    cur_file = overlaps[0]
    for i, overlap in enumerate(overlaps[1:]):
        out_file = os.path.join(tmp_dir, "%s.o%s" % (os.path.basename(overlaps[0]), i))
        cl = ["bed_intersect.py", cur_file, overlap]
        with open(out_file, "w") as out_handle:
            subprocess.check_call(cl, stdout=out_handle)
        cur_file = out_file
    for i, avoid in enumerate(cconfig.get("avoid", [])):
        out_file = os.path.join(tmp_dir, "%s.a%s" % (os.path.basename(overlaps[0]), i))
        cl = ["bed_intersect.py", "--reverse", cur_file, avoid]
        with open(out_file, "w") as out_handle:
            subprocess.check_call(cl, stdout=out_handle)
        cur_file = out_file
    os.rename(cur_file, final_file)
    return out_file

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
            screen_peak = self._screener.check(node.interval, self._oligos)
            if screen_peak:
                self.final.append(screen_peak)

    def _current_overlap(self, ol, cmp_bases):
        ol_bases = set(range(ol['start'], ol['end']))
        return len(ol_bases & cmp_bases) / float(len(cmp_bases))

def _combine_peaks(peak_files, oligos, bam_files, peak_screener,
                   overlap_percent):
    """Retrieve only peaks that overlap in all three experiments.
    """
    itrees = [_peaks_to_intervaltree(f) for f in peak_files]
    for chrom, base_itree in itrees[0].iteritems():
        overlapper = PeakOverlapper([t[chrom] for t in itrees[1:]], chrom,
                                    overlap_percent, peak_screener, oligos[0])
        base_itree.traverse(overlapper.check_region)
        for peak in overlapper.final:
            yield peak

def output_combine_peaks(peak_files, oligos, bam_files, peak_screener,
                         overlap_percent, config):
    with peak_writer_prep(peak_files[0], config, "combined") as write_peak:
        if write_peak:
            for peak in _combine_peaks(peak_files, oligos, bam_files,
                                       peak_screener, overlap_percent):
                write_peak(peak)

def filter_peaks(peak_file, oligos, bam_file, peak_screener, config):
    with peak_writer_prep(peak_file, config, "filter") as write_peak:
        if write_peak:
            peak_screener.cur_bam = bam_file
            for peak in _read_macs_out(peak_file):
                screen_peak = peak_screener.check(peak, oligos)
                if screen_peak:
                    write_peak(screen_peak)
            peak_screener.cur_bam = None

@contextlib.contextmanager
def peak_writer_prep(peak_file, config, suffix):
    """Handle writing out peak information in multiple formats.

    Supports 'bed' and 'detailed'
    'detailed' is an extended tab delimited format with additional stats.
    """
    out_type = config["algorithm"].get("output", "bed")
    if out_type == "bed":
        header = None
        ext = "bed"
    elif out_type == "detailed":
        header = ["chrom", "start", "end", "width", "count", "blastoligo", "readcount"]
        ext = "tsv"
    else:
        raise NotImplementedError(out_type)
    out_file = "%s-%s.%s" % (os.path.splitext(peak_file)[0], suffix, ext)
    if not os.path.exists(out_file):
        out_handle = open(out_file, "w")
        writer = csv.writer(out_handle, dialect="excel-tab")
        if header:
            writer.writerow(header)
        def _write_peak(peak):
            if out_type == "bed":
                writer.writerow([peak["chrom"], peak["start"], peak["end"]])
            elif out_type == "detailed":
                writer.writerow([peak.get(h, "") for h in header])
        yield _write_peak
        out_handle.close()
    else:
        yield None

def _read_macs_out(peak_file):
    if peak_file.endswith("_peaks.xls"):
        fn = _read_macs_xls
    elif peak_file.endswith("subpeaks.bed"):
        fn = _read_macs_subpeaks
    elif peak_file.endswith(".bed"):
        fn = _read_bed_peaks
    else:
        raise ValueError("Unexpected file %s" % peak_file)
    with open(peak_file) as in_handle:
        for p in fn(in_handle):
            yield p

def _read_bed_peaks(in_handle):
    reader = csv.reader(in_handle, dialect="excel-tab")
    for chrom, start, end, _, score in reader:
        width = int(end) - int(start)
        yield dict(chrom=chrom, start=max(int(start), 1), end=int(end), width=width,
                   score=float(score))

def _read_macs_subpeaks(in_handle):
    """Read peak details from a subpeaks MACS file.
    """
    reader = csv.reader(in_handle, dialect="excel-tab")
    reader.next() # header
    for chrom, start, end, height, _ in reader:
        width = int(end) - int(start)
        yield dict(chrom=chrom, start=max(int(start), 1), end=int(end), width=width,
                   height=int(height))

def _read_macs_xls(in_handle):
    """Read information from a MACS xls output file.
    """
    # read off header
    for line in in_handle:
        if line.strip() and not line.startswith("#"):
            break
    # read the file
    reader = csv.reader(in_handle, dialect="excel-tab")
    for chrom, start, end, width, _, count, score in (l[:7] for l in reader):
        yield dict(chrom=chrom, start=max(int(start), 1), end=int(end), width=int(width),
                   count=int(count), score=float(score))

def _read_bed(in_handle):
    reader = csv.reader(in_handle, dialect="excel-tab")
    for chrom, start, end in ((l[0], int(l[1]), int(l[2])) for l in reader):
        yield dict(chrom=chrom, start=max(start, 1), end=end)

def _peaks_to_intervaltree(peak_file):
    itree = collections.defaultdict(IntervalTree)
    for peak in _read_macs_out(peak_file):
        itree[peak["chrom"]].insert(peak["start"], peak["end"], peak)
    return itree

# ## Remove peak regions where DNA sequence matches binding oligos
# A potential source of pull down noise are DNA segments that overlap with
# the oligos used in the pull down experiment. This detects and removes these
# regions to yield a set of clean peaks for overlap experiments.

class PeakCrossReactScreener:
    def __init__(self, ref, xreact_bases, min_width, min_height,
                 min_count, min_score, min_reads, quick_screen,
                 tmp_dir):
        self._2bit = self._get_twobit(ref)
        self._xreact = xreact_bases
        self._min_width = min_width
        self._min_height = min_height
        self._min_count = min_count
        self._min_score = min_score
        self._min_reads = min_reads
        self._quick = quick_screen
        self._tmpdir = tmp_dir
        self.cur_bam = None

        print "Screening with %s oligo, %s width, %s height, %s count; %s score; quick %s" % \
              (xreact_bases, min_width, min_height, min_count, min_score, quick_screen)

    def _get_twobit(self, ref):
        """Retrieve bx-python 2bit class from /reference hierarchy.
        """
        fname = os.path.join(os.path.dirname(ref), os.pardir, "ucsc",
                             "%s.2bit" % os.path.basename(ref))
        return TwoBitFile(open(fname))

    def check(self, peak, oligos):
        peak["blastoligo"] = 0
        peak["readcount"] = sys.maxint
        # shortcut the BLAST step if we can screen without it
        if self._quick and self._screen_this(peak):
            return None
        # check for score and end if we screen by this
        if self.cur_bam:
            peak["readcount"] = self.cur_bam.read_count(peak["chrom"], peak["start"],
                                                        peak["end"])
        if self._quick and self._screen_this(peak):
            return None
        # check for regions containing oligo matches
        seq = self._2bit[peak["chrom"]].get(peak["start"], peak["end"])
        blast_db = self._make_blastdb(seq)
        input_file = self._make_input_file(oligos)
        blast_out = os.path.join(self._tmpdir, "out.blast")
        cl = NcbiblastnCommandline(query=input_file, db=blast_db,
                                   out=blast_out, outfmt=5, num_alignments=1,
                                   num_descriptions=0, task="blastn-short")
        subprocess.check_call(str(cl).split())
        peak["blastoligo"] = self._max_identity(blast_out)
        if self._screen_this(peak):
            return None
        else:
            return peak

    def _max_identity(self, blast_out):
        identities = [0]
        with open(blast_out) as blast_handle:
            for rec in NCBIXML.parse(blast_handle):
                for align in rec.alignments:
                    for hsp in align.hsps:
                        identities.append(hsp.identities)
        return max(identities)

    def _screen_this(self, peak):
        """Determine if we should screen the current peak based on settings.
        """
        if self._min_width > 0 and peak["end"] - peak["start"] < self._min_width:
            return True
        if self._min_height > 0 and peak["height"] < self._min_height:
            return True
        if self._min_count > 0 and peak["count"] < self._min_count:
            return True
        if self._min_score > 0 and peak["score"] < self._min_score:
            return True
        if self._xreact > 0 and peak["blastoligo"] >= self._xreact:
            return True
        if self._min_reads > 0 and peak["readcount"] < self._min_reads:
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

def call_bayespeak_peaks(exp_bed, back_bed, out_base, ref, config):
    """Peak calling using BayesPeak.
    """
    peak_dir = config["peak_dir"]
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    out_file = os.path.join(peak_dir, "%s-%s-bayespeak.bed" %
                            (out_base.replace(" ", "_"), os.path.basename(ref)))
    if not os.path.exists(out_file):
        rpy.r.assign("treat.bed", exp_bed)
        rpy.r.assign("back.bed", back_bed)
        rpy.r.assign("out.bed", out_file)
        rpy.r('''
        library(multicore)
        library(BayesPeak)
        library(rtracklayer)

        raw.output <- bayespeak(treat.bed, back.bed,
                                use.multicore=TRUE, mc.cores=4)
        output <- summarize.peaks(raw.output, method="lowerbound")
        output$score <- output$PP
        export(output, out.bed)
        ''')
    return out_file

def call_macs_peaks(exp_bed, back_bed, out_base, ref, config):
    """Peak calling using MACS.
    """
    peak_dir = config["peak_dir"]
    subpeaks = config["algorithm"].get("subpeaks", False)
    shiftsize = config["algorithm"].get("shiftsize", None)
    largelambda = config["algorithm"].get("largelambda", None)
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    genome_size = _size_from_fai(_get_ref_fai(ref))
    out_name = os.path.join(peak_dir, "%s-%s-macs" %
                            (out_base.replace(" ", "_"), os.path.basename(ref)))
    ext = "peaks.subpeaks.bed" if subpeaks else "peaks.xls"
    peak_file = "%s_%s" % (out_name, ext)
    if not os.path.exists(peak_file):
        cl = ["macs14", "-t", _full_path(exp_bed), "-c", _full_path(back_bed),
              "--name=%s" % os.path.basename(out_name),
              "--format=BED", "--gsize=%s" % genome_size]
        if largelambda:
            cl += ["--llocal", str(largelambda)]
        if shiftsize:
            cl += ["--nomodel", "--shiftsize=%s" % shiftsize]
        if subpeaks:
            cl += ["--call-subpeaks", "--wig"]
        with chdir(os.path.dirname(out_name)):
            print " ".join(cl)
            subprocess.check_call(cl)
    return peak_file

def call_pics_peaks(exp_bed, back_bed, out_base, ref, config):
    """Peak calling with PICS, using external R script.
    """
    # get mappability files
    mfiles = glob.glob(os.path.join(os.path.dirname(ref), os.pardir, "mappability",
                                    "%s*bed" % os.path.basename(ref)))
    assert len(mfiles) == 1, "Could not find mappability files"
    mfile = mfiles[0]
    pics_script = os.path.join(os.path.dirname(__file__), "chip_seq_w_pics.R")
    peak_dir = config["peak_dir"]
    if not os.path.exists(peak_dir):
        os.makedirs(peak_dir)
    out_file = os.path.join(peak_dir, "%s-%s-PICS.bed" % (out_base.replace(" ", "_"),
                                                          os.path.basename(ref)))
    if not os.path.exists(out_file):
        cl = ["Rscript", pics_script, exp_bed, back_bed, mfile, out_file]
        subprocess.check_call(cl)
    return out_file

def _full_path(fname):
    if not fname.startswith("/"):
        fname = os.path.join(os.getcwd(), fname)
    return fname

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
    fastq_dir = cur_fastq.get("fastq_dir", None)
    if not fastq_dir:
        fastq_dir = config["fastq_dir"]
    fnames = [os.path.join(fastq_dir, f) for f in cur_fastq["files"]]
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
