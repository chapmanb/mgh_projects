#!/usr/bin/env python
"""Combine trimmed and non-trimmed reads, providing alignments and stats.

Usage:
  cg_combine_map.py <YAML config>
"""
import os
import sys
import glob
import shutil
import collections
import subprocess

import yaml
import pysam
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

from BCBio import GFF
from bcbio.picard import run, PicardRunner

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    picard = PicardRunner(config["programs"]["picard"])
    errors = config["algorithm"].get("errors", None)
    outdir = config["outdir"]
    scriptdir = os.path.dirname(__file__)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    map_trim = trim_three_end(config["fastq"]["align"], 2)
    full_fastq = join_map_and_trim(map_trim, config["fastq"]["trim"])
    if errors is None:
        test_errors(full_fastq, config["reference"])
    else:
        do_aligns(full_fastq, config["reference"], errors, outdir, picard,
                  scriptdir)
    final_dir = os.path.join(outdir, os.pardir, "final")
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
    prepare_final(full_fastq, config["reference"], outdir, final_dir)

def trim_three_end(in_fastq, num_bases):
    out_fastq = "%s-trimend%s" % os.path.splitext(in_fastq)
    if not os.path.exists(out_fastq):
        three_end_distribution(in_fastq, num_bases, 10)
        with open(in_fastq) as in_handle:
            with open(out_fastq, "w") as out_handle:
                for (name, seq, qual) in FastqGeneralIterator(in_handle):
                    out_handle.write("@%s\n%s\n+\n%s\n" % (name, seq[:-num_bases],
                                                           qual[:-num_bases]))
    return out_fastq

def three_end_distribution(in_fastq, num_bases, to_display):
    """Print 3' end distribution of reads to check for small adapter.
    """
    counts = collections.defaultdict(int)
    with open(in_fastq) as in_handle:
        for (_, seq, _) in FastqGeneralIterator(in_handle):
            counts[seq[-num_bases:]] += 1
    counts = [(v, k) for (k, v) in counts.iteritems()]
    counts.sort(reverse=True)
    for count, seq in counts[:to_display]:
        print seq, count

def join_map_and_trim(mapped, trimmed):
    """Generate combined fastq file of mapped and trimmed reads.
    """
    out_file = "%sfinal.fastq" % os.path.commonprefix([mapped, trimmed])
    if not os.path.exists(out_file):
        with open(out_file, "w") as out_handle:
            for fname in [mapped, trimmed]:
                with open(fname) as in_handle:
                    for line in in_handle:
                        out_handle.write(line)
    return out_file

def test_errors(in_fastq, references):
    """Examine align rates with different error choices.
    """
    error_checks = [0, 1, 2, 3]
    for error in error_checks:
        for ref in references:
            print ref["file"], error
            run_bowtie(in_fastq, ref["file"], None, error, 1e6)

def to_strand_bam(in_file, forward=True):
    """Convert to BAM file with either forward or reverse strand.
    """
    fext = "for" if forward else "rev"
    out_file = "%s-%s.bam" % (os.path.splitext(in_file)[0], fext)
    if not os.path.exists(out_file):
        in_sam = pysam.Samfile(in_file, "r")
        out_bam = pysam.Samfile(out_file, "wb", template=in_sam)
        for read in in_sam:
            if not read.is_unmapped:
                if ((forward and not read.is_reverse) or
                    (not forward and read.is_reverse)):
                    out_bam.write(read)
        in_sam.close()
        out_bam.close()
    return out_file

def bam_to_wig(in_bam):
    out_file = "%s.bigwig" % os.path.splitext(in_bam)[0]
    if not os.path.exists(out_file):
        cl = ["bam_to_wiggle.py", in_bam, "--outfile=%s" % out_file]
        subprocess.check_call(cl)
    return out_file

def plot_chr_dist(bam_files, ref, picard, scriptdir):
    """Calculate stats on chromosome distribution of reads.
    """
    combine_bam = run.picard_merge(picard, bam_files)
    run.picard_index(picard, combine_bam)
    cl = ["Rscript", os.path.join(scriptdir, "chr_bias.R"), combine_bam,
          ref["genome"], ref["build"]]
    subprocess.check_call(cl)

def _filter_features(in_recs, max_size):
    """Remove large features from record which hide sub features in IGV.
    """
    for rec in in_recs:
        final = []
        for f in rec.features:
            if len(f.location) < max_size:
                final.append(f)
            else:
                for sub in f.sub_features:
                    if len(sub.location) < max_size:
                        final.append(sub)
        rec.annotations = {}
        rec.features = final
        yield rec

def genbank_to_gff(gb_file):
    """Convert GenBank file to GFF for IGV display.
    """
    max_size = 1e4
    gff_file = "%s.gff3" % os.path.splitext(gb_file)[0]
    if not os.path.exists(gff_file):
        with open(gb_file) as in_handle:
            with open(gff_file, "w") as out_handle:
                gb_iterator = SeqIO.parse(in_handle, "genbank")
                GFF.write(_filter_features(gb_iterator, max_size),
                          out_handle)

def do_aligns(in_fastq, references, errors, outdir, picard, scriptdir):
    """Perform final alignments, preparing sense/anti BAM and BigWig files.
    """
    for ref in references:
        out_sam = run_bowtie(in_fastq, ref["file"], os.path.basename(ref["file"]),
                             errors, out_dir=outdir)
        cur_bams = []
        for strand in [True, False]:
            out_bam = to_strand_bam(out_sam, strand)
            sort_bam = run.picard_sort(picard, out_bam)
            cur_bams.append(sort_bam)
            bam_to_wig(sort_bam)
        if ref.get("chr_dist", False):
            #plot_chr_dist(cur_bams, ref, picard, scriptdir)
            pass
        if ref.get("feature_prep", False):
            genbank_to_gff("%s.gb" % ref["file"])

def run_bowtie(fastq_in, ref_genome, file_ext, errors="2",
        limit=None, extra_params=None, out_dir=None):
    if file_ext:
        out_file = "%s-%s.sam" % (os.path.splitext(fastq_in)[0], file_ext)
        if out_dir:
            out_file = os.path.join(out_dir, os.path.basename(out_file))
    else:
        out_file = "/dev/null"
    if file_ext is None or not os.path.exists(out_file):
        cl = ["bowtie", "--best", "-S", "-m 1", "-v", str(errors), ref_genome,
              fastq_in, out_file]
        if limit:
            cl += ["-u", str(int(limit))]
        if extra_params:
            cl += extra_params
        print " ".join(cl)
        subprocess.check_call(cl)
    return out_file

def _clean_filename(fname):
    to_rem = ["-final", "-sort", "-merge"]
    for rem in to_rem:
        fname = fname.replace(rem, "")
    return fname

def _make_outdir(outdir_base, base_ref):
    outdir = os.path.join(outdir_base, base_ref)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

def _copy_galaxy_files(base, base_ref, outdir_base):
    outdir = _make_outdir(outdir_base, base_ref)
    exts = ["-sort.bam", "-sort.bigwig", ".pdf"]
    for ext in exts:
        for initial in glob.glob("%s-%s*%s" % (base, base_ref, ext)):
            final = os.path.join(outdir, _clean_filename(os.path.basename(initial)))
            shutil.copy(initial, final)

def _copy_igv_files(base, base_ref, ref_file, outdir_base):
    outdir = _make_outdir(outdir_base, base_ref)
    exts = ["-sort.bam", "-sort.bam.bai"]
    for ext in exts:
        for initial in glob.glob("%s-%s*%s" % (base, base_ref, ext)):
            final = os.path.join(outdir, _clean_filename(os.path.basename(initial)))
            shutil.copy(initial, final)
    rexts = [".fasta", ".gff3"]
    for ext in rexts:
        initial = "%s%s" % (ref_file, ext)
        final = os.path.join(outdir, os.path.basename(initial))
        shutil.copy(initial, final)

def prepare_final(fastq_in, refs, aligndir, outdir):
    """Organize output files for re-distribution.
    """
    base = os.path.join(aligndir, os.path.splitext(os.path.basename(fastq_in))[0])
    for ref in refs:
        base_ref = os.path.basename(ref["file"])
        if ref.get("chr_dist", False):
            _copy_galaxy_files(base, base_ref, outdir)
        if ref.get("feature_prep", False):
            _copy_igv_files(base, base_ref, ref["file"], outdir)

if __name__ == "__main__":
    main(*sys.argv[1:])
