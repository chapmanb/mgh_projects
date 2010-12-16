"""Examine non-matching reads to determine causes of alignment failures.

This looks at end alignment stats to determine if there is a possible pattern
of end adaptors preventing alignment or cross-contamination:

- Align 5' and 3' non-matching ends to reference genome.
- Align to potential contaminating sources (ribosomal, E coli, yeast)

Usage:
    examine_nomatch.py <config_file>
"""
import sys
import os
import subprocess

import yaml

from bcbio import picard

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    picard_run = picard.PicardRunner(config["programs"]["picard"])
    fastq_in = convert_to_fastq(config["files"]["nomatch"], picard_run)
    ref_genome = config["files"]["reference"]
    trim_size = str(fastq_half_size(fastq_in))
    run_bowtie(fastq_in, ref_genome, "3trim", "0", 1e6,
            ["-3", trim_size])
    run_bowtie(fastq_in, ref_genome, "5trim", "0", 1e6,
            ["-5", trim_size])

    for contam_ref in config["files"]["contam_refs"]:
        name = os.path.basename(contam_ref)
        run_bowtie(fastq_in, contam_ref, name, limit=1e6)

def run_bowtie(fastq_in, ref_genome, file_ext, errors = "2",
        limit = None, extra_params = None):
    out_file = "%s-%s.sam" % (os.path.splitext(fastq_in)[0], file_ext)
    if not os.path.exists(out_file):
        cl = ["bowtie", "--solexa1.3-quals", "--best",
              "-S", "-M 1", "-v", errors, ref_genome, fastq_in, out_file]
        if limit:
            cl += ["-u", str(int(limit))]
        if extra_params:
            cl += extra_params
        print " ".join(cl)
        subprocess.check_call(cl)

def fastq_half_size(in_file):
    with open(in_file) as in_handle:
        name = in_handle.readline()
        seq = in_handle.readline()
    return len(seq.rstrip()) // 2

def convert_to_fastq(in_file, picard_run):
    out_file = "%s.fastq" % os.path.splitext(in_file)[0]
    if not os.path.exists(out_file):
        opts = [("INPUT", in_file),
                ("FASTQ", out_file)]
        picard_run.run("SamToFastq", opts)
    return out_file

if __name__ == "__main__":
    main(*sys.argv[1:])
