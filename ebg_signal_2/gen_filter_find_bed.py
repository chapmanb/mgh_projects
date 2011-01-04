#!/usr/bin/env python
"""Generate ready to query BED files of regions to filter and analyze.

Usage:
    gen_filter_find_bed.py <input config>

Where the configuration file is a YAML formatted file specifying each of
the files to be processed.
"""
import sys
import csv
import os
import pprint

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    write_bed(config["filter_out"],
              filter_regions(config["gene_gff"]))
    write_bed(config["find_out"],
              target_regions(config["gene_gff"], config["target_csv"]),
              shortrna_regions(config["mirna_gff"], config["star_csv"],
                  config["seq_genome"]))

def write_bed(out_file, *region_generators):
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, dialect="excel-tab")
            for region_gen in region_generators:
                for info in region_gen:
                    writer.writerow(info)

def target_regions(gene_gff, target_csv):
    """Read gene coordinates for target genes of interest.
    """
    targets = dict()
    with open(target_csv) as in_handle:
        reader = csv.reader(in_handle)
        for gene, full_name in reader:
            targets[gene.upper()] = full_name.strip()

    interest_types = ["gene", "mRNA", "ncRNA", "pseudogene",
            "pseudogenic_transcript"]
    limit = dict(gff_type=interest_types)
    for rec in GFF.parse(gene_gff, limit_info=limit):
        chrom = _arab_name_fix(rec.id)
        for g_feat in rec.features:
            name = g_feat.qualifiers["ID"][0]
            if targets.has_key(name):
                tx_feat = g_feat.sub_features[-1]
                yield (chrom, tx_feat.location.nofuzzy_start,
                        tx_feat.location.nofuzzy_end, targets[name])
                del targets[name]
    assert len(targets) == 0, targets

def shortrna_regions(mirna_gff, star_csv, seq_file):
    """Return miRNA sequences with corresponding guide and star regions.
    """
    seq_index = SeqIO.index(seq_file, "fasta")
    mirna_seqs = dict()
    with open(star_csv) as in_handle:
        for name, guide, star in csv.reader(in_handle):
            mirna_seqs[name] = (guide.strip(), star.strip())

    for rec in GFF.parse(mirna_gff):
        cur_seq = str(seq_index[rec.id].seq)
        for f in rec.features:
            name = f.qualifiers["ID"][0]
            start, end = (f.location.nofuzzy_start, f.location.nofuzzy_end)
            yield (rec.id, start, end, name)
            #guide, star = mirna_seqs.get(name, ("", ""))
            for seq_name, guide, star in [(n, g, s) for n, (g, s) in
                    mirna_seqs.iteritems() if n.startswith(name)]:
                for find_seq, ext in [(guide, "guide"), (star, "star")]:
                    if find_seq:
                        if f.strand == -1:
                            find_seq = str(Seq(find_seq).reverse_complement())
                        region = cur_seq[start:end]
                        pos = region.find(find_seq)
                        if pos > -1:
                            yield (rec.id, start + pos, start + pos + len(find_seq),
                                    "%s_%s" % (seq_name, ext))
                        else:
                            print f.strand, name, ext, pos, find_seq, region
                            raise NotImplementedError

def filter_regions(gff_file):
    """Generate (chrom, start, end, name) of RNA regions to filter.
    """
    base_gene = ["transposable_element_gene", "gene"]
    to_use = ["rRNA", "tRNA", "snoRNA", "snRNA"]
    extra_rnas = ["mRNA", "miRNA", "ncRNA"]
    sub_features = ["exon"]
    limit = dict(gff_type = base_gene + to_use + extra_rnas + sub_features)

    for rec in GFF.parse(gff_file, limit_info=limit, target_lines=100):
        chrom = _arab_name_fix(rec.id)
        for g_feat in [f for f in rec.features if f.type.find("gene") >= 0]:
            for tx_feat in [f for f in g_feat.sub_features if f.type in to_use]:
                name = tx_feat.qualifiers["ID"][0]
                for i, cds_feat in enumerate(tx_feat.sub_features):
                    yield (chrom, cds_feat.location.nofuzzy_start,
                            cds_feat.location.nofuzzy_end, "%s_%s" % (name, i))

def _arab_name_fix(base):
    if base == "ChrC":
        return "Pt"
    elif base == "ChrM":
        return "Mt"
    else:
        return base.replace("Chr", "")


if __name__ == "__main__":
    main(*sys.argv[1:])

