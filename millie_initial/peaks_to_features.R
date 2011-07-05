#!/usr/bin/env RScript
#
# Associate peaks from an input BED file with gene features.
# A CSV file is generated of the output table, which contains
# the BED inputs plus an extra feature type column.
#
# Usage:
#   peaks_to_features.R <input file> <output file>

library(IRanges)
library(ChIPpeakAnno)
library(rtracklayer)
library(biomaRt)

options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]

input.rd <- import(infile)
ens.mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Retrieve ranged data from a BioMart instance, saving for backup.
getMartData <- function(mart, ftype) {
  save.file <- paste("annotations/hsapiens-ensembl-", ftype, ".gff3", sep="")
  if (file.exists(save.file)) {
    mart.data <- import(save.file)
  } else {
    if (ftype %in% c("ncRNA", "lincRNA")) {
      mart.data <- getNcData(ftype)
    } else {
      mart.data <- getAnnotation(ens.mart, ftype)
    }
    export(mart.data, save.file)
  }
  mart.data
}

# Retrieve non-coding annotations not covered by ChIPpeakAnno
# http://uswest.ensembl.org/info/docs/genebuild/ncrna.html
getNcData <- function(ftype) {
  mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  other.nc <- c("IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
                "Mt_tRNA", "Mt_tRNA_pseudogene", "TR_V_pseudogene")
  if (ftype == "ncRNA") {
    nc.types <- c("lincRNA", "miRNA_pseudogene", "misc_RNA_pseudogene",
                  "polymorphic_pseudogene", "pseudogene",
                  "snoRNA", "snoRNA_pseudogene", "snRNA_pseudogene",
                  "processed_transcript")
  } else if (ftype == "lincRNA") {
    nc.types <- c("lincRNA")
  }
  attrs <- c("ensembl_gene_id", "chromosome_name", "start_position",
             "end_position", "strand", "description", "gene_biotype")
  result <- getBM(attributes=attrs, filters=c("biotype"), values=nc.types, mart=mart)
  # remove duplicated
  result <- unique(result)
  result <- result[order(result[,3]), ]              
  duplicated.id <- result[duplicated(result[,1]), 1]
  result <- result[!duplicated(result[,1]), ]
  result.rd <- RangedData(IRanges(start=as.numeric(result[,3]), end=as.numeric(result[,4]),
                 names= as.character(result[,1])), strand=result[,5],
                 description=as.character(result[,6]), space=as.character(result[,2]),
	         biotype=as.character(result[,7]))
}

# Associate features with peaks in the input dataframe
# Returns a list of:
#   remain - Remaining peaks which have not been associated
#   associated - Data frame of peaks associated with this ftype
#   annotations - Data frame with details on annotations
annotateWithFeature <- function(mart, ftype, input.df) {
  rownames(input.df) <- seq(1, length(input.df$name))
  mart.data <- getMartData(mart, ftype)
  ann.rd <- annotatePeakInBatch(RangedData(input.df), AnnotationData=mart.data,
                                output="overlapping", multiple=TRUE)
  ann.df <- as.data.frame(ann.rd)

  cur.found <- unique(as.integer(ann.df$peak))
  if (length(cur.found) > 0) {
    input.this <- input.df[cur.found, ]
    if (ftype == "TSS") ftype <- "intron"
    input.this$ftype <- ftype
    print(ftype)
    print(summary(input.this))
  } else {
    input.this <- NULL
  }
  if (length(cur.found) < length(input.df$name)) {
    input.remain <- input.df[!(rownames(input.df) %in% cur.found), ]
  } else {
    input.remain <- NULL
  }
  list(remain=input.remain, associated=input.this, annotations=ann.df)
}

final.table <- NULL
#want <- c("MACS_peak_99175", "MACS_peak_100031")
#cur.table <- as.data.frame(input.rd)
#cur.table <- cur.table[cur.table$name %in% want, ]
#features <- c("ncRNA")
cur.table <- as.data.frame(input.rd)
features <- c("lincRNA", "ncRNA", "miRNA", "Exon", "5utr", "3utr", "TSS")
print(summary(cur.table))
for (ftype in features) {
  annotate <- annotateWithFeature(ens.mart, ftype, cur.table)
  final.table <- rbind(annotate$associated, final.table)
  cur.table <- annotate$remain
}
# add on the remaining peaks which haven't been annotated
if (!is.null(cur.table)) cur.table$ftype <- NA
final.table <- rbind(cur.table, final.table)
print(summary(final.table))
write.csv(final.table, file=outfile)
