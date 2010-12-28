#!/usr/bin/env RScript
#
# Assess chromosome bias in a BAM file.
#
# Based on Jeremy's writeup and discussion:
# http://jermdemo.blogspot.com/2010/12/chromosome-bias-in-r-my-notebook.html
#
# Usage:
#   chr_bias.R <input file> <genome name> <build name>
#     <genome name> is like Mmusculus, Hsapiens
#     <build name> is like mm8, hg19

library(BSgenome)
library(IRanges)
library(GenomicRanges)
library(Rsamtools)
library(ggplot2)

options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
names(args) <- c("infile", "genome", "build")

library(paste("BSgenome", args["genome"], "UCSC", args["build"], sep="."),
        character.only=TRUE)
min.records <- 5000
genome <- get(args["genome"])
plot.file <- gsub('.bam', '-sizes.pdf', args["infile"])
reg.plot.file <- gsub('.bam', '-reg.pdf', args["infile"])

# Generate data frame of masked chromosome sizes
unmaskedWidth <- function(chr) {length(chr) - maskedwidth(chr)}
dfSizes <- function(chr) {
  data.frame(space=chr, size=length(genome[[chr]]),
             unmaskedsize=unmaskedWidth(genome[[chr]]))
}
maskedSizes <- ldply(.data=seqnames(genome), .fun=dfSizes,
                     .progress="text", .parallel=TRUE)

# Retrieve count of reads per chromosome for the input BAM file
seqlengths2gr <- function(x, strand="*") {
  GRanges(names(x), IRanges(1, x), strand=strand)
}
gr <- seqlengths2gr(seqlengths(genome))
bamCnt <- countBam(args["infile"], param=ScanBamParam(which=gr))

# Merge, subset based on a minimum number of mapped records, and plot
allCombined <- merge(maskedSizes, bamCnt, sort=FALSE)
chrSizesReads <- subset(allCombined, records > min.records)
print(chrSizesReads)
p <- ggplot(data=chrSizesReads, aes(x=unmaskedsize, y=records, label=space)) +
     geom_point() + geom_text(vjust=2, size=3) +
     stat_smooth(method="lm", se=TRUE, level=0.95) +
     ylab("Reads aligned") +
     xlab("Unmasked chromosome size") +
     opts(title = "Reads vs Chromosome Size")
ggsave(plot.file, p, width=6, height=6)

# Check for outliers with linear regression
size.lm <- lm(records~unmaskedsize, data=chrSizesReads)
print(summary(size.lm))
p <- qplot(chrSizesReads$space, rstandard(size.lm)) +
     aes(label=chrSizesReads$space) +
     geom_text(vjust=2, size=3) +
     xlab("Chromosome") +
     ylab("Std Residual from lm (reads)") +
     geom_abline(slope=0, intercept=0) +
     opts(axis.text.x = theme_text(angle=45, hjust=1)) +
     opts(title = "Linear Regression Residuals")
ggsave(reg.plot.file, p, width=6, height=6)
