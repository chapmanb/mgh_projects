#Perform ChIP-seq peak calling using the Bioconductor PICS package.
#
# http://wiki.rglab.org/index.php?title=Public:PICS
#
# Usage:
#    Rscript chip_seq_w_pics.R <treatment BED> <background BED> <mappability BED>
#                              [<out BED file>]

library(multicore)
library(PICS)
library(rtracklayer)

options(cores=5)
TESTING_CHR <- "chrX"
INTERACTIVE <- FALSE
MINFRAGSIZE <- 50
MAXFRAGSIZE <- 5000

args <- (commandArgs(TRUE))
print(args)
treatFile <- args[1]
backFile <- args[2]
mapFile <- args[3]
outFile <- args[4]
if (is.na(outFile)) {
  sop <- match(".", rev(strsplit(treatFile, NULL)[[1]]))[1]
  baseName <- substring(treatFile, 1, nchar(treatFile) - sop)
  outFile <- paste(baseName, "-PICS.bed", sep="")
}
print(outFile)

# Convert a full BED file with strand into a GenomeData object
fileToGenomeData <- function(inFile) {
  print(inFile)
  data <- read.delim(inFile, header=FALSE,
  		   colClass=c("character", "integer", "integer", "NULL",
  			      "NULL", "character"))
  names(data)[1] <- "space"
  names(data)[2] <- "start"
  names(data)[3] <- "end"
  names(data)[4] <- "strand"
  if (!is.null(TESTING_CHR)) {
    data <- subset(data, space==TESTING_CHR)
  }
  data <- as(data, "GenomeData")
  data <- unique(data)
}

# Convert an input simple BED file into a ranged data object
fileToRangedData <- function(inFile) {
  data <- read.delim(inFile, header=FALSE,
                     colClass=c("character", "integer", "integer", "NULL"))
  names(data)[1] <- "space"
  names(data)[2] <- "start"
  names(data)[3] <- "end"
  data <- as(data, "RangedData")
}
# load input files
print("Import BED file data")
treatData <- fileToGenomeData(treatFile)
backData <- fileToGenomeData(backFile)
mapData <- fileToRangedData(mapFile)
print("Run PICS")
# run PICS on the treatment and control sets
seg <- segmentReads(treatData, backData, mapData, minReads=NULL)
pics <- PICS(seg, dataType="TF")
segC <- segmentReads(backData, treatData, mapData, minReads=NULL)
picsC <- PICS(segC, dataType="TF")

print("Select FDR and filter")

# Select a FDR of 10% and use the score for output filtering
stdFilter <- list(se=c(0, 50), sigmaSqF=c(0, 22500), sigmaSqR=c(0, 22500))
fdrFilter <- c(stdFilter, list(delta=c(MINFRAGSIZE, Inf)))
fdr <- picsFDR(pics, picsC, filter=fdrFilter)
fdrFilter <- subset(fdr, FDR < 0.05)
minScore <- fdrFilter$score[1]
print(fdr)
print(minScore)

if (INTERACTIVE) {
  plot(pics, picsC, xlim=c(2, 8), ylim=c(0, 0.2),
       filter=c(stdFilter, list(delta=c(MINFRAGSIZE, MAXFRAGSIZE))))
  readline("Press <Enter> to continue")
}

picsFilter <- c(stdFilter, list(score=c(minScore, Inf), 
                                delta=c(MINFRAGSIZE, MAXFRAGSIZE)))
outBed <- makeRangedDataOutput(pics, type="bed", filter=picsFilter)
export(outBed, outFile)
