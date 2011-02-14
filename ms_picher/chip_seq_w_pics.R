#Perform ChIP-seq peak calling using the Bioconductor PICS package.
#
# http://wiki.rglab.org/index.php?title=Public:PICS
#
# Usage:
#    Rscript chip_seq_w_pics.R <treatment BED> <background BED> <mappability BED>

library(multicore)
library(PICS)
library(rtracklayer)

options(cores=8)
FRAGSIZE <- 300
MINFRAGSIZE <- 150
MAXFRAGSIZE <- 500

args <- (commandArgs(TRUE))
print(args)
treatFile <- args[1]
backFile <- args[2]
mapFile <- args[3]
sop <- match(".", rev(strsplit(treatFile, NULL)[[1]]))[1]
baseName <- substring(treatFile, 1, nchar(treatFile) - sop)
outFile <- paste(baseName, "-PICS.bed", sep="")
print(outFile)

# Convert a full BED file with strand into a GenomeData object
fileToGenomeData <- function(inFile) {
  data <- read.table(inFile, header=FALSE,
  		   colClass=c("character", "integer", "integer", "NULL",
  			      "NULL", "character"))
  names(data)[1] <- "space"
  names(data)[2] <- "start"
  names(data)[3] <- "end"
  names(data)[4] <- "strand"
  data <- data[order(data$space, data$start, data$end),]
  data <- as(data, "RangedData")
  data <- as(data, "GenomeData")
  data <- unique(data)
}

# Convert an input simple BED file into a ranged data object
fileToRangedData <- function(inFile) {
  data <- read.table(inFile, header=TRUE,
                     colClass=c("character", "integer", "integer", "NULL"))
  names(data)[1] <- "space"
  names(data)[2] <- "start"
  names(data)[3] <- "end"
  data <- data[order(data$space, data$start, data$end),]
  data <- as(data, "RangedData")
}
# load input files
treatData <- fileToGenomeData(treatFile)
backData <- fileToGenomeData(backFile)
mapData <- fileToRangedData(mapFile)
# run PICS on the treatment and control sets
setParaPrior(xi=FRAGSIZE,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
seg <- segmentReads(treatData, backData, mapData, minReads=NULL)
pics <- PICS(seg, dataType="TF")
segC <- segmentReads(backData, treatData, mapData, minReads=NULL)
picsC <- PICS(segC, dataType="TF")


# Select a FDR of 10% and use the score for output filtering
stdFilter <- list(se=c(0, 50), sigmaSqF=c(0, 22500), sigmaSqR=c(0, 22500))
fdrFilter = c(stdFilter, list(delta=c(MINFRAGSIZE, Inf)))
fdr <- picsFDR(pics, picsC, filter=fdrFilter)
fdrFilter <- subset(fdr, FDR < 0.10)
minScore <- fdrFilter$score[1]
print(fdr)
print(minScore)

plot(pics, picsC, xlim=c(2, 8), ylim=c(0, 0.2),
     filter=c(stdFilter, list(delta=c(MINFRAGSIZE, MAXFRAGSIZE))))
readline("Press <Enter> to continue")

picsFilter <- c(stdFilter, list(score=c(minScore, Inf), 
                                delta=c(MINFRAGSIZE, MAXFRAGSIZE)))
outBed <- makeRangedDataOutput(pics, type="bed", filter=picsFilter)
export(outBed, outFile)
