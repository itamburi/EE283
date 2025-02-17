# load the library
# it seems like some needed libraries are not loaded automatically
library(ATACseqQC)
library(Rsamtools)
library(cowplot)
library(ggplot2)
# define the bamfile in the working directory
# this code will run MUCH faster if the bam for testing is a single chromosome
getwd()
setwd("/Users/ian/Library/Application Support/CRSP Desktop/Volumes.noindex/pub.localized/ee283/week6")

## 0.0 - Read in Bam files ####

# QCing 6 files
bamfiles = list.files("bam.processed")
bamfiles = subset(bamfiles, grepl(".chrX.dd.bam$",bamfiles))
bampaths = file.path("bam.processed", bamfiles)

# Apply indexBam() to each file if we didnt already
# lapply(bampaths, indexBam)


## 1.0 - Check Frag size and Lib Complexity ####

mkqcplots = function(bampath){

  basename(bampath)
  # generate fragement size distribution
  pdf(paste0("qcplots/FragSizeDist_",basename(bampath),".pdf"))
  fragSize = fragSizeDist(bampath, basename(bampath))
  dev.off()
  
  ## estimate complexity
  pdf(paste0("qcplots/LibComplexity_",basename(bampath),".pdf"))
  estimateLibComplexity(readsDupFreq(bampath, index=bampath))
  dev.off()
  
}
  
lapply(bampaths, mkqcplots)




## 2.0 -Shift G alignments and split bam by fragment size (nucleosome free, mono, di-, etc) ####
# shift the GAlignmentsLists by 5' ends.
# All reads aligning to the positive strand will be offset by +4bp,
# and all reads aligning to the negative strand will be offset -5bp

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txs = transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
TSS = promoters(txs, upstream=0, downstream=1)
TSS = unique(TSS)

shiftandsplit = function(bampath){
  
  #bampath = bampaths[1] #debug
  dir = paste0("qcfiles/", gsub(".chrX.dd.bam","",basename(bampath)))
  dir.create(dir)
  
  gal = readBamFile(bampath, asMates=TRUE, bigFile=TRUE)
  gal1 = shiftGAlignmentsList(gal)
  
  dir2=paste0(dir,"/",basename(dir),".shifted.bam")
  BiocIO::export(gal1, dir2)
  
  objs = splitGAlignmentsByCut(gal1, txs = txs, outPath = dir)
  null = writeListOfGAlignments(objs, dir)
  
}
lapply(bampaths, shiftandsplit)




## 4.0 -calc near transcr start sites, and Feature aligned heatmap of log2 transformed signals ####

# estimate the library size for normalization
library(ChIPpeakAnno)

makeqcheatmap = function(bf){
  
  #bf=bamfiles[2] # debug
  
  # estimate lib size
  file = paste0("bam.processed/",bf)
  librarySize = estLibSize(file)
  
  # Only consider X chromosome for expediency
  # check chrom annotations for TSS and bamfile
  head(seqlevels(BamFile(file)))
  head(seqlevels(TSS))
  seqlevels(TSS) <- sub("^chr", "", seqlevels(TSS))
  head(seqlevels(TSS))
  seqlev = "X"
  # seqlev = c("chrX","chr2L","chr2R","chr3L","chr3R")
  NTILE = 101
  dws = ups = 1010
  
  # make paths to qc files for enrichedFragments()
  dir=paste0("qcfiles/", gsub(".chrX.dd.bam","",bf))
  qcfiles = file.path(dir,c("NucleosomeFree.bam","mononucleosome.bam","dinucleosome.bam","trinucleosome.bam"))
  
  sigs = enrichedFragments(qcfiles, TSS=TSS, librarySize=librarySize, seqlev=seqlev,
                           TSS.filter=0.5, n.tile = NTILE, upstream = ups, downstream = dws)
  
  names(sigs) <- gsub(".bam", "", basename(names(sigs)))
  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
  #plot heatmap
  pdf(paste0("qcplots/FeatureAlignedHeatmap_", bf,".pdf"))
  featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws), zeroAt=.5, n.tile=NTILE)
  dev.off()
  
  
}

lapply(bamfiles, makeqcheatmap)



# It is possible that reference genomes for many species in SCGB are
# missing or that those libraries are not up2date in terms of releases.
# I am not sure if the package includes instructions for rolling your own files
# given assemblies and a GTF file (these two files should be sufficient)? It does seem
# possible to define the TSS for a genome given a GTF file.

# library(GenomicFeatures)
# txdbmm = makeTxDbFromGFF(file = "ref/dm6.ncbiRefSeq.gtf")
# txs2 = transcripts(txdbmm)
# TSS2 = promoters(txs2, upstream=0, downstream=1)
# TSS2 = unique(TSS) 





