#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<4) stop("Usage: peak_annotation.R peaks.bed species assembly out.csv")
peakfile <- args[1]; species <- args[2]; assembly <- args[3]; outfile <- args[4]
suppressPackageStartupMessages({library(ChIPseeker)})
if(species=="human" && assembly=="hg38") {suppressPackageStartupMessages({library(TxDb.Hsapiens.UCSC.hg38.knownGene);library(org.Hs.eg.db)}); txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene; annodb<-"org.Hs.eg.db"} else if(species=="human" && assembly=="hg19") {suppressPackageStartupMessages({library(TxDb.Hsapiens.UCSC.hg19.knownGene);library(org.Hs.eg.db)}); txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene; annodb<-"org.Hs.eg.db"} else if(species=="mouse" && assembly=="mm10") {suppressPackageStartupMessages({library(TxDb.Mmusculus.UCSC.mm10.knownGene);library(org.Mm.eg.db)}); txdb<-TxDb.Mmusculus.UCSC.mm10.knownGene; annodb<-"org.Mm.eg.db"} else stop("Automatic ChIPseeker annotation currently supports hg38/hg19/mm10. For mm39/custom, provide a TxDb-based custom annotation script.")
p <- readPeakFile(peakfile)
a <- annotatePeak(p, tssRegion=c(-3000,3000), TxDb=txdb, annoDb=annodb)
annotated <- as.data.frame(a)
write.csv(annotated, outfile, row.names=FALSE)
original <- read.delim(peakfile, header=FALSE, stringsAsFactors=FALSE)
if(ncol(original)>=4 && "V4" %in% names(annotated)) {
  missing <- original[!original[[4]] %in% annotated$V4,,drop=FALSE]
  write.table(missing,paste0(outfile,".unannotated.bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  if(nrow(missing)>0) warning(nrow(missing)," peaks lacked annotation; see .unannotated.bed (check contig/annotation compatibility)")
}
