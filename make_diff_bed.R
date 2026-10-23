#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop("Usage: make_diff_bed.R diff.csv consensus.bed output.bed")
d <- read.csv(args[1], stringsAsFactors=FALSE)
p <- read.delim(args[2], header=FALSE, stringsAsFactors=FALSE)
colnames(p)[1:4] <- c("chr","start","end","peak_id")
m <- merge(d,p[,1:4],by="peak_id")
write.table(m[,c("chr","start","end","peak_id")], args[3], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
