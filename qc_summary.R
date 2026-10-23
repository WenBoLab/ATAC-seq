#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<6) stop("Usage: qc_summary.R workdir samples.tsv warn_map warn_dup warn_mt warn_frip")
work <- args[1]; sample_file <- args[2]
warn_map <- as.numeric(args[3]); warn_dup <- as.numeric(args[4]); warn_mt <- as.numeric(args[5]); warn_frip <- as.numeric(args[6])
meta <- read.delim(sample_file, stringsAsFactors=FALSE)
qcdir <- if(length(args)>=7) args[7] else file.path(work,"06.qc"); aligndir <- if(length(args)>=8) args[8] else file.path(work,"04.align")
lib <- read.delim(file.path(qcdir,"library_qc.tsv"), stringsAsFactors=FALSE)
fr <- read.delim(file.path(qcdir,"frip.tsv"), stringsAsFactors=FALSE)

parse_flagstat <- function(f){
  x <- readLines(f, warn=FALSE)
  total <- suppressWarnings(as.numeric(sub(" .*","", x[grepl("in total",x)][1])))
  mapped_line <- x[grepl(" mapped ",x) & grepl("%",x)][1]
  mapped <- suppressWarnings(as.numeric(sub(" .*","",mapped_line)))
  rate <- ifelse(is.finite(total) && total>0, mapped/total, NA_real_)
  c(total=total,mapped=mapped,mapping_rate=rate)
}
parse_dup <- function(f){
  if(!file.exists(f)) return(NA_real_)
  x <- readLines(f, warn=FALSE)
  hdr <- grep("^LIBRARY\t",x)
  if(!length(hdr) || hdr[1]>=length(x)) return(NA_real_)
  h <- strsplit(x[hdr[1]],"\t")[[1]]; v <- strsplit(x[hdr[1]+1],"\t")[[1]]
  i <- match("PERCENT_DUPLICATION",h); if(is.na(i)) return(NA_real_)
  as.numeric(v[i])
}

rows <- lapply(meta$sample, function(s){
  fs <- parse_flagstat(file.path(aligndir,paste0(s,".flagstat.txt")))
  dup <- parse_dup(file.path(qcdir,paste0(s,".dup_metrics.txt")))
  l <- lib[lib$sample==s,,drop=FALSE]; f <- fr[fr$sample==s,,drop=FALSE]
  mt_rate <- if(nrow(l) && l$dedup_records_before_mt[1]>0) l$mt_records_before_removal[1]/l$dedup_records_before_mt[1] else NA_real_
  data.frame(sample=s, group=meta$group[meta$sample==s][1], library_type=l$library_type[1],
             raw_records=fs["total"], mapped_records=fs["mapped"], mapping_rate=fs["mapping_rate"],
             duplication_rate=dup, mitochondrial_fraction=mt_rate,
             final_records=l$final_records[1], peaks=f$peaks[1], FRiP=f$FRiP[1], FRiP_unit=if("unit" %in% names(f)) f$unit[1] else "legacy_read_proxy", stringsAsFactors=FALSE)
})
out <- do.call(rbind,rows)
out$status <- "PASS"
out$warnings <- ""
addwarn <- function(i,msg){ out$status[i] <<- "WARN"; out$warnings[i] <<- paste(c(out$warnings[i],msg)[c(out$warnings[i],msg)!=""],collapse="; ") }
for(i in seq_len(nrow(out))){
  if(any(!is.finite(as.numeric(unlist(out[i,c("mapping_rate","duplication_rate","mitochondrial_fraction","FRiP")]))))) addwarn(i,"missing/nonfinite QC metric")
  if(out$FRiP_unit[i]=="legacy_read_proxy") addwarn(i,"legacy FRiP proxy; recompute fragment FRiP")
  if(!is.na(out$mapping_rate[i]) && out$mapping_rate[i] < warn_map) addwarn(i,paste0("mapping<",warn_map))
  if(!is.na(out$duplication_rate[i]) && out$duplication_rate[i] > warn_dup) addwarn(i,paste0("duplication>",warn_dup))
  if(!is.na(out$mitochondrial_fraction[i]) && out$mitochondrial_fraction[i] > warn_mt) addwarn(i,paste0("mt>",warn_mt))
  if(!is.na(out$FRiP[i]) && out$FRiP[i] < warn_frip) addwarn(i,paste0("FRiP<",warn_frip))
}
write.csv(out,file.path(qcdir,"ATAC_QC_summary.csv"),row.names=FALSE)

pdf(file.path(qcdir,"ATAC_QC_summary_plots.pdf"),width=8,height=6)
if(nrow(out)){
  par(mar=c(8,4,3,1)); barplot(out$mapping_rate,names.arg=out$sample,las=2,ylab="Mapping rate",main="Alignment mapping rate")
  par(mar=c(8,4,3,1)); barplot(out$duplication_rate,names.arg=out$sample,las=2,ylab="Duplication rate",main="PCR duplication")
  par(mar=c(8,4,3,1)); barplot(out$mitochondrial_fraction,names.arg=out$sample,las=2,ylab="Mitochondrial fraction",main="Mitochondrial reads")
  par(mar=c(8,4,3,1)); barplot(out$FRiP,names.arg=out$sample,las=2,ylab="FRiP",main="Fraction of reads/fragments in peaks")
}
dev.off()
