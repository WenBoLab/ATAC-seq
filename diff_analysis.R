#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) stop("Usage: diff_analysis.R counts.txt samples.tsv outdir min_total padj abs_log2fc")
count_file <- args[1]; sample_file <- args[2]; outdir <- args[3]
min_total <- as.numeric(args[4]); padj_cut <- as.numeric(args[5]); lfc_cut <- ifelse(length(args)>=6, as.numeric(args[6]), 1)
suppressPackageStartupMessages(library(DESeq2))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

fc <- read.delim(count_file, check.names=FALSE, comment.char="#")
meta <- read.delim(sample_file, check.names=FALSE, stringsAsFactors=FALSE)
count_cols <- 7:ncol(fc)
mat <- as.matrix(fc[,count_cols,drop=FALSE]); rownames(mat) <- fc$Geneid
colnames(mat) <- sub("\\.clean\\.sorted\\.bam$", "", basename(colnames(mat)))
missing <- setdiff(meta$sample, colnames(mat)); if(length(missing)) stop("Samples missing from featureCounts matrix: ", paste(missing, collapse=","))
mat <- mat[, meta$sample, drop=FALSE]
meta <- meta[match(colnames(mat), meta$sample),,drop=FALSE]; rownames(meta) <- meta$sample
meta$group <- factor(meta$group)
if("batch" %in% names(meta) && any(!is.na(meta$batch) & meta$batch!="") && any(is.na(meta$batch) | meta$batch=="")) stop("Partially missing batch labels")
use_batch <- "batch" %in% colnames(meta) && all(!is.na(meta$batch) & meta$batch!="") && length(unique(meta$batch))>1
if(any(table(meta$group)<2)) stop("At least two samples per group are required; verify biological independence separately")
if (use_batch) { meta$batch <- factor(meta$batch); design_formula <- ~ batch + group } else design_formula <- ~ group

if(any(!is.finite(mat)) || any(mat<0) || any(mat!=round(mat))) stop("Counts must be finite nonnegative integers")
if(anyDuplicated(meta$sample) || anyDuplicated(rownames(mat))) stop("Duplicate sample or peak identifiers")
mm <- model.matrix(design_formula,meta)
if(qr(mm)$rank<ncol(mm)) stop("Design is not full rank: check group/batch confounding")
dds <- DESeqDataSetFromMatrix(mat, meta, design_formula)
dds <- dds[rowSums(counts(dds)) >= min_total,]
if(nrow(dds)==0) stop("No peaks remain after count filtering")
dds <- DESeq(dds)
vsd <- if(sum(rowMeans(counts(dds, normalized=TRUE)) > 5) >= 1000) vst(dds, blind=FALSE) else varianceStabilizingTransformation(dds, blind=FALSE)
pdf(file.path(outdir,"PCA.pdf"), width=7, height=6); print(plotPCA(vsd, intgroup=if(use_batch)c("group","batch") else "group")); dev.off()
write.csv(as.data.frame(assay(vsd)), file.path(outdir,"VST_matrix.csv"))

groups <- levels(meta$group)
if(length(groups)<2) stop("Need at least two groups for differential analysis")
for(i in 2:length(groups)) for(j in 1:(i-1)) {
  test <- groups[i]; ref <- groups[j]; nm <- paste0(test,"_vs_",ref)
  res <- results(dds, contrast=c("group",test,ref))
  out <- as.data.frame(res); out$peak_id <- rownames(out); out <- out[,c("peak_id",setdiff(names(out),"peak_id"))]
  out <- out[order(out$padj, na.last=TRUE),]
  write.csv(out, file.path(outdir,paste0(nm,".csv")), row.names=FALSE)
  sig <- subset(out, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut)
  gained <- subset(sig, log2FoldChange > 0); lost <- subset(sig, log2FoldChange < 0)
  write.csv(sig, file.path(outdir,paste0(nm,".significant.csv")), row.names=FALSE)
  write.csv(gained, file.path(outdir,paste0(nm,".gained.csv")), row.names=FALSE)
  write.csv(lost, file.path(outdir,paste0(nm,".lost.csv")), row.names=FALSE)
}
