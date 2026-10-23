#!/usr/bin/env bash

set -Eeuo pipefail

# ---- Project ----
WORKDIR="/data/wujiaxin/ATAC-seq"
SAMPLE_SHEET="${WORKDIR}/samples.tsv"
RAW_DIR="${WORKDIR}/01.rawdata"
THREADS=12

# ---- Reference ----
# SPECIES: auto | human | mouse
SPECIES="human"
# ASSEMBLY: auto | hg38 | hg19 | mm39 | mm10 | custom
ASSEMBLY="hg38"
REFERENCE_FASTA="${WORKDIR}/00.reference/genome.fa"
GTF="${WORKDIR}/00.reference/genes.gtf"                 
BLACKLIST_BED="${WORKDIR}/00.reference/blacklist.bed"   
BOWTIE2_INDEX_PREFIX="${WORKDIR}/00.reference/bowtie2/genome"

# ---- Library / preprocessing ----
# LIBRARY_TYPE: auto | PE | SE ; auto determines per sample from samples.tsv R2 column
LIBRARY_TYPE="PE"
FASTP_QUAL=20
FASTP_MINLEN=20
BOWTIE2_MAX_INSERT=2000
MAPQ=30
REMOVE_DUPLICATES=true
REMOVE_MT=true
REMOVE_BLACKLIST=true

# ---- Signal / peak ----
USE_TN5_SHIFT_FOR_SIGNAL=true     # use deepTools alignmentSieve --ATACshift for BigWig/TSS signal only
BIGWIG_NORMALIZATION="CPM"        # CPM | RPGC
BIGWIG_BIN_SIZE=10
MACS_QVALUE=0.01
# Paired-end: BAMPE is used. Single-end: shift/extsize below are used.
SE_MACS_SHIFT=-100
SE_MACS_EXTSIZE=200
CALL_SUMMITS=true

# ---- Consensus peaks ----
# Keep merged intervals observed in at least this many samples.
CONSENSUS_MIN_SAMPLES=2
# 0 = use full merged intervals; >0 = make fixed-width windows around MACS summits (e.g. 250 => summit +/-250 bp)
CONSENSUS_SUMMIT_HALF_WIDTH=250

# ---- Differential accessibility ----
MIN_TOTAL_COUNT=10
DESEQ2_PADJ=0.05
DESEQ2_ABS_LOG2FC=1.0
# Optional batch column in samples.tsv. If present and non-empty, DESeq2 uses ~ batch + group.

# ---- Motif / annotation ----
RUN_CHIPSEEKER=true
RUN_HOMER=true
RUN_TOBIAS=false
HOMER_GENOME="/data/wujiaxin/ATAC-seq/00.reference/genome.fa"

# ---- QC thresholds ----
WARN_MAPPING_RATE=0.70
WARN_DUP_RATE=0.50
WARN_MT_RATE=0.30
WARN_FRIP=0.20

# ---- Derived directories----
REF_DIR="${WORKDIR}/00.reference"
FASTQC_DIR="${WORKDIR}/02.fastqc"
CLEAN_DIR="${WORKDIR}/03.clean"
ALIGN_DIR="${WORKDIR}/04.align"
FILTER_DIR="${WORKDIR}/05.filter"
QC_DIR="${WORKDIR}/06.qc"
BIGWIG_DIR="${WORKDIR}/07.bigwig"
PEAK_DIR="${WORKDIR}/08.peaks"
CONSENSUS_DIR="${WORKDIR}/09.consensus"
COUNT_DIR="${WORKDIR}/10.count"
DIFF_DIR="${WORKDIR}/11.diff"
ANNO_DIR="${WORKDIR}/12.annotation"
MOTIF_DIR="${WORKDIR}/13.motif"
FOOTPRINT_DIR="${WORKDIR}/14.footprint"
LOG_DIR="${WORKDIR}/logs"
SCRIPT_DIR="${WORKDIR}/scripts"
