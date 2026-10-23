#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
if [[ "$REMOVE_BLACKLIST" == true ]]; then require_file "$BLACKLIST_BED"; fi
require_cmd bedtools; require_cmd featureCounts; mkdirs "$CONSENSUS_DIR" "$COUNT_DIR"

require_file "${REFERENCE_FASTA}.fai"
all="${CONSENSUS_DIR}/all_sample_peaks.bed"
: > "$all"
if (( CONSENSUS_SUMMIT_HALF_WIDTH > 0 )); then
  while IFS=$'\t' read -r sample group r1 r2 batch; do
    sfile="${PEAK_DIR}/${sample}/${sample}_summits.bed"
    require_file "$sfile"
    use_summits="$sfile"
    if [[ "$REMOVE_BLACKLIST" == true ]]; then
      bedtools intersect -v -a "$sfile" -b "$BLACKLIST_BED" > "${CONSENSUS_DIR}/${sample}.summits.filtered.bed"
      use_summits="${CONSENSUS_DIR}/${sample}.summits.filtered.bed"
    fi
    awk -v w="$CONSENSUS_SUMMIT_HALF_WIDTH" 'BEGIN{OFS="\t"}{mid=$2; s=mid-w; if(s<0)s=0; e=mid+w; print $1,s,e,"'"$sample"'"}' "$use_summits" >> "$all"
  done < <(iter_samples)
else
  while IFS=$'\t' read -r sample group r1 r2 batch; do
    p="${PEAK_DIR}/${sample}/${sample}_peaks.filtered.narrowPeak"; require_file "$p"
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"'"$sample"'"}' "$p" >> "$all"
  done < <(iter_samples)
fi
sort -k1,1 -k2,2n "$all" -o "$all"

# Merge while retaining unique sample names, then keep intervals supported by >= N samples.
bedtools merge -i "$all" -c 4 -o distinct > "${CONSENSUS_DIR}/merged_with_samples.bed"
awk -F'\t' -v n="$CONSENSUS_MIN_SAMPLES" 'BEGIN{OFS="\t"}{x=$4; gsub(/[^,]/,"",x); ns=length(x)+1; if(ns>=n) print $1,$2,$3,"peak_"++k,ns,$4}' \
  "${CONSENSUS_DIR}/merged_with_samples.bed" > "${CONSENSUS_DIR}/consensus_peaks.meta.bed"
if [[ "$REMOVE_BLACKLIST" == true ]]; then
  bedtools intersect -v -a "${CONSENSUS_DIR}/consensus_peaks.meta.bed" -b "$BLACKLIST_BED" > "${CONSENSUS_DIR}/consensus_peaks.blacklist_filtered.bed"
  cp "${CONSENSUS_DIR}/consensus_peaks.blacklist_filtered.bed" "${CONSENSUS_DIR}/consensus_peaks.meta.bed"
fi
# Clip to the reference boundaries after expansion/merging; discard empty intervals.
awk 'BEGIN{FS=OFS="\t"} NR==FNR{len[$1]=$2;next} ($1 in len){if($2<0)$2=0;if($3>len[$1])$3=len[$1];if($3>$2)print}' "${REFERENCE_FASTA}.fai" "${CONSENSUS_DIR}/consensus_peaks.meta.bed" > "${CONSENSUS_DIR}/consensus_peaks.bounded.bed"
cp "${CONSENSUS_DIR}/consensus_peaks.bounded.bed" "${CONSENSUS_DIR}/consensus_peaks.meta.bed"
cut -f1-4 "${CONSENSUS_DIR}/consensus_peaks.meta.bed" > "${CONSENSUS_DIR}/consensus_peaks.bed"
require_file "${CONSENSUS_DIR}/consensus_peaks.bed"

awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}{print $4,$1,$2+1,$3,"."}' \
  "${CONSENSUS_DIR}/consensus_peaks.bed" > "${COUNT_DIR}/consensus_peaks.saf"

mapfile -t bams < <(while IFS=$'\t' read -r sample group r1 r2 batch; do echo "${FILTER_DIR}/${sample}.clean.sorted.bam"; done < <(iter_samples))
# If all samples are PE, count fragments; otherwise run read counting to avoid invalid mixed mode.
allpe=1
while IFS=$'\t' read -r sample group r1 r2 batch; do [[ "$(library_type_for "$r2")" == PE ]] || allpe=0; done < <(iter_samples)
args=( -F SAF -a "${COUNT_DIR}/consensus_peaks.saf" -o "${COUNT_DIR}/ATAC_peak_counts.txt" -T "$THREADS" )
if (( allpe == 1 )); then args+=( -p --countReadPairs ); fi
featureCounts "${args[@]}" "${bams[@]}"
