#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd bamCoverage; require_cmd samtools; mkdirs "$BIGWIG_DIR" "${QC_DIR}/tss"
if [[ "$USE_TN5_SHIFT_FOR_SIGNAL" == true ]]; then require_cmd alignmentSieve; fi

# deepTools may concatenate thousands of temporary BAM chunks.
soft=$(ulimit -Sn); hard=$(ulimit -Hn)
if [[ "$soft" != unlimited ]] && (( soft < 65536 )); then
  target=65536
  if [[ "$hard" != unlimited ]] && (( hard < target )); then target=$hard; fi
  ulimit -Sn "$target" || die "Cannot raise open-file soft limit"
fi
log "Open-file soft limit: $(ulimit -Sn)"
egs=$(effective_genome_size)
forrow(){ :; }
while IFS=$'\t' read -r sample group r1 r2 batch; do
  bam="${FILTER_DIR}/${sample}.clean.sorted.bam"
  signalbam="$bam"
  if [[ "$USE_TN5_SHIFT_FOR_SIGNAL" == true ]]; then
    lt=$(library_type_for "$r2")
    if [[ "$lt" == PE ]]; then
      alignmentSieve --ATACshift -b "$bam" -o "${FILTER_DIR}/${sample}.tn5shift.bam" -p "$THREADS"
      samtools sort -@ "$THREADS" -o "${FILTER_DIR}/${sample}.tn5shift.sorted.bam" "${FILTER_DIR}/${sample}.tn5shift.bam"
      samtools index -@ "$THREADS" "${FILTER_DIR}/${sample}.tn5shift.sorted.bam"
      signalbam="${FILTER_DIR}/${sample}.tn5shift.sorted.bam"
    else
      warn "alignmentSieve --ATACshift is intended for paired-end ATAC; using unshifted BAM for SE sample ${sample}."
    fi
  fi

  args=( -b "$signalbam" -o "${BIGWIG_DIR}/${sample}.${BIGWIG_NORMALIZATION}.bw" --binSize "$BIGWIG_BIN_SIZE" --normalizeUsing "$BIGWIG_NORMALIZATION" -p "$THREADS" )
  if [[ "$BIGWIG_NORMALIZATION" == RPGC ]]; then
    [[ -n "$egs" ]] || die "RPGC requested but effective genome size cannot be inferred. Set EFFECTIVE_GENOME_SIZE."
    args+=( --effectiveGenomeSize "$egs" )
  fi
  bamCoverage "${args[@]}"
done < <(iter_samples)

TSS_BED="${REF_DIR}/TSS.bed"
if [[ -s "$TSS_BED" ]]; then
  require_cmd computeMatrix; require_cmd plotProfile
  mapfile -t bws < <(find "$BIGWIG_DIR" -maxdepth 1 -name "*.${BIGWIG_NORMALIZATION}.bw" | sort)
  computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 -R "$TSS_BED" -S "${bws[@]}" \
    --skipZeros -p "$THREADS" -o "${QC_DIR}/tss/TSS_matrix.gz"
  plotProfile -m "${QC_DIR}/tss/TSS_matrix.gz" -out "${QC_DIR}/tss/TSS_profile.pdf" --perGroup
else
  warn "TSS.bed unavailable; TSS profile skipped. Provide GTF to auto-generate it."
fi
