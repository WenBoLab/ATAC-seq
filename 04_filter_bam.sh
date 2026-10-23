#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
if [[ "$REMOVE_BLACKLIST" == true ]]; then require_file "$BLACKLIST_BED"; fi
require_cmd python
python -c "import pysam"
require_cmd samtools; require_cmd picard; require_cmd bedtools; mkdirs "$FILTER_DIR" "$QC_DIR"

while IFS=$'\t' read -r sample group r1 r2 batch; do
  lt=$(library_type_for "$r2")
  inbam="${ALIGN_DIR}/${sample}.sorted.bam"
  require_file "$inbam"
  if ! samtools view -H "$inbam" | awk '$1=="@RG"{found=1} END{exit !found}'; then
    log "Adding read group to a separate copy for ${sample}"
    samtools addreplacerg -@ "$THREADS" -r "ID:${sample}" -r "SM:${sample}" -r "LB:${sample}" -r "PL:ILLUMINA" -o "${FILTER_DIR}/${sample}.rg.bam" "$inbam"
    inbam="${FILTER_DIR}/${sample}.rg.bam"
  fi
  log "Filtering ${sample} (${lt}) MAPQ>=${MAPQ}"

  if [[ "$lt" == PE ]]; then
    # exclude unmapped(4), mate unmapped(8), secondary(256), QC fail(512), supplementary(2048); require proper pair(2)
    samtools view -@ "$THREADS" -bh -f 2 -F 2828 -q "$MAPQ" "$inbam" > "${FILTER_DIR}/${sample}.predup.bam"
  else
    # exclude unmapped, secondary, QC fail, supplementary
    samtools view -@ "$THREADS" -bh -F 2820 -q "$MAPQ" "$inbam" > "${FILTER_DIR}/${sample}.predup.bam"
  fi

  predup="${FILTER_DIR}/${sample}.predup.bam"
  if [[ "$lt" == PE ]]; then
    complete_pairs "$predup" "${FILTER_DIR}/${sample}.predup.paired.bam" "${QC_DIR}/${sample}.predup_pair_audit.json"
    predup="${FILTER_DIR}/${sample}.predup.paired.bam"
  fi
  picard MarkDuplicates \
    I="$predup" \
    O="${FILTER_DIR}/${sample}.markdup.bam" \
    M="${QC_DIR}/${sample}.dup_metrics.txt" \
    REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

  if [[ "$REMOVE_DUPLICATES" == true ]]; then
    samtools view -@ "$THREADS" -bh -F 1024 "${FILTER_DIR}/${sample}.markdup.bam" > "${FILTER_DIR}/${sample}.dedup.bam"
  else
    cp "${FILTER_DIR}/${sample}.markdup.bam" "${FILTER_DIR}/${sample}.dedup.bam"
  fi
  samtools sort -@ "$THREADS" -o "${FILTER_DIR}/${sample}.dedup.sorted.bam" "${FILTER_DIR}/${sample}.dedup.bam"
  samtools index -@ "$THREADS" "${FILTER_DIR}/${sample}.dedup.sorted.bam"

  cur="${FILTER_DIR}/${sample}.dedup.sorted.bam"
  if [[ "$REMOVE_MT" == true ]]; then
    mt=$(mt_contig "$cur" || true)
    if [[ -n "$mt" ]]; then
      log "Removing mitochondrial contig ${mt} from ${sample}"
      samtools idxstats "$cur" | awk -v mt="$mt" '$1!=mt && $1!="*" {print $1}' | \
        xargs samtools view -@ "$THREADS" -bh "$cur" > "${FILTER_DIR}/${sample}.nuclear.bam"
      samtools sort -@ "$THREADS" -o "${FILTER_DIR}/${sample}.nuclear.sorted.bam" "${FILTER_DIR}/${sample}.nuclear.bam"
      samtools index -@ "$THREADS" "${FILTER_DIR}/${sample}.nuclear.sorted.bam"
      cur="${FILTER_DIR}/${sample}.nuclear.sorted.bam"
    else
      warn "No chrM/MT contig found in ${sample}; mitochondrial removal skipped."
    fi
  fi

  if [[ "$REMOVE_BLACKLIST" == true ]]; then
    if [[ -s "$BLACKLIST_BED" ]]; then
      bedtools intersect -v -abam "$cur" -b "$BLACKLIST_BED" > "${FILTER_DIR}/${sample}.clean.unsorted.bam"
      samtools sort -@ "$THREADS" -o "${FILTER_DIR}/${sample}.clean.stage.bam" "${FILTER_DIR}/${sample}.clean.unsorted.bam"
    else
      warn "BLACKLIST_BED missing; copying current BAM without blacklist filtering."
      cp "$cur" "${FILTER_DIR}/${sample}.clean.stage.bam"
    fi
  else
    cp "$cur" "${FILTER_DIR}/${sample}.clean.stage.bam"
  fi

  if [[ "$lt" == PE ]]; then
    complete_pairs "${FILTER_DIR}/${sample}.clean.stage.bam" "${FILTER_DIR}/${sample}.clean.sorted.bam" "${QC_DIR}/${sample}.final_pair_audit.json"
  else
    cp "${FILTER_DIR}/${sample}.clean.stage.bam" "${FILTER_DIR}/${sample}.clean.sorted.bam"
  fi
  samtools index -@ "$THREADS" "${FILTER_DIR}/${sample}.clean.sorted.bam"
  samtools flagstat -@ "$THREADS" "${FILTER_DIR}/${sample}.clean.sorted.bam" > "${QC_DIR}/${sample}.final.flagstat.txt"

done < <(iter_samples)
