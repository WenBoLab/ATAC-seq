#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd bowtie2; require_cmd samtools; mkdirs "$ALIGN_DIR"
while IFS=$'\t' read -r sample group r1 r2 batch; do
  lt=$(library_type_for "$r2")
  log "Aligning ${sample} (${lt})"
  if [[ "$lt" == PE ]]; then
    bowtie2 --rg-id "$sample" --rg "SM:$sample" --rg "LB:$sample" --rg "PL:ILLUMINA" -x "$BOWTIE2_INDEX_PREFIX" \
      -1 "${CLEAN_DIR}/${sample}_R1.clean.fastq.gz" -2 "${CLEAN_DIR}/${sample}_R2.clean.fastq.gz" \
      --very-sensitive -X "$BOWTIE2_MAX_INSERT" --no-mixed --no-discordant -p "$THREADS" \
      2> "${ALIGN_DIR}/${sample}.bowtie2.log" | \
      samtools sort -@ "$THREADS" -o "${ALIGN_DIR}/${sample}.sorted.bam" -
  else
    bowtie2 --rg-id "$sample" --rg "SM:$sample" --rg "LB:$sample" --rg "PL:ILLUMINA" -x "$BOWTIE2_INDEX_PREFIX" -U "${CLEAN_DIR}/${sample}.clean.fastq.gz" \
      --very-sensitive -p "$THREADS" 2> "${ALIGN_DIR}/${sample}.bowtie2.log" | \
      samtools sort -@ "$THREADS" -o "${ALIGN_DIR}/${sample}.sorted.bam" -
  fi
  samtools index -@ "$THREADS" "${ALIGN_DIR}/${sample}.sorted.bam"
  samtools flagstat -@ "$THREADS" "${ALIGN_DIR}/${sample}.sorted.bam" > "${ALIGN_DIR}/${sample}.flagstat.txt"
done < <(iter_samples)
