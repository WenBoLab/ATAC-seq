#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd fastp; mkdirs "$CLEAN_DIR"
while IFS=$'\t' read -r sample group r1 r2 batch; do
  lt=$(library_type_for "$r2")
  log "fastp ${sample} (${lt})"
  if [[ "$lt" == PE ]]; then
    fastp -i "${RAW_DIR}/${r1}" -I "${RAW_DIR}/${r2}" \
      -o "${CLEAN_DIR}/${sample}_R1.clean.fastq.gz" -O "${CLEAN_DIR}/${sample}_R2.clean.fastq.gz" \
      --detect_adapter_for_pe --qualified_quality_phred "$FASTP_QUAL" --length_required "$FASTP_MINLEN" \
      --thread "$THREADS" --html "${CLEAN_DIR}/${sample}.fastp.html" --json "${CLEAN_DIR}/${sample}.fastp.json"
  else
    fastp -i "${RAW_DIR}/${r1}" -o "${CLEAN_DIR}/${sample}.clean.fastq.gz" \
      --qualified_quality_phred "$FASTP_QUAL" --length_required "$FASTP_MINLEN" \
      --thread "$THREADS" --html "${CLEAN_DIR}/${sample}.fastp.html" --json "${CLEAN_DIR}/${sample}.fastp.json"
  fi
done < <(iter_samples)
