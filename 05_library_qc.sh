#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd samtools; require_cmd bedtools; mkdirs "$QC_DIR"

printf 'sample\tlibrary_type\traw_aligned_records\tdedup_records_before_mt\tfinal_records\tmt_records_before_removal\tfragment_file\n' > "${QC_DIR}/library_qc.tsv"
while IFS=$'\t' read -r sample group r1 r2 batch; do
  lt=$(library_type_for "$r2")
  rawbam="${ALIGN_DIR}/${sample}.sorted.bam"
  finalbam="${FILTER_DIR}/${sample}.clean.sorted.bam"
  rawn=$(samtools view -@ "$THREADS" -c "$rawbam")
  finaln=$(samtools view -@ "$THREADS" -c "$finalbam")
  dedupn=$(samtools view -@ "$THREADS" -c "${FILTER_DIR}/${sample}.dedup.sorted.bam")
  mt=0
  mtname=$(mt_contig "${FILTER_DIR}/${sample}.dedup.sorted.bam" || true)
  if [[ -n "$mtname" ]]; then mt=$(samtools idxstats "${FILTER_DIR}/${sample}.dedup.sorted.bam" | awk -v x="$mtname" '$1==x{print $3}'); fi
  fragfile="NA"
  if [[ "$lt" == PE ]]; then
    fragfile="${QC_DIR}/${sample}.fragment_length.txt"
    samtools view -@ "$THREADS" "$finalbam" | awk '$9>0{print $9}' > "$fragfile"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$lt" "$rawn" "$dedupn" "$finaln" "$mt" "$fragfile" >> "${QC_DIR}/library_qc.tsv"
done < <(iter_samples)
