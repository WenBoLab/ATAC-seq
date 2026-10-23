#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
if [[ "$REMOVE_BLACKLIST" == true ]]; then require_file "$BLACKLIST_BED"; fi
require_cmd python
python -c "import pysam"
require_cmd macs3; require_cmd samtools; require_cmd bedtools; mkdirs "$PEAK_DIR" "$QC_DIR"
mg=$(macs_genome)
printf 'sample\tpeaks\ttotal_units\tin_peak_units\tFRiP\tunit\n' > "${QC_DIR}/frip.tsv"

while IFS=$'\t' read -r sample group r1 r2 batch; do
  lt=$(library_type_for "$r2")
  bam="${FILTER_DIR}/${sample}.clean.sorted.bam"
  out="${PEAK_DIR}/${sample}"; mkdir -p "$out"
  common=( callpeak -t "$bam" -g "$mg" -n "$sample" --outdir "$out" -q "$MACS_QVALUE" --keep-dup all )
  [[ "$CALL_SUMMITS" == true ]] && common+=( --call-summits )
  if [[ "$lt" == PE ]]; then
    macs3 "${common[@]}" -f BAMPE
  else
    macs3 "${common[@]}" -f BAM --nomodel --shift "$SE_MACS_SHIFT" --extsize "$SE_MACS_EXTSIZE"
  fi

  peak="${out}/${sample}_peaks.narrowPeak"
  require_file "$peak"
  fpeak="${out}/${sample}_peaks.filtered.narrowPeak"
  if [[ "$REMOVE_BLACKLIST" == true ]]; then bedtools intersect -v -a "$peak" -b "$BLACKLIST_BED" > "$fpeak"; else cp "$peak" "$fpeak"; fi
  npeaks=$(wc -l < "$fpeak")

  # Count each complete physical fragment once; BAM filtering guarantees mate integrity.
  python "${SCRIPT_DIR}/utils/calculate_frip.py" "$bam" "$fpeak" "$sample" "$lt" >> "${QC_DIR}/frip.tsv"
done < <(iter_samples)
