#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd fastqc; require_cmd multiqc; sample_header_check
mkdirs "${FASTQC_DIR}/raw"
files=()
while IFS=$'\t' read -r sample group r1 r2 batch; do
  files+=("${RAW_DIR}/${r1}")
  [[ "$(library_type_for "$r2")" == PE ]] && files+=("${RAW_DIR}/${r2}")
done < <(iter_samples)
fastqc -t "$THREADS" -o "${FASTQC_DIR}/raw" "${files[@]}"
multiqc -f "${FASTQC_DIR}/raw" -o "${FASTQC_DIR}/raw/multiqc"
