#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd Rscript; require_cmd multiqc
Rscript "${SCRIPT_DIR}/R/qc_summary.R" "$WORKDIR" "$SAMPLE_SHEET" \
  "$WARN_MAPPING_RATE" "$WARN_DUP_RATE" "$WARN_MT_RATE" "$WARN_FRIP" "$QC_DIR" "$ALIGN_DIR"
multiqc -f "$WORKDIR" -o "${QC_DIR}/multiqc_final" || die "Final MultiQC failed; outputs retained for diagnosis."
require_file "${QC_DIR}/multiqc_final/multiqc_report.html"
log "QC summary: ${QC_DIR}/ATAC_QC_summary.csv"
