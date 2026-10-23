#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_HOME="$(cd "$(dirname "$0")" && pwd)"
export CONFIG="${CONFIG:-${SCRIPT_HOME}/config.sh}"
source "${SCRIPT_HOME}/utils/common.sh"
load_config

mkdirs "$REF_DIR" "$FASTQC_DIR" "$CLEAN_DIR" "$ALIGN_DIR" "$FILTER_DIR" "$QC_DIR" "$BIGWIG_DIR" "$PEAK_DIR" "$CONSENSUS_DIR" "$COUNT_DIR" "$DIFF_DIR" "$ANNO_DIR" "$MOTIF_DIR" "$FOOTPRINT_DIR" "$LOG_DIR"

steps=(
  00_prepare_reference.sh
  01_fastqc.sh
  02_fastp.sh
  03_align.sh
  04_filter_bam.sh
  05_library_qc.sh
  06_bigwig_tss.sh
  07_call_peaks.sh
  08_consensus_counts.sh
  09_diff_annotation_motif.sh
  10_qc_summary.sh
)

"${SCRIPT_HOME}/check_env.sh" 2>&1 | tee "${LOG_DIR}/00_check_env.log"
for s in "${steps[@]}"; do
  log "===== Running ${s} ====="
  bash "${SCRIPT_HOME}/steps/${s}" 2>&1 | tee "${LOG_DIR}/${s%.sh}.log"
done
log "ATAC-seq pipeline completed successfully."
