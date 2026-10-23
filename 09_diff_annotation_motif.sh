#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd Rscript; mkdirs "$DIFF_DIR" "$ANNO_DIR" "$MOTIF_DIR"

Rscript "${SCRIPT_DIR}/R/diff_analysis.R" \
  "${COUNT_DIR}/ATAC_peak_counts.txt" "$SAMPLE_SHEET" "$DIFF_DIR" \
  "$MIN_TOTAL_COUNT" "$DESEQ2_PADJ" "$DESEQ2_ABS_LOG2FC"

sp=$(infer_species); asm=$(infer_assembly)
if [[ "$HOMER_GENOME" == auto ]]; then hgen="$asm"; else hgen="$HOMER_GENOME"; fi

shopt -s nullglob
for csv in "${DIFF_DIR}"/*.gained.csv "${DIFF_DIR}"/*.lost.csv; do
  base=$(basename "$csv" .csv)
  bed="${DIFF_DIR}/${base}.bed"
  Rscript "${SCRIPT_DIR}/R/make_diff_bed.R" "$csv" "${CONSENSUS_DIR}/consensus_peaks.bed" "$bed"
  [[ -s "$bed" ]] || { warn "No peaks in ${base}; skipping annotation/motif."; continue; }

  if [[ "$RUN_CHIPSEEKER" == true ]]; then
    if Rscript "${SCRIPT_DIR}/R/peak_annotation.R" "$bed" "$sp" "$asm" "${ANNO_DIR}/${base}.annotation.csv"; then
      require_file "${ANNO_DIR}/${base}.annotation.csv"
      log "Annotated ${base}"
    else
      die "ChIPseeker annotation failed for ${base}; inspect the R error before continuing."
    fi
  fi

  if [[ "$RUN_HOMER" == true ]]; then
    require_cmd findMotifsGenome.pl
    if [[ "$hgen" == custom || "$hgen" == auto ]]; then
      warn "HOMER genome could not be inferred for ${base}; set HOMER_GENOME in config.sh."
    else
      mkdir -p "${MOTIF_DIR}/${base}"
      findMotifsGenome.pl "$bed" "$hgen" "${MOTIF_DIR}/${base}" -size given -p "$THREADS"
      require_file "${MOTIF_DIR}/${base}/knownResults.txt"
      require_file "${MOTIF_DIR}/${base}/homerResults.html"
    fi
  fi
done

if [[ "$RUN_TOBIAS" == true ]]; then
  require_cmd TOBIAS
  warn "TOBIAS is enabled, but motif-file/database choice is project-specific. Run the template in README.md after setting a JASPAR motif file."
fi
