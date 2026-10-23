#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "$0")/.." && pwd)/utils/common.sh"; load_config
require_cmd samtools; require_cmd bowtie2-build; require_file "$REFERENCE_FASTA"
mkdirs "$(dirname "$BOWTIE2_INDEX_PREFIX")" "$REF_DIR"
[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"
complete=0; partial=0
for ext in bt2 bt2l; do
  n=0
  for part in 1 2 3 4 rev.1 rev.2; do
    if [[ -s "${BOWTIE2_INDEX_PREFIX}.${part}.${ext}" ]]; then n=$((n+1)); fi
  done
  if (( n==6 )); then complete=1; fi
  if (( n>0 && n<6 )); then partial=1; fi
done
if (( complete==0 )); then
  (( partial==0 )) || die "Incomplete Bowtie2 index; preserve and relocate partial shards before rebuilding"
  bowtie2-build "$REFERENCE_FASTA" "$BOWTIE2_INDEX_PREFIX"
else
  log "All six Bowtie2 index shards found."
fi

# Generate TSS BED from GTF if possible. BED is 0-based half-open, one-base TSS intervals.
TSS_BED="${REF_DIR}/TSS.bed"
if [[ -s "$GTF" && ! -s "$TSS_BED" ]]; then
  log "Generating TSS BED from GTF..."
  awk 'BEGIN{OFS="\t"} $3=="transcript" {
    if($7=="+") {s=$4-1; if(s<0)s=0; e=$4} else {s=$5-1; if(s<0)s=0; e=$5}
    print $1,s,e,".",".",$7
  }' "$GTF" | sort -k1,1 -k2,2n -k3,3n -k6,6 | uniq | awk 'BEGIN{FS=OFS="\t"}{$4="TSS_"NR; print}' > "$TSS_BED"
fi
log "Reference ready. species=$(infer_species), assembly=$(infer_assembly)"
