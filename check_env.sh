#!/usr/bin/env bash
set -Eeuo pipefail
SCRIPT_HOME="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${CONFIG:-${SCRIPT_HOME}/config.sh}"
source "${SCRIPT_HOME}/utils/common.sh"
load_config
if [[ "$REMOVE_BLACKLIST" == true ]]; then require_file "$BLACKLIST_BED"; fi

req=(fastqc multiqc fastp bowtie2 bowtie2-build samtools bedtools picard macs3 bamCoverage computeMatrix plotProfile featureCounts Rscript python)
for x in "${req[@]}"; do require_cmd "$x"; done
if [[ "$USE_TN5_SHIFT_FOR_SIGNAL" == true ]]; then require_cmd alignmentSieve; fi
if [[ "$RUN_HOMER" == true ]]; then require_cmd findMotifsGenome.pl; fi
if [[ "$RUN_TOBIAS" == true ]]; then require_cmd TOBIAS; fi

python -c "import pysam" || die "Python pysam is required"
Rscript -e 'stopifnot(requireNamespace("DESeq2", quietly=TRUE))'
if [[ "$RUN_CHIPSEEKER" == true ]]; then
  Rscript -e 'stopifnot(requireNamespace("ChIPseeker", quietly=TRUE))'
fi
require_file "$REFERENCE_FASTA"
sample_header_check
[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"

sp=$(infer_species); asm=$(infer_assembly)
if [[ "$RUN_CHIPSEEKER" == true ]]; then
  case "$sp:$asm" in
    human:hg38) db=TxDb.Hsapiens.UCSC.hg38.knownGene; org=org.Hs.eg.db;;
    human:hg19) db=TxDb.Hsapiens.UCSC.hg19.knownGene; org=org.Hs.eg.db;;
    mouse:mm10) db=TxDb.Mmusculus.UCSC.mm10.knownGene; org=org.Mm.eg.db;;
    *) die "ChIPseeker automatic annotation unsupported for $sp:$asm";;
  esac
  Rscript -e 'for(p in commandArgs(TRUE)) if(!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)' "$db" "$org"
fi
log "Detected species=${sp}; assembly=${asm}"
[[ "$sp" != unknown ]] || warn "Species could not be inferred; set SPECIES=human/mouse or MACS_GENOME_SIZE manually."

while IFS=$'\t' read -r sample group r1 r2 batch; do
  [[ -n "$sample" && -n "$group" && -n "$r1" ]] || die "Invalid row in samples.tsv"
  require_file "${RAW_DIR}/${r1}"
  lt=$(library_type_for "$r2")
  if [[ "$lt" == PE ]]; then require_file "${RAW_DIR}/${r2}"; fi
  log "Sample ${sample}: group=${group}, library=${lt}${batch:+, batch=${batch}}"
done < <(iter_samples)

log "Environment and input checks passed."
