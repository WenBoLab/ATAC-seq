# ATAC-seq executable pipeline

Target location: `/data/username/ATAC-seq/scripts/`

## 1\. Expected project layout



/data/username/ATAC-seq/
├── 00.reference/
│   ├── genome.fa
│   ├── genes.gtf                 # recommended
│   └── blacklist.bed             # recommended
├── 01.rawdata/
│   └── \*.fastq.gz
├── samples.tsv
└── scripts/                      # copy this package here



The pipeline creates `02.fastqc` through `14.footprint`, plus `logs`.

## 2\. Install software



cd /data/username/ATAC-seq/scripts
conda env create -f environment.yml
conda activate atacseq



## 3\. Prepare samples.tsv

Copy `samples.tsv.example` to `/data/usename/ATAC-seq/samples.tsv` and edit it.

Required columns: `sample`, `group`, `R1`, `R2`.
Optional column: `batch`.

For single-end data, set R2 to `NA` or `.`. PE/SE is automatically detected per row when `LIBRARY\\\_TYPE=auto`.

## 4\. Edit config.sh

The defaults are designed for standard bulk ATAC-seq.

Important choices:

* Human or mouse: `SPECIES=auto` detects common assembly names from FASTA paths and falls back to chr1 length. Set it explicitly if your FASTA name is unconventional.
* Assembly: `ASSEMBLY=auto` recognizes hg38/hg19/mm39/mm10 from the FASTA path. For other genomes, set `MACS\\\_GENOME\\\_SIZE` and optionally `EFFECTIVE\\\_GENOME\\\_SIZE`.
* PE ATAC: Bowtie2 uses `--very-sensitive -X 2000`; MACS3 uses `-f BAMPE`.
* SE ATAC: Bowtie2 uses `-U`; MACS3 uses `--nomodel --shift -100 --extsize 200` by default.
* Alignment filter: MAPQ >= 30. Set 20 for shallow data only when justified.
* Duplicates: removed by default after Picard marking.
* Mitochondrial reads: removed by default; contig name chrM/MT/M/chrMT is auto-detected.
* Blacklist: removed if `00.reference/blacklist.bed` exists.
* BigWig: CPM by default. RPGC requires a supported assembly or `EFFECTIVE\\\_GENOME\\\_SIZE`.
* Tn5 offset: applied only to PE signal tracks/TSS profiles through `alignmentSieve --ATACshift`; the primary clean BAM is not permanently modified.
* Consensus peaks: MACS summits are expanded +/-250 bp and merged, then retained if supported by at least 2 samples. Set `CONSENSUS\\\_MIN\\\_SAMPLES=1` when only one replicate per group exists.
* Differential peaks: DESeq2; if a non-empty `batch` column exists with >1 level, design is `\\\~ batch + group`; otherwise `\\\~ group`.
* Differential cutoffs: padj < 0.05 and |log2FC| >= 1 by default.

## 5\. Check, then run



cd /data/username/ATAC-seq/scripts
bash check\_env.sh
bash run\_atacseq.sh



For a long server job:



nohup bash run\_atacseq.sh > /data/username/ATAC-seq/logs/master.nohup.log 2>\&1 \&



## 6\. Main outputs

* Clean BAM: `05.filter/<sample>.clean.sorted.bam`
* Fragment length: `06.qc/<sample>.fragment\\\_length.txt` for PE
* QC summary: `06.qc/ATAC\\\_QC\\\_summary.csv`
* QC plots: `06.qc/ATAC\\\_QC\\\_summary\\\_plots.pdf`
* BigWig: `07.bigwig/<sample>.CPM.bw`
* Peaks: `08.peaks/<sample>/<sample>\\\_peaks.filtered.narrowPeak`
* Consensus peaks: `09.consensus/consensus\\\_peaks.bed`
* Count matrix: `10.count/ATAC\\\_peak\\\_counts.txt`
* DESeq2 contrasts: `11.diff/\\\*\\\_vs\\\_\\\*.csv`
* Gained/lost BEDs: `11.diff/\\\*.gained.bed`, `\\\*.lost.bed`
* ChIPseeker annotations: `12.annotation/`
* HOMER motif results: `13.motif/`

## 7\. Parameter choices by data type

|Data|Alignment|Peak calling|Counting|
|-|-|-|-|
|PE50/75/100/150 bulk ATAC|`--very-sensitive -X 2000 --no-mixed --no-discordant`|`MACS3 -f BAMPE`|`featureCounts -p --countReadPairs`|
|SE bulk ATAC|`--very-sensitive -U`|`--nomodel --shift -100 --extsize 200`|read counting, no PE flags|
|Human|auto -> `-g hs`|hg38/hg19 recognized|ChIPseeker supported for hg38/hg19|
|Mouse|auto -> `-g mm`|mm10/mm39 recognized|ChIPseeker bundled script supports mm10; mm39 requires custom TxDb|
|Shallow library|optionally MAPQ 20|keep q=0.01 initially|report sensitivity analysis if relaxing thresholds|
|No biological replicates|set consensus min to 1|peak calling works|formal DESeq2 inference is not reliable; treat as exploratory|

## 8\. Important interpretation notes

1. Biological replicates are never merged before DESeq2.
2. `-X 2000` is maximum paired fragment length, not read length; PE50 and PE150 can use the same value.
3. PE MACS3 uses observed fragments with `BAMPE`; do not mechanically combine it with SE shift/extsize settings.
4. Tn5 shift is used for signal tracks/cut-site-type analyses, not as a destructive rewrite of the main analysis BAM.
5. Nearest-gene peak annotation is not proof of enhancer-target regulation.
6. Gained and lost peaks are sent to motif analysis separately.

## 9\. Optional TOBIAS template

After choosing a motif file, e.g. JASPAR motifs:



TOBIAS ATACorrect \\
--bam 05.filter/sample.clean.sorted.bam \\
--genome 00.reference/genome.fa \\
--peaks 09.consensus/consensus\_peaks.bed \\
--outdir 14.footprint/sample \\
--cores 12

TOBIAS FootprintScores \\
--signal 14.footprint/sample/sample\_corrected.bw \\
--regions 09.consensus/consensus\_peaks.bed \\
--output 14.footprint/sample/sample\_footprints.bw \\
--cores 12



`BINDetect` should then be run with a project-appropriate motif database and named conditions.

