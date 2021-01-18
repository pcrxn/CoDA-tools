# CoDA-tools
Scripts for helping with compositional data analysis.

## `sorted-bam-to-fold-coverage.sh`

This script takes a sorted BAM file as input and uses `SAMtools` (v1.11+) and shell commands to generate a tab-separated table of features with their respective fold-coverages (`_fold-coverages.tsv`), along with alignment summary statistics (`_idxstats.tsv`) and read depth \[options: -a\] for each alignment base pair position (`_depth.tsv`). For more information, see the SAMtools [idxstats](https://www.htslib.org/doc/samtools-idxstats.html) and [depth](http://www.htslib.org/doc/samtools-depth.html) manual pages.
