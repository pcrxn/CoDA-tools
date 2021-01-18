# CoDA-tools

Scripts for helping with the compositional data analysis of microbiome data sets.

## `sorted-bam-to-fold-coverage.sh`

This script takes a sorted BAM file as input and uses [SAMtools](http://samtools.sourceforge.net/) (v1.11+) and shell commands to generate a tab-separated table of features with their respective fold-coverages (`_fold-coverages.tsv`), along with alignment summary statistics (`_idxstats.tsv`) and read depth \[options: -a\] for each alignment base pair position (`_depth.tsv`). For more information, see the SAMtools [idxstats](https://www.htslib.org/doc/samtools-idxstats.html) and [depth](http://www.htslib.org/doc/samtools-depth.html) manual pages.

### Usage

```bash
$ ./sorted-bam-to-fold-coverage.sh input.bam /path/to/output/dir/
```

A **warning**: This script will remove any files ending with .tmp in the output directory provided. If you're currently running other scripts on this directory, I'd recommend choosing another output directory or wait until other scripts are finished running.

Run the script with no arguments (`./sorted-bam-to-fold-coverage.sh`) for a help message.

### Output

If the input BAM file is named `input.bam`:

* `input_idxstats.tsv`: A tab-separated table describing the alignment summary statistics. See SAMtools [`idxstats`](https://www.htslib.org/doc/samtools-idxstats.html) for more details.
* `input_depth.tsv`: A tab-separated table describing read depth for each position in each alignment; positions with no alignments are included as '0'. See SAMtools [`depth`](http://www.htslib.org/doc/samtools-depth.html) for more details.
* `input_fold-coverages.tsv`: A tab-separated table with the features in the BAM file as the first column and fold-coverage as the second column. Fold-coverage is calculated by calculating the cumulative read depth for each feature and dividing by the length of the reference sequence of that feature.

### Example

To run the script on an input file named `Control_55_S1_sorted.bam` in `bowtie2/`, and to send the output to `results/`:

```bash
$ ./sorted-bam-to-fold-coverage.sh bowtie2/Control_55_S1_sorted.bam results/
```

The output files will be named `Control_55_S1_sorted_idxstats.tsv`, `Control_55_S1_sorted_depth.tsv`, and `Control_55_S1_sorted_fold-coverages.tsv`.
