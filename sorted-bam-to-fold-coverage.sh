#!/bin/bash

<<COMMENT

----------------------------
SORTED-BAM-TO-FOLD-COVERAGE
----------------------------

DESCRIPTION:

    Takes a sorted BAM file as input and uses SAMtools (v1.11+) and shell
    commands to generate a tab-separated table of features with their respective
    fold-coverages (_fold-coverages.tsv), along with alignment summary
    statistics (_idxstats.tsv) and read depth for each alignment base pair
    position (_depth.tsv).

    WARNING: This script will remove any files ending with
    ".tmp" in the output directory provided.

USAGE:

    ./sorted-bam-to-fold-coverage.sh input.bam /path/to/output/dir/

EXAMPLE:

    ./sorted-bam-told-coverage.sh bowtie2/Control_55_S1_sorted.bam results/

COMMENT

#-------------------------------------------------------------------------------
# Set global variables
#-------------------------------------------------------------------------------

# Set the input BAM file to the first argument provided on the command-line
BAM="$1"

# Set the output directory to the second argument provided on the command-line
OUTDIR="$2"

#-------------------------------------------------------------------------------
# usage()
#-------------------------------------------------------------------------------

usage(){

  echo
  echo "----------------------------"
  echo "SORTED-BAM-TO-FOLD-COVERAGE"
  echo "----------------------------"
  echo
  echo "DESCRIPTION:"
  echo
  echo "  Takes a sorted BAM file as input and uses SAMtools (v1.11+) and shell"
  echo "  commands to generate a tab-separated table of features with their"
  echo "  respective fold-coverages (_fold-coverages.tsv), along with alignment"
  echo "  summary statistics (_idxstats.tsv) and read depth for each alignment"
  echo "  base pair position (_depth.tsv)."
  echo
  echo "  WARNING: This script will remove any files ending with \".tmp\" in"
  echo "  the output directory provided."
  echo
  echo "USAGE:"
  echo
  echo "  ./sorted-bam-to-fold-coverage.sh input.bam /path/to/output/dir/"
  echo
  echo "EXAMPLE:"
  echo
  echo "  ./sorted-bam-to-fold-coverage.sh bowtie2/Control_55_S1_sorted.bam \\"
  echo "     results/"
  echo

}

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

main(){

  # Exit when any command fails
  set -e

  # Obtain the base file name of the input BAM file
  # The base file name will be everything that comes before the first "."
  FILE_NAME="$(echo $BAM | sed 's/.*\///' | sed 's/\..*//')"

  # echo "Main function is running on input BAM file \"$BAM\" and sending output to \"$OUTDIR\""
  # echo $FILE_NAME

  # Generate idxstats file
  echo "Generating samtools alignment summary statistics file..."
  samtools idxstats $BAM > \
    $OUTDIR/${FILE_NAME}_idxstats.tsv

  # Generate depth file
  echo "Generating samtools depth file..."
  samtools depth -a $BAM > \
    $OUTDIR/${FILE_NAME}_depth.tsv

  echo "Processing depth and alignment summary statistics..."

  # Sum depths for each reference sequence and output to new TSV
  awk '{a[$1]+=$3}END{for(i in a) print i,"\t",a[i]}' \
    $OUTDIR/${FILE_NAME}_depth.tsv \
    > $OUTDIR/${FILE_NAME}_depth-sum.tsv.tmp

  # Sort by first column of each table
  sort -k 1 $OUTDIR/${FILE_NAME}_depth-sum.tsv.tmp \
      > $OUTDIR/${FILE_NAME}_depth-sum-sorted.tsv.tmp
  sort -k 1 $OUTDIR/${FILE_NAME}_idxstats.tsv \
      > $OUTDIR/${FILE_NAME}_idxstats.tsv.tmp

  # Join together the sorted tables by the first columns and retain ONLY the rows
  # that are present in both files
  join \
    <(dos2unix <$OUTDIR/${FILE_NAME}_depth-sum-sorted.tsv.tmp) \
    <(dos2unix <$OUTDIR/${FILE_NAME}_idxstats.tsv.tmp) \
    > $OUTDIR/${FILE_NAME}_table.tsv.tmp

  # Calculate fold-coverages and output to final TSV file
  awk '{print $1,"\t",$2/$3}' $OUTDIR/${FILE_NAME}_table.tsv.tmp \
    > $OUTDIR/${FILE_NAME}_fold-coverages.tsv

  echo "Fold-coverage table generated"
  rm $OUTDIR/*.tmp
  echo "Removed all .tmp files in $OUTDIR"
  exit 0

}

#-------------------------------------------------------------------------------
# Execute
#-------------------------------------------------------------------------------

# If less than two arguments provided:
if [ $# -lt 2 ]; then
  echo
  echo "$# arguments provided on the command-line (two required)"
  usage
  exit 1
fi

main
