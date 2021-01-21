#!/bin/bash
#SBATCH --account=def-vtai4
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=kmer-trim_array
#SBATCH --output=kmer-trim_array_%j.out
#SBATCH --error=kmer-trim_array_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lbrow55@uwo.ca

# This script can be run using the following command:
# `sbatch --array=1-9 kmer-trim-submit.sh`
# 9 is the number of lines in the ids.txt document and
# the number of sample IDs

#-------------------------------------------------------------------------------
# Initialize environment
#-------------------------------------------------------------------------------

# Activate conda environment
module --force unload StdEnv/2020
module load nixpkgs/16.09 miniconda3
source activate sourmash

#-------------------------------------------------------------------------------
# State global variables
#-------------------------------------------------------------------------------

# Directory containing input interleaved .fastq.gz files
DATA_DIR="/home/brownl/scratch/MAC-HiSeq/data/combined/interleaved-fastq"

# Directory to output interleaved .fastq.gz files and error logs
OUT_DIR="/home/brownl/scratch/MAC-HiSeq/sourmash/data/kmer-trimmed"

#-------------------------------------------------------------------------------
# Retrieve sample IDs
#-------------------------------------------------------------------------------

cd $DATA_DIR

# Obtain the base file names (sample IDs) and write to a temporary file
ls -l *_001_P.fastq.gz | cut -d ' ' -f 9 | sed 's/\(.*\)\_001_P.fastq.gz/\1/g' > $SLURM_TMPDIR/ids.txt

#-------------------------------------------------------------------------------
# K-mer trim interleaved read pairs
#-------------------------------------------------------------------------------

# Acquire the sample ID from the line number (within ids.txt) specified by the current SLURM_ARRAY_TASK_ID
SAMPLE_ID=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SLURM_TMPDIR/ids.txt)

cd $OUT_DIR

# Run interleave-reads.py
time trim-low-abund.py \
  --variable-coverage \
  --cutoff 3 \
  --trim-at-coverage 18 \
  --gzip \
  -M 24e9 \
  --summary-info tsv \
  -o $OUT_DIR/${SAMPLE_ID}_001_P.ktrim.fastq.gz \
  $DATA_DIR/${SAMPLE_ID}_001_P.fastq.gz \
  >> $OUT_DIR/${SAMPLE_ID}_001_P.ktrim.log 2>> $OUT_DIR/${SAMPLE_ID}_001_P.ktrim.err
