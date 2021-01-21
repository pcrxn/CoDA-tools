#!/bin/bash
#SBATCH --account=def-vtai4
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=bowtie2-mged_array
#SBATCH --output=bowtie2-mged_array_%j.out
#SBATCH --error=bowtie2-mged_array_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lbrow55@uwo.ca

# This script can be run using the following command:
# `sbatch --array=1-9 bowtie2-mged-submit.sh`
# 9 is the number of lines in the ids.txt document and
# the number of sample IDs

#-------------------------------------------------------------------------------
# Initialize environment
#-------------------------------------------------------------------------------

# Activate conda environment
module --force unload StdEnv/2020
module load nixpkgs/16.09 miniconda3
source activate bowtie2

#-------------------------------------------------------------------------------
# State global variables
#-------------------------------------------------------------------------------

# Absolute path to reference database (bowtie2-build output)
REF="/home/brownl/scratch/MAC-HiSeq/mged/bowtie2/MGEs_FINAL_99perc_trim.btindex"

# Directory of sorted-bam-to-fold-coverage.sh
CODEDIR="/home/brownl/scratch/MAC-HiSeq/code"

# Directory of fastq.gz files to be mapped
DATADIR="/home/brownl/scratch/MAC-HiSeq/data/combined/fastq"

# Directory where output files and subdirectories will be created
OUTDIR="/home/brownl/scratch/MAC-HiSeq/mged/bowtie2"

# Number of threads to use for called scripts
THREADS=8

#-------------------------------------------------------------------------------
# Retrieve sample IDs of files to be mapped
#-------------------------------------------------------------------------------

cd $DATADIR

# Obtain the base file names (sample IDs) and write to a temporary file
ls -l *R1_001_P.fastq.gz | cut -d ' ' -f 9 | sed 's/\(.*\)\_R1_001_P.fastq.gz/\1/g' \
  > $SLURM_TMPDIR/ids.txt

#-------------------------------------------------------------------------------
# Map reads to the reference database
#-------------------------------------------------------------------------------

# Acquire the sample ID from the line number (within ids.txt) specified by the
# current SLURM_ARRAY_TASK_ID
SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SLURM_TMPDIR/ids.txt)

cd $OUTDIR

echo "Now mapping $SAMPLE to $REF..."

bowtie2 \
  --time \
  --threads $THREADS \
  -D 20 \
  -R 3 \
  -N 1 \
  -L 20 \
  -i S,1,0.5 \
  -x $REF \
  -1 $DATADIR/${SAMPLE}_R1_001_P.fastq.gz \
  -2 $DATADIR/${SAMPLE}_R2_001_P.fastq.gz \
  -S $SLURM_TMPDIR/${SAMPLE}_mged-mapped.sam

echo "Mapping completed"

#-------------------------------------------------------------------------------
# Convert SAM to BAM
#-------------------------------------------------------------------------------

# Create an output directory based on the sample ID
mkdir $OUTDIR/$SAMPLE

echo "Now converting SAM to BAM..."

time samtools view -bS \
  --threads $THREADS \
  $SLURM_TMPDIR/${SAMPLE}_mged-mapped.sam > $OUTDIR/$SAMPLE/${SAMPLE}_mged-mapped.bam

echo "SAM file created"

#-------------------------------------------------------------------------------
# Sort the BAM file
#-------------------------------------------------------------------------------

echo "Now sorting BAM file..."

time samtools sort $OUTDIR/${SAMPLE}/${SAMPLE}_mged-mapped.bam \
  -o $OUTDIR/${SAMPLE}/${SAMPLE}_mged-mapped_sorted.bam \
  --threads $THREADS

echo "BAM file sorted"

#-------------------------------------------------------------------------------
# Index the sorted BAM file
#-------------------------------------------------------------------------------

echo "Indexing the sorted BAM file..."

time samtools index -@ $THREADS \
  $OUTDIR/${SAMPLE}/${SAMPLE}_mged-mapped_sorted.bam

echo "Sorted BAM file indexed"

#-------------------------------------------------------------------------------
# Generate alignment summary statistics, report read depths, and create table of
# fold-coverages
#-------------------------------------------------------------------------------

$CODEDIR/sorted-bam-to-fold-coverage.sh \
  $OUTDIR/${SAMPLE}/${SAMPLE}_mged-mapped_sorted.bam \
  $OUTDIR/${SAMPLE}
