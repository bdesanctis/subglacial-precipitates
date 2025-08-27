#!/bin/bash

#SBATCH -J mapsubg
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 ## this is usually num threads if you want communication b/w threads
#SBATCH --time=10:00:00
#SBATCH --array=30-36 # test first
#SBATCH --output=array_job_%A_task_%a.out

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

### Set up tools ###
samtools="/usr/bin/samtools"

# There will be 90 mapping jobs in total and we can do it in 3 arrays, each array will map a single sample
 # the sample are everything ending in .fq in sampledir!
sampledir="/private/groups/corbettlab/bianca/subglacial/raw/pilot" # fq-fastp-sga-trim"
outputdir="/private/groups/corbettlab/bianca/subglacial/bam"
hpcc="/private/groups/corbettlab/bianca/databases/contaminants/HumanPigChickenCow.fna"

mkdir -p $outputdir
cd $outputdir

sample_index=$SLURM_ARRAY_TASK_ID    # $(( (SLURM_ARRAY_TASK_ID - 1) / (15 / refs_per_job) + 1 ))

# Get the sample file
sample=$(ls $sampledir/S${sample_index}.cleantrim.fq | head -n 1)
sample_id=$(basename $sample .cleantrim.fq)

# First map against possible contaminants and remove anything that maps
mapped_hpcc_bam="${outputdir}/S${sample_index}_mappedhpcc.bam"
unmapped_hpcc_fq="${sampledir}/S${sample_index}_unmappedhpcc.fq"
bowtie2 -p 16 -k 1 -x $hpcc -U $sample | $samtools view -bS - > $mapped_hpcc_bam

# Extract aligned reads
$samtools view -b -f 4 $mapped_hpcc_bam | $samtools fastq - > $unmapped_hpcc_fq
rm $mapped_hpcc_bam  # clean up
