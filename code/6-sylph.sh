#!/bin/bash

#SBATCH -J sylph_array
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --array=1-29
#SBATCH --output=sylph_%A_task_%a.out

# Number of nodes and tasks per node allocated by SLURM
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

# Set up tools
sylph="/private/groups/corbettlab/bianca/software/sylph"
sylphdb="/private/groups/corbettlab/bianca/databases/sylph/gtdb-r220-c200-dbv1.syldb"

# Directory containing your FASTQ files and results
FASTQ_DIR="/private/groups/corbettlab/bianca/subglacial/fq-fastp-sga-trim"
RESULTS_DIR="/private/groups/corbettlab/bianca/subglacial/sylph"

# Construct the filename based on the task array ID
SAMPLE_ID=$SLURM_ARRAY_TASK_ID
FASTQ_FILE="S${SAMPLE_ID}_unmappedhpcc.fq"

# Navigate to the PhyloFlash software directory
cd $RESULTS_DIR

# Run PhyloFlash with merged reads option
$sylph profile $sylphdb ${FASTQ_DIR}/${FASTQ_FILE} -t 16 > ${RESULTS_DIR}/S${SAMPLE_ID}profile.tsv
