#!/bin/bash

#SBATCH -J bamdamcontrols
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ## this is usually num threads if you want communication b/w threads
#SBATCH --time=6:00:00
#SBATCH --array=1-5

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS


# Define the common paths
BASE_DIR="/private/groups/corbettlab/bianca/subglacial"
LCA_DIR="$BASE_DIR/ngslca-gtdb"
BAM_DIR="$BASE_DIR/bam-gtdb"
LOG_DIR="$BASE_DIR/logs"
ngslca="/private/groups/corbettlab/bianca/software/ngsLCA/ngsLCA"  # Correct path to ngslca
bamdam="/private/groups/corbettlab/bianca/software/bamdam/bamdam"
BAMDAM_DIR="$BASE_DIR/bamdam-gtdb"

NAMES="/private/groups/corbettlab/bianca/databases/gtdb/names.dmp"
NODES="/private/groups/corbettlab/bianca/databases/gtdb/nodes.dmp"
ACC2TAX="/private/groups/corbettlab/bianca/databases/gtdb/gtdb.acc2tax"

# Get the current sample ID
SAMPLE_IDS=("S14" "S28" "S29" "S34" "S35")
SAMPLE_ID=${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID-1]}


# Define input and output file paths
IN_LCA="$LCA_DIR/${SAMPLE_ID}.lca"
IN_BAM="$BAM_DIR/${SAMPLE_ID}_sorted.bam"
OUT_LCA="$BAMDAM_DIR/${SAMPLE_ID}_filtered.lca"
OUT_BAM="$BAMDAM_DIR/${SAMPLE_ID}_filtered.bam"
OUT_STATS="$BAMDAM_DIR/${SAMPLE_ID}.tsv"
OUT_SUBS="$BAMDAM_DIR/${SAMPLE_ID}.subs.txt"

cd $BAMDAM_DIR
#rm -f $OUT_BAM
#rm -f $OUT_LCA

# $ngslca -fix-ncbi 0 \
#  -names $NAMES \
#  -nodes $NODES \
#  -acc2tax $ACC2TAX \
#  -simscorelow 0.95 \
#  -bam $IN_BAM \
#  -outnames $LCA_DIR/${SAMPLE_ID}

# Check if the input LCA and BAM files exist before running bamdam shrink
if [ -f "$IN_LCA" ] && [ -f "$IN_BAM" ]; then
  echo "Input LCA and BAM files exist. Running bamdam shrink..."
  $bamdam shrink \
    --in_lca $IN_LCA \
    --in_bam $IN_BAM \
    --out_lca $OUT_LCA \
    --out_bam $OUT_BAM \
    --stranded ss \
    --upto phylum \
    --mincount 10 \
    --minsim 0.95
else
  echo "Error: Input LCA or BAM file is missing. Skipping bamdam shrink."
fi

# Check if the output LCA and BAM files from shrink exist before running bamdam compute
if [ -f "$OUT_LCA" ] && [ -f "$OUT_BAM" ]; then
  echo "Output LCA and BAM files from shrink exist. Running bamdam compute..."
  $bamdam compute \
    --in_lca $OUT_LCA \
    --in_bam $OUT_BAM \
    --out_tsv $OUT_STATS \
    --out_subs $OUT_SUBS \
    --stranded ss \
    --upto phylum
else
  echo "Error: Output LCA or BAM file from shrink is missing. Skipping bamdam compute."
fi


# after they all run, then do
# cat S14.stats.txt S28.stats.txt S29.stats.txt | awk '{print $1'} | sort | uniq > control_tax_list.txt
# which gets you a list of tax node ids in the controls

