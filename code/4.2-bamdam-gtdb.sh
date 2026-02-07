
#!/bin/bash

#SBATCH -J bamdam
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ## this is usually num threads if you want communication b/w threads
#SBATCH --time=1:00:00
#SBATCH --array=0-25 # non control samples

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
sample_list=("AIS1" "BV1" "BV2" "BV3" "EM1" "EM2" "EM3" "EM4" "LG1" "LG2" "LG3" "LG4" "LG5" "LG6" "LG7" "MA1" "MA2" "MBL1" "MBL2" "MV1" "MV2" "NC1" "NC2" "NC3" "NC4" "NC5" "NC6" "NC7" "PC1" "PM1" "PM2" "RM1" "RM2")
SAMPLE_ID=${sample_list[$SLURM_ARRAY_TASK_ID]}

# Define input and output file paths
IN_LCA="$LCA_DIR/${SAMPLE_ID}.lca"
IN_BAM="$BAM_DIR/${SAMPLE_ID}_sorted.bam"

OUT_LCA="$BAMDAM_DIR/${SAMPLE_ID}_filtered.lca"
OUT_BAM="$BAMDAM_DIR/${SAMPLE_ID}_filtered.bam"
OUT_STATS="$BAMDAM_DIR/${SAMPLE_ID}.tsv"
OUT_SUBS="$BAMDAM_DIR/${SAMPLE_ID}.subs.txt"

OUT_LCA_C="$BAMDAM_DIR/${SAMPLE_ID}_filtered_includecontroltax.lca"
OUT_BAM_C="$BAMDAM_DIR/${SAMPLE_ID}_filtered_includecontroltax.bam"
OUT_STATS_C="$BAMDAM_DIR/${SAMPLE_ID}_includecontroltax.tsv"
OUT_SUBS_C="$BAMDAM_DIR/${SAMPLE_ID}.subs_includecontroltax.txt"

$ngslca -fix-ncbi 0 \
  -names $NAMES \
  -nodes $NODES \
  -acc2tax $ACC2TAX \
  -simscorelow 0.95 \
  -bam $IN_BAM \
  -outnames $LCA_DIR/${SAMPLE_ID}

# Check if the input LCA and BAM files exist before running bamdam shrink
if [ -f "$IN_LCA" ] && [ -f "$IN_BAM" ]; then
  echo "Input LCA and BAM files exist. Running bamdam shrink..."
  $bamdam shrink \
    --in_lca $IN_LCA \
    --in_bam $IN_BAM \
    --out_lca $OUT_LCA_C \
    --out_bam $OUT_BAM_C \
    --stranded ss \
    --upto phylum \
    --annotate_pmd \
    --minsim .95
else
  echo "Error: Input LCA or BAM file is missing. Skipping bamdam shrink."
fi

# Check if the output LCA and BAM files from shrink exist before running bamdam compute
if [ -f "$OUT_LCA_C" ] && [ -f "$OUT_BAM_C" ]; then
  echo "Output LCA and BAM files from shrink exist. Running bamdam compute..."
  $bamdam compute \
    --in_lca $OUT_LCA_C \
    --in_bam $OUT_BAM_C \
    --out_tsv $OUT_STATS_C \
    --out_subs $OUT_SUBS_C \
    --stranded ss \
    --upto phylum
else
  echo "Error: Output LCA or BAM file from shrink is missing. Skipping bamdam compute."
fi


$bamdam shrink \
    --in_lca $IN_LCA \
    --in_bam $IN_BAM \
    --out_lca $OUT_LCA \
    --out_bam $OUT_BAM \
    --stranded ss \
    --upto phylum \
    --annotate_pmd \
    --minsim 0.95 \
    --exclude_keyword_file /private/groups/corbettlab/bianca/subglacial/bamdam-gtdb/control_tax_list.txt

$bamdam compute \
    --in_lca $OUT_LCA \
    --in_bam $OUT_BAM \
    --out_tsv $OUT_STATS \
    --out_subs $OUT_SUBS \
    --stranded ss \
    --upto phylum
