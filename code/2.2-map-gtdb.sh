#!/bin/bash

#SBATCH -J mapsubg
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 ## this is usually num threads if you want communication b/w threads
#SBATCH --time=11:00:00
#SBATCH --array=30-36 # test first!!
#SBATCH --output=array_job_%A_task_%a.out

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

### Set up tools ###
samtools="/usr/bin/samtools"
ngslca="/private/groups/corbettlab/bianca/software/ngsLCA/ngsLCA"
metadmg="/private/groups/corbettlab/bianca/software/metaDMG-cpp/metaDMG-cpp"

# Set up directories
sampledir="/private/groups/corbettlab/bianca/subglacial/fq/pilot"
outputdir="/private/groups/corbettlab/bianca/subglacial/bam-gtdb"
refdir="/private/groups/corbettlab/bianca/databases/gtdb"

# Database files
names="/private/groups/corbettlab/bianca/databases/gtdb/names.dmp"
nodes="/private/groups/corbettlab/bianca/databases/gtdb/nodes.dmp"
acc2tax="/private/groups/corbettlab/bianca/databases/gtdb/gtdb.acc2tax"

mkdir -p $outputdir
cd $outputdir

# Get the sample files
total_samples=7
sample_index=$SLURM_ARRAY_TASK_ID
unmapped_hpcc_fq="${sampledir}/S${sample_index}_unmappedhpcc.fq"

for (( i=1; i<=15; i++ )); do
	ref="${refdir}/gtdb${i}.fna"
    	output_file="${outputdir}/S${sample_index}_gtdb${i}.bam"
	bowtie2 -p 16 -k 100 -N 1 -x $ref -U $unmapped_hpcc_fq --no-unal | $samtools view -bS > $output_file
done

# Process BAM files to remove unnecessary headers
for (( i=1; i<=15; i++ )); do
    ref="${refdir}/gtdb${i}.fna"
    ref_basename=$(basename $ref .fna)
    input_bam="${outputdir}/S${sample_index}_gtdb${i}.bam"
    output_bam="${outputdir}/S${sample_index}_gtdb${i}_cleaned.bam"
    sample_base=$(basename "$input_bam" .bam)

    ogheaders_tmp="${outputdir}/ogh_${sample_base}.tmp"
    mapped_headers_tmp="${outputdir}/mh_${sample_base}.tmp"
    tempsam="${outputdir}/tmp_${sample_base}.sam"

    # Get the original header
    $samtools view -H "$input_bam" > "$ogheaders_tmp"

    # Grab all the header lines you want to keep
    $samtools view "$input_bam" | cut -f 3 | sort | uniq > "$mapped_headers_tmp"
    grep '^@HD' "$ogheaders_tmp"  >> "$tempsam"
    grep '^@SQ' "$ogheaders_tmp" | grep -F -f "$mapped_headers_tmp" >> "$tempsam"
    grep '^@PG' "$ogheaders_tmp"  >> "$tempsam"

    # Add non-header lines
    $samtools view "$input_bam" >> "$tempsam"

    # Turn it into a BAM file
    $samtools view -b -o "$output_bam" "$tempsam"

    # Remove temporary files
    rm "$tempsam"
    rm "$mapped_headers_tmp"
    rm "$ogheaders_tmp"
done

# Merge bams
cleaned_bams=($(ls ${outputdir}/S${sample_index}_*cleaned.bam))
# Stop if things didn't work for whatever reason
if [ ${#cleaned_bams[@]} -ne 15 ]; then
    echo "Error: Expected 15 cleaned BAM files, but found ${#cleaned_bams[@]}. Exiting."
    exit 1
fi
merged_bam="${outputdir}/S${sample_index}_merged.bam"
$samtools merge --threads 15 -f "$merged_bam" "${cleaned_bams[@]}"

sorted_bam="${outputdir}/S${sample_index}_sorted.bam"
$samtools sort -n -@ 15 $merged_bam > $sorted_bam

# If things seem to have worked, remove intermediate files
if [ -s "$sorted_bam" ]; then
    rm "$merged_bam"
    rm "${cleaned_bams[@]}"
else
    echo "Error: Sorted BAM file does not exist or is empty. Exiting."
    exit 1
fi

# Run run lca analyses!

ngslca_outputdir="/private/groups/corbettlab/bianca/subglacial/ngslca-gtdb"
mkdir -p $ngslca_outputdir

# Run ngslca
cd $ngslca_outputdir
$ngslca -fix-ncbi 0 -names $names -nodes $nodes -acc2tax $acc2tax -simscorelow 0.95 -bam $sorted_bam -outnames S$sample_index


