#!/bin/bash

#SBATCH -J kegg
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 ## this is usually num threads if you want communication b/w threads
#SBATCH --time=40:00:00
#SBATCH --array=0-30 # needs to be 0-30

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

source ~/.bashrc
source /private/groups/corbettlab/bianca/software/miniforge3/etc/profile.d/conda.sh
conda activate functional_annotate

# Define the common paths
idba_ud="/private/groups/corbettlab/bianca/software/idba/bin/idba_ud"
sample_list=($(ls /private/groups/corbettlab/bianca/subglacial/fq/all_unmapped_to_contams_merged | sed 's/_unmappedhpcc\.fq//g'))
SAMPLE_ID=${sample_list[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: $SAMPLE_ID"
in_fq="/private/groups/corbettlab/bianca/subglacial/fq/all_unmapped_to_contams_merged/${SAMPLE_ID}_unmappedhpcc.fq"
idba_out="/private/groups/corbettlab/bianca/subglacial/kegg/idba/$SAMPLE_ID"
mkdir -p $idba_out
cd $idba_out
in_fa="${idba_out}/input.fa"
seqtk seq -A $in_fq > $in_fa
# turn the fq into an fa for idba input

echo "doing the idba"
$idba_ud -r $in_fa -o $idba_out
contig="${idba_out}/contig.fa"
filtered_contig="${idba_out}/contig_filtered.fa"
contig_cds="${idba_out}/contig_CDS.fna"
seqtk seq -L 500 $contig > $filtered_contig
# keep only the contigs that are more than 500bp long

if [ ! -f "$filtered_contig" ] || [ ! -s "$filtered_contig" ]; then
    echo "Error: $filtered_contig is missing or empty. Exiting."
    echo "fail" > fail_$SAMPLE_ID.out
    exit 1
fi

# prodigal will optionally turn it into a protein sequence
# i read that keeping it as dna is better for ancient guys
prodigal -i $filtered_contig -a ${idba_out}/contig_CDS.faa -d $contig_cds -o ${idba_out}/genes.gbk -p meta
# this predicts open reading frames

echo "doing the emapper"

# annotate using the kegg and cog databases using eggNOG-mapper and using the diamond database
emapper.py --itype CDS -i $contig_cds -o ${idba_out}/eggnog_output --override --cpu 8
# -m diamond if it's slow


echo "doing the mapping or read recruitment"
# use bwa aln w/ params to remap reads against filtered contigs and fnas. compute damage on both.
# the eggnog annotations are on the CDS files, but i think they contain the contig it came from, and it feels fair enough to me to compute damage on the whole contig, as it's one chunk.
bwa index $filtered_contig
mapped_bam="${idba_out}/mapped_to_contig_filtered.bam"
bwa aln -l 1024 -n 0.01 -o 2 $filtered_contig $in_fq > aln_output.sai
bwa samse $filtered_contig aln_output.sai $in_fq | samtools view -q 25 -b | samtools sort -n > $mapped_bam

# now i want to make a fake lca file so i can run bamdam and get damage stats per contig
# the lca format goes
# READ ID : READ: NUMBER : NUMBER \t NUMBER : NAME : PHYLUM \t 1:root:no rank
# so i gotta loop through the sam file and do
mapped_sam="$idba_out/mapped_to_contig_filtered.sam"
samtools view $mapped_bam > $mapped_sam
fake_lca1="$idba_out/fake_filtered_contig.lca"
> $fake_lca1 # initialize it

total_lines=$(wc -l < "$mapped_sam")
echo "Total lines to make an lca for: $total_lines"
echo "Processing file and generating fake LCA..."
awk -v total_lines="$total_lines" '
BEGIN { OFS = "\t"; processed = 0 }
!/^@/ {
    processed++
    # Print progress for every 10% of the file
    if (processed % int(total_lines / 10) == 0) {
        printf("%.0f%% done...\n", (processed / total_lines) * 100) > "/dev/stderr"
    }

    read_id = $1                         # First column (read ID)
    read_seq = $10                       # 10th column (read sequence)
    read_length = length($10)            # Calculate the length of the sequence
    ref_name = $3                        # 3rd column (reference name)

    split(ref_name, ref_parts, "_")      # Split reference name by "_"
    ref_number = ref_parts[2] + 0        # Extract numeric part, ensuring it is treated as a number

    number = ref_number + 2              # Increment the reference number by 2

    # Construct the new line
    new_line = read_id ":" read_seq ":" read_length ":1"
    new_line = new_line OFS number ":" ref_name ":superkingdom"
    new_line = new_line OFS "1:root:no rank"

    print new_line
}
END {
    print "100% done." > "/dev/stderr"
}' "$mapped_sam" > "$fake_lca1"


echo "running bamdam once"
# now we have to run bamdam but we don't need shrink
bamdam="/private/groups/corbettlab/bianca/software/bamdam/bamdam"
$bamdam compute --in_bam mapped_to_contig_filtered.bam --in_lca fake_filtered_contig.lca --out_tsv bd_contig_filtered.tsv --out_subs bd_contig_filtered.subs --stranded ss --upto root

echo "Sample $SAMPLE_ID processing complete."



