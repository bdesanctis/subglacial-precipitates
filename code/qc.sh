#!/bin/bash

#SBATCH -J qc
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=76
###SBATCH --mem=120000
##SBATCH --cpus-per-task=1
#SBATCH --time=35:00:00
##SBATCH --mail-type=FAIL
##SBATCH --no-requeue
#SBATCH --array=1-31


SAMPLE_LIST="/rds/project/ew482/rds-ew482-geogenetics/bdd28/subglacial/code/configs/SE7878_samplelist.txt"
adap_list_path="/rds/project/ew482/rds-ew482-geogenetics/bdd28/subglacial/code/adapter_list.fa"
fastp="/rds/project/ew482/rds-ew482-geogenetics/bdd28/software/fastp"
sga="/rds/project/ew482/rds-ew482-geogenetics/bdd28/software/sga/src/SGA/sga"
outdir="/rds/project/ew482/rds-ew482-geogenetics/bdd28/subglacial/fq/after-fastp"
outdir2="/rds/project/ew482/rds-ew482-geogenetics/bdd28/subglacial/fq/fastp-then-sga"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)

# Extract the R1 file path
R1_fq=$(echo $line | cut -d ' ' -f 2)
temp=${R1_fq%_R1_001.fastq.gz}
lib=$(basename ${temp})
# Construct R2 file path
R2_fq="${temp}_R2_001.fastq.gz"

# Define output file paths using the base filename
merged_out="${outdir}/${lib}.ppm.fq"
html_report="${outdir}/${lib}.fastp.report.html"

# Run fastp
$fastp -i ${R1_fq} -I ${R2_fq} -m --merged_out ${merged_out} -V --adapter_fasta $adap_list_path -D --dup_calc_accuracy 5 -g -x --poly_g_min_len 5 --poly_x_min_len 5 -q 30 -e 25 -l 30 -y -c -p -h ${html_report} -w 50

# SGA outputs
prepfile="${outdir2}/${lib}.prep.fq"
cleanfile="${outdir2}/${lib}.clean.fq"
# Run SGA
cd $outdir2
$sga preprocess --dust-threshold=30 -m 30 ${merged_out} -o ${prepfile}
$sga index --algorithm=ropebwt -t 76 ${prepfile}
$sga filter --threads=76 --no-kmer-check ${prepfile} -o ${cleanfile}
# Check for errors in SGA filtering
if [ $? -ne 0 ]; then
    echo "Error with SGA filtering"
    exit 1
fi

# Run FastQC on the cleaned file
fastqc -noextract -o ${outdir2} ${cleanfile}
# Check for errors in FastQC
if [ $? -ne 0 ]; then
    echo "Error with FastQC"
    exit 1
fi
