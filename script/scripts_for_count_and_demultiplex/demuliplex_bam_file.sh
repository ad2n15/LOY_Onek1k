#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=12:00:00
#SBATCH --array=1-14

# Load the sample variable
sample=$1

# Change to the appropriate directory
cd /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/

# Define paths to the necessary tools and files
subset_bam_linux=/mainfs/ddnb/Ahmed/images/subset-bam/subset-bam_linux
bam=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/possorted_genome_bam.bam
barcodes=/mainfs/ddnb/Ahmed/Data/OneK1K/scRNA_raw/${sample}_OneK1K_scRNA_Sample*_Individual_Barcodes.csv.gz
out_bam=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/

# Generate the individuals list if it doesn't already exist
if [ ! -f individuals_list.txt ]; then
  zcat $barcodes | cut -f1 -d, | sort | uniq | grep -v ID > individuals_list.txt
fi

# Get the individual for this array task
indv=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR == id' individuals_list.txt)
echo $sample
echo $indv
# Filter the barcodes for this individual
zcat $barcodes | grep -w $indv  | cut -f2 -d, > barcode_${indv}.tsv

# Run the subset-bam tool
${subset_bam_linux} --bam $bam \
  --cell-barcodes barcode_${indv}.tsv \
  --out-bam ${out_bam}/${sample}_subset_${indv}.bam \
  --cores $SLURM_CPUS_PER_TASK

