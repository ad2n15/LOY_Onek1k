#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=4:00:00
#SBATCH --array=1-14

# Load the sample variable
sample=$1

# Change to the appropriate directory
cd /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/

# Define paths to the necessary tools and files
subset_bam_linux=/mainfs/ddnb/Ahmed/images/subset-bam/subset-bam_linux
bam=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/possorted_genome_bam.bam
barcodes=/mainfs/ddnb/Ahmed/Data/OneK1K/scRNA_raw/${sample}_OneK1K_scRNA_Sample*_Individual_Barcodes.csv.gz
out_bam=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/out/

# Generate the individuals list if it doesn't already exist
if [ ! -f individuals_list.txt ]; then
  zcat $barcodes | cut -f1 -d, | sort | uniq | grep -v ID > individuals_list.txt
fi

# Get the individual for this array task
indv=$(awk -v id="$SLURM_ARRAY_TASK_ID" 'NR == id' individuals_list.txt)
echo $sample
echo $indv


module load conda

source activate /mainfs/ddnb/Ahmed/conda/r_env



data_dir=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/filtered_feature_bc_matrix/

all_barcode=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

zcat $all_barcode | grep -f barcode_${indv}.tsv  | gzip > /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/barcodes_${indv}_filtered.tsv.gz


barcode_indv=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/barcodes_${indv}_filtered.tsv.gz

output_dir=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/${sample}_cellranger/outs/filtered_feature_bc_matrix_${indv}/
mkdir $output_dir

Rscript /mainfs/ddnb/Ahmed/call_somatic_variant_Onek1k_2/filter_expresion_data_by_individual.R $data_dir $barcode_indv $output_dir

