#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --time=60:00:00


cat $1 | while read a b
do
grep -w $b onek1k_barcode_table_by_donor_id.csv > ${b}_specific_barcod

file2=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls/${a}_${b}_LOY_total_calls.csv
output=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls_with_barcodes2/${b}_LOY_total_calls_withbarcode.csv


module load apptainer

apptainer exec -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/metacells_project_py:/mainfs/ddnb/Ahmed/metacells_project_py \
--bind /mainfs/ddnb/Ahmed/images/copyVAE:/mainfs/ddnb/Ahmed/images/copyVAE \
--bind /mainfs/ddnb/Ahmed/call_mCA_genotype/:/mainfs/ddnb/Ahmed/call_mCA_genotype \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K/analysis:/mainfs/ddnb/Ahmed/Data/OneK1K/analysis \
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 merge_files.py --file1 ${b}_specific_barcod --file2 $file2 --output_file  $output




done
