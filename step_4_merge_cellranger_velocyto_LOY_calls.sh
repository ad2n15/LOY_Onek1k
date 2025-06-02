# -----------------------------------------------------------------------------
# Script Name: step_4_merge_cellranger_velocyto_LOY_calls.sh
#
# Description:
#   This SLURM batch script merges Loss Of Y chromosome (LOY) calls from two
#   different sources: CellRanger and Velocyto. It utilizes an Apptainer
#   container to ensure a reproducible Python environment for running the
#   merge_cellranger_velocyto_LOY_calls.py script. The script takes two input
#   CSV files containing LOY status calls from Velocyto and CellRanger,
#   respectively, and produces a merged output in the specified directory.
#
# Usage:
#   Submit this script to a SLURM cluster. It requires access to the specified
#   data directories and the datascience-notebook Apptainer image.
#
# SLURM Resources:
#   - 1 node
#   - 150 GB memory
#   - 60 hours walltime
#
# Main Steps:
#   1. Load the Apptainer module.
#   2. Execute the merge script inside the container with appropriate bind mounts.
#   3. Merge LOY calls from CellRanger and Velocyto into a single output.
#
# Inputs:
#   --file1 : Path to Velocyto LOY calls CSV file.
#   --file2 : Path to CellRanger LOY calls CSV file.
#   --output_dir : Directory to store the merged LOY calls.
#
# Author: [Your Name or Team]
# -----------------------------------------------------------------------------
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --time=60:00:00



cat $1 | while read a b
do

file1=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/velocyto_LOY_calls2/${a}_subset_${b}_*_LOY_status_velocyto.csv
file2=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/cellranger_LOY_calls/filtered_feature_bc_matrix_${b}_LOY_status.csv
output_file=/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls/${a}_${b}_LOY_total_calls.csv



module load apptainer


apptainer exec -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/metacells_project_py:/mainfs/ddnb/Ahmed/metacells_project_py \
--bind /mainfs/ddnb/Ahmed/images/copyVAE:/mainfs/ddnb/Ahmed/images/copyVAE \
--bind /mainfs/ddnb/Ahmed/call_mCA_genotype/:/mainfs/ddnb/Ahmed/call_mCA_genotype \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K/analysis:/mainfs/ddnb/Ahmed/Data/OneK1K/analysis \
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 scripts_for_preparing_scRNA-seq_LOY_calls/merge_cellranger_velocyto_LOY_calls.py \
--file1 $file1 \
--file2 $file2  \
--output_file $output_file

done

