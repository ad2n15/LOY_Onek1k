#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --time=60:00:00




module load apptainer


apptainer exec -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/metacells_project_py:/mainfs/ddnb/Ahmed/metacells_project_py \
--bind /mainfs/ddnb/Ahmed/images/copyVAE:/mainfs/ddnb/Ahmed/images/copyVAE \
--bind /mainfs/ddnb/Ahmed/call_mCA_genotype/:/mainfs/ddnb/Ahmed/call_mCA_genotype \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K/analysis:/mainfs/ddnb/Ahmed/Data/OneK1K/analysis \
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 merge_cellranger_velocyto_LOY_calls.py \
--file1 /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/velocyto_LOY_calls/GSM5899875_subset_669_670_60PJD_LOY_status_velocyto.csv \
--file2 /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/cellranger_LOY_calls/filtered_feature_bc_matrix_669_670_LOY_status.csv \
--output_dir /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls
