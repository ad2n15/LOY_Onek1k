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
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 merge_ids_and_CA_summary.py
