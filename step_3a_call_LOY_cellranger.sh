#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00
#SBATCH --mem=150G


module load apptainer

output_dir="/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/cellranger_LOY_calls"
#for dir in $(find /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/filtered_counts/ -mindepth 2 -maxdepth 2 -type d); do

cat $1 | while read dir;
do

apptainer run -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/metacells_project_py:/mainfs/ddnb/Ahmed/metacells_project_py \
--bind /mainfs/ddnb/Ahmed/images:/mainfs/ddnb/Ahmed/images \
--bind /mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes:/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes \
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 run_detect_LOY_genes_cellranger.py \
$dir \
/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
$output_dir



done
