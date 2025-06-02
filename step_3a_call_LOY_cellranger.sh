# -----------------------------------------------------------------------------
# Script Name: step_3a_call_LOY_cellranger.sh
#
# Description:
#   This Bash script is designed to run a Python script for detecting Loss Of Y 
#   chromosome (LOY) events from CellRanger filtered count data. It is intended 
#   to be executed as a SLURM batch job on a high-performance computing cluster.
#   The script reads a list of directories (each containing filtered CellRanger 
#   count data) from an input file, and for each directory, it launches an 
#   Apptainer (Singularity) container to run the LOY detection Python script.
#
# Usage:
#   sbatch step_3a_call_LOY_cellranger.sh <input_directories.txt>
#
# Arguments:
#   <input_directories.txt> : A text file containing a list of directories, 
#                             each on a separate line, with filtered CellRanger 
#                             count data to process.
#
# SLURM Directives:
#   --nodes=1         : Request 1 compute node.
#   --time=60:00:00   : Set a maximum runtime of 60 hours.
#   --mem=150G        : Allocate 150 GB of memory.
#
# Main Steps:
#   1. Load the Apptainer module.
#   2. Set the output directory for LOY call results.
#   3. For each directory listed in the input file:
#      - Run the Python LOY detection script inside an Apptainer container,
#        binding necessary data and reference directories.
#      - Pass the directory, GTF reference, and output directory as arguments.
#
# Requirements:
#   - Apptainer (Singularity) installed and available as a module.
#   - Python script 'run_detect_LOY_genes_cellranger.py' available inside the 
#     container image.
#   - Properly structured input directories and reference files.
#
# Author: [Your Name]
# Date: [Date]
# -----------------------------------------------------------------------------
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
python3 scripts_for_calling_LOY_from_scRNA-seq/run_detect_LOY_genes_cellranger.py \
$dir \
/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
$output_dir



done
