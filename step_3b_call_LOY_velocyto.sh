# This bash script is designed to run a Python script (`run_detect_LOY_genes_velocyto.py`)
# for detecting Loss Of Y chromosome (LOY) from velocyto counts in a high-performance
# computing environment managed by SLURM. The script reads a list of input items from a
# file (provided as the first argument), and for each item, it launches an Apptainer
# container with the necessary directories bound. The Python script is executed inside
# the container, processing each input and saving results to the specified output directory.
# The script requests 1 node, 150GB memory, and a maximum runtime of 60 hours.
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00
#SBATCH --mem=150G


module load apptainer


output_dir="/mainfs/ddnb/Ahmed/Data/OneK1K/analysis/velocyto_LOY_calls2"

cat $1 | while read a
do

apptainer run -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/metacells_project_py:/mainfs/ddnb/Ahmed/metacells_project_py \
--bind /mainfs/ddnb/Ahmed/images:/mainfs/ddnb/Ahmed/images \
--bind /mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes:/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes \
--unsquash \
/mainfs/ddnb/Ahmed/images/datascience-notebook_latest.sif \
python3 scripts_for_calling_LOY_from_scRNA-seq/run_detect_LOY_genes_velocyto.py \
$a $output_dir

done
