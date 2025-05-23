#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00



## the master file for the code

export dir=/mainfs/ddnb/Ahmed/Data/OneK1K/


sbatch --wait --job-name=cellranger_$1 --output=cellranger_$1.log script/scripts_for_count_and_demultiplex/run_cellranger.sh $1

sbatch --wait --job-name=demuliplex_$1 --output=demuliplex_$1.log script/scripts_for_count_and_demultiplex/demuliplex_bam_file.sh $1

sbatch --wait --job-name=subset_$1 --output=subset_$1.log script/scripts_for_count_and_demultiplex/subset_expression_by_indv.sh $1


