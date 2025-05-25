: '
This script is the master file for running a series of SLURM jobs of the OneK1K project.
It performs three main steps for a given pool ID ($1):

1. Runs Cell Ranger to process sequencing data for the specified pool.
2. Demultiplexes the resulting BAM file for the same pool.
3. Subsets the expression data by individual within the pool.

Arguments:
    $1 - The pool ID to process.

Each step is submitted as a separate SLURM job using sbatch, with job names and log files
that include the pool ID for easy tracking. The script waits for each job to finish before
submitting the next one.
'
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:00:00



## the master file for the code

export dir=/mainfs/ddnb/Ahmed/Data/OneK1K/


sbatch --wait --job-name=cellranger_$1 --output=cellranger_$1.log script/scripts_for_count_and_demultiplex/run_cellranger.sh $1

sbatch --wait --job-name=demuliplex_$1 --output=demuliplex_$1.log script/scripts_for_count_and_demultiplex/demuliplex_bam_file.sh $1

sbatch --wait --job-name=subset_$1 --output=subset_$1.log script/scripts_for_count_and_demultiplex/subset_expression_by_indv.sh $1


