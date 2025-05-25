#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40 # Number of CPUs you want to use
#SBATCH --time=12:00:00

# Load necessary modules (if applicable)
# module load cellranger # Uncomment if you need to load Cell Ranger as a module


sample=$1
echo $sample
# Set the working directory
cd /mainfs/ddnb/Ahmed/Data/OneK1K/analysis

# Define the path to the Cell Ranger executable
CELLRANGER_PATH=/mainfs/ddnb/Ahmed/images/cellranger/cellranger-8.0.1/cellranger

# Define the sample IDs and the transcriptome reference
cd $sample
SAMPLE_IDS=$(ls *.gz | cut -f1 -d"_" | sort | uniq | tr '\n' ',' | paste -sd "," | sed 's/,$//g')
cd ..
TRANSCRIPTOME_REF=/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A
FASTQ_DIR=${sample}

# Create a unique output directory
OUTPUT_DIR=${sample}_cellranger

# Run Cell Ranger count with specified number of threads
$CELLRANGER_PATH count --id=${sample}_cellranger \
                       --sample=$SAMPLE_IDS \
                       --transcriptome=$TRANSCRIPTOME_REF \
                       --fastqs=$FASTQ_DIR \
                       --create-bam=true \
                       --localcores=$SLURM_CPUS_PER_TASK \
                       --localmem=150 # Adjust based on the memory allocated


rm $FASTQ_DIR/*fastq.gz
