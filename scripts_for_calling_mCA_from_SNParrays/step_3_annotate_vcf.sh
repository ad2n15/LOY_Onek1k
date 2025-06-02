# -----------------------------------------------------------------------------
# Script Name: step_3_annotate_vcf.sh
#
# Description:
#   This SLURM batch script is designed to process and annotate VCF files 
#   generated from SNP array data in a high-performance computing environment.
#   It loads required modules, sets up environment variables, and defines 
#   file paths for manifest and cluster files. The script then converts a BCF 
#   file to VCF format and annotates the resulting VCF file using bcftools, 
#   retaining only specific INFO and FORMAT fields relevant for downstream 
#   analysis (e.g., ALLELE_A, ALLELE_B, GC, GT, BAF, LRR).
#
# Usage:
#   sbatch step_3_annotate_vcf.sh
#
# Requirements:
#   - SLURM workload manager
#   - samtools/1.20 module
#   - bcftools with plugins
#   - Properly set file paths for input and reference files
#
# Output:
#   - Annotated, filtered BCF file with index: <out_prefix>.unphased.bcf
# -----------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=gtcvcf1
#SBATCH --nodes=1
#SBATCH --time=6:00:00



module load samtools/1.20
export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/

#cd /mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/




bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.bpm
csv_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.csv

egt_cluster_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1_ClusterFile.egt

path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_gtc_verify
ref="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/human_g1k_v37.fasta"
out_prefix="Onek1k-withclusterfile"


bcftools convert -O v -o $out_prefix.vcf $out_prefix.bcf

vcf=$out_prefix.vcf

bcftools annotate --no-version -o $out_prefix.unphased.bcf -Ob --write-index $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR

