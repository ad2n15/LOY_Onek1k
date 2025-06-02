# -----------------------------------------------------------------------------
# Script Name: step_5_extract_genotypes.sh
#
# Description:
#   This SLURM batch script processes genotype data from SNP arrays to generate
#   a filtered and annotated VCF file. It loads necessary modules, sets up
#   environment variables and file paths, and then uses a series of bcftools
#   commands to:
#     - Exclude variants with multiple alternate alleles and those present in
#       an exclusion list.
#     - Filter out uncalled genotypes and variants with no alternate alleles.
#     - Remove unnecessary VCF fields, retaining only genotype (GT) information.
#     - Split the output VCF by chromosome using the scatter plugin.
#
#   The script is tailored for the OneK1K project and expects specific input
#   files and directory structures.
#
# SLURM Options:
#   --job-name=gtcvcf1      : Name of the SLURM job.
#   --nodes=1               : Number of nodes to use.
#   --time=6:00:00          : Maximum runtime for the job.
#
# Requirements:
#   - samtools/1.20 module
#   - bcftools with scatter plugin
#   - Properly set input file paths and output prefix
#
# Output:
#   - Chromosome-scattered, filtered, and annotated VCF files with the specified prefix.
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



vcf=$out_prefix.vcf

bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 $out_prefix.unphased.bcf $out_prefix.xcl.bcf | \
  bcftools view --no-version -Ou --min-ac 0 --exclude-uncalled | \
  bcftools annotate --no-version -Ou --remove ID,QUAL,INFO,^FMT/GT | \
  bcftools +/local/software/samtools/1.20/./libexec/bcftools/scatter.so --no-version -Ob --output ./ --scatter $(echo {{1..22},X} | tr ' ' ',') --prefix $out_prefix.

