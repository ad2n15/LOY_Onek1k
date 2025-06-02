: '
This script converts Illumina GTC genotype files to a sorted and normalized BCF file using bcftools and the gtc2vcf plugin.
It is intended to be run as a SLURM batch job.

Description:
- Loads required modules and sets environment variables for bcftools plugins.
- Defines file paths for manifest, cluster, GTC input, and reference genome.
- Runs bcftools +gtc2vcf to convert GTC files to VCF, using manifest and cluster files for accurate genotype calling.
- Adjusts clusters using the provided EGT cluster file.
- Outputs extra information to a TSV file.
- Pipes the VCF output to bcftools sort and then to bcftools norm for sorting, normalization, and indexing.
- Produces a compressed, indexed BCF file as final output.

Requirements:
- bcftools with gtc2vcf plugin installed and accessible.
- Properly formatted manifest, cluster, and GTC files.
- SLURM workload manager for job scheduling.
'
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

bcftools +gtc2vcf \
  --no-version -Ou \
  --adjust-clusters \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --gtcs $path_to_gtc_folder \
  --egt $egt_cluster_file \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index


