#!/bin/bash
#SBATCH --job-name=cellranger_job
#SBATCH --nodes=1



module load samtools/1.20
export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/

#cd /mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_idat




bpm_manifest_file=bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.bpm
csv_manifest_file=bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.csv

egt_cluster_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1_ClusterFile.egt
path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/testgtc
ref="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/human_g1k_v37.fasta"
out_prefix="Onek1k"
bcftools +gtc2vcf \
  --no-version -Ou \
  --adjust-clusters \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index
