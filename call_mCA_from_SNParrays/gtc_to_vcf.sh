#!/bin/bash
#SBATCH --job-name=gtcvcf1
#SBATCH --nodes=1
#SBATCH --time=6:00:00



module load samtools/1.20
export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/

#cd /mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/




bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSAMD-24v2-0_20024620_A1.bpm
csv_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSAMD-24v2-0_20024620_A1.csv

#egt_cluster_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1_ClusterFile.egt

path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_gtc-GSMAD
ref="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/human_g1k_v37.fasta"
out_prefix="Onek1k-noclusterfile"
bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index
