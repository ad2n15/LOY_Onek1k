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






ref="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/human_g1k_v37.fasta"
mhc_reg="6:27486711-33448264"
kir_reg="19:54574747-55504099"
map="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/genetic_map_hg19_withX.txt.gz"
panel_pfx="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/ALL.chr"
panel_sfx=".phase3_integrated.20130502.genotypes"
assembly="GRCh37"
cnp="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/cnps.bed"
dup="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/segdups.bed.gz"
cyto="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/cytoBand.txt.gz"





crt="call_rate.txt"
dir="/mainfs/ddnb/Ahmed/call_mCA_genotype_3"

vcf=$out_prefix.vcf






awk -F"\t" '$2<.97 {print $1}' $crt > samples_xcl_list.txt; \
echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin $dir/$out_prefix.unphased.bcf | \
  bcftools view --no-version -Ou -S ^samples_xcl_list.txt | \
  bcftools +/local/software/samtools/1.20/./libexec/bcftools/fill-tags.so --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
  bcftools view --no-version -Ou -G | \
  bcftools annotate --no-version -o $dir/${out_prefix}.xcl.bcf -Ob --write-index \
    -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
    -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING
#/bin/rm samples_xcl_list.txt

