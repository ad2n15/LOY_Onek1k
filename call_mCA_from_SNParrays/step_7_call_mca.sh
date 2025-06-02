#!/bin/bash
#SBATCH --job-name=gtcvcf1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
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

##########





vcf=$out_prefix.vcf # input VCF file with phased GT, LRR, and BAF
pfx="Onek1k-withclusterfile" # output prefix
thr=$SLURM_CPUS_PER_TASK # number of threads to use
crt="call_rate.txt" # tab delimited file with call rate information (first column sample ID, second column call rate)
sex="calculated_gender.txt" # tab delimited file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="/mainfs/ddnb/Ahmed/call_mCA_genotype_3" # directory where output files will be generated
mkdir -p $dir


phase_common=/mainfs/ddnb/Ahmed/images/shapeit5/phase_common_static



bcftools concat --no-version -o $dir/$pfx.pgt.bcf -Ob --write-index $dir/$pfx.chr{{1..22},X}.pgt.bcf

bcftools annotate --no-version -o $dir/$pfx.bcf -Ob --annotations $dir/$pfx.pgt.bcf --columns -FMT/GT --write-index $dir/$pfx.unphased.bcf
###########


cut -f2 call_rate.txt | paste calculated_gender.txt - > tsv.txt

tsv=/mainfs/ddnb/Ahmed/call_mCA_genotype-3/tsv.txt # file with sample statistics (sample_id, computed_gender, call_rate)





bcftools +mocha \
  --genome $assembly \
  --no-version \
  --output $dir/$pfx.as.bcf \
  --output-type b \
  --variants ^$dir/$pfx.xcl.bcf \
  --calls $dir/$pfx.calls.tsv \
  --stats $dir/$pfx.stats.tsv \
  --ucsc-bed $dir/$pfx.ucsc.bed \
  --write-index \
  --cnp $cnp \
  --mhc $mhc_reg \
  --kir $kir_reg \
  $dir/$pfx.bcf

