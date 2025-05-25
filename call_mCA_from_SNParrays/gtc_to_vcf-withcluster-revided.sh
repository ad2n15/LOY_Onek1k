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

path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_gtc
ref="/mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/human_g1k_v37.fasta"
out_prefix="Onek1k-withclusterfile"

: '
bcftools +gtc2vcf \
  --no-version -Ou \
  --adjust-clusters \
  --bpm $bpm_manifest_file \
  --csv $csv_manifest_file \
  --egt $egt_cluster_file \
  --gtcs $path_to_gtc_folder \
  --fasta-ref $ref \
  --extra $out_prefix.tsv | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -o $out_prefix.bcf -Ob -c x -f $ref --write-index

bcftools convert -O v -o $out_prefix.vcf $out_prefix.bcf
'
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

cat Onek1k-withclusterfile.tsv | cut -f1,21 | sed 's/.gtc//g' | sed 's/gtc/sample_id/g' > call_rate.txt
cat Onek1k-withclusterfile.tsv | cut -f1,22 | sed 's/.gtc//g' | sed 's/M/1/g' | sed 's/F/2/g' | sed 's/gtc/sample_id/g' > calculated_gender.txt




vcf=$out_prefix.vcf # input VCF file with phased GT, LRR, and BAF
pfx="OneK1K" # output prefix
thr=$SLURM_CPUS_PER_TASK # number of threads to use
crt="call_rate.txt" # tab delimited file with call rate information (first column sample ID, second column call rate)
sex="calculated_gender.txt" # tab delimited file with computed gender information (first column sample ID, second column gender: 1=male; 2=female)
xcl="..." # VCF file with additional list of variants to exclude (optional)
ped="..." # pedigree file to use if parent child duos are present
dir="/mainfs/ddnb/Ahmed/call_mCA_genotype/results" # directory where output files will be generated
mkdir -p $dir



#bcftools annotate --no-version -o $dir/$pfx.unphased.bcf -Ob --write-index $vcf \
#  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR




#awk -F"\t" '$2<.97 {print $1}' $crt > samples_xcl_list.txt; \
#echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
#  bcftools annotate --no-version -Ou -a $dup -c CHROM,FROM,TO,JK -h /dev/stdin $dir/$pfx.unphased.bcf | \
#  bcftools view --no-version -Ou -S ^samples_xcl_list.txt | \
#  bcftools +/local/software/samtools/1.20/./libexec/bcftools/fill-tags.so --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
#  bcftools view --no-version -Ou -G | \
#  bcftools annotate --no-version -o $dir/$pfx.xcl.bcf -Ob --write-index \
#    -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
#    -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING
#/bin/rm samples_xcl_list.txt


#bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 $dir/$pfx.unphased.bcf $dir/$pfx.xcl.bcf | \
#  bcftools view --no-version -Ou --min-ac 0 --exclude-uncalled | \
#  bcftools annotate --no-version -Ou --remove ID,QUAL,INFO,^FMT/GT | \
#  bcftools +/local/software/samtools/1.20/./libexec/bcftools/scatter.so --no-version -Ob --output $dir --scatter $(echo {{1..22},X} | tr ' ' ',') --prefix $pfx.chr


#phase_common=/mainfs/ddnb/Ahmed/images/shapeit5/phase_common_static
#for chr in {1..22} X; do
#  bcftools +/local/software/samtools/1.20/./libexec/bcftools/fill-tags.so \
#  $dir/$pfx.chr$chr.bcf -Ob -o $dir/$pfx.chr$chr.AN.bcf -- -t AN,AC

#  bcftools index --force $dir/$pfx.chr$chr.AN.bcf
#  zcat $map | sed 's/^23/X/' | awk -v chr=$chr '$1==chr {print $2,$3,$4}' > $dir/genetic_map.chr$chr.txt
#  $phase_common \
#    --thread $thr \
#    --input $dir/$pfx.chr$chr.AN.bcf \
#    --reference $panel_pfx${chr}$panel_sfx.bcf \
#    --map $dir/genetic_map.chr$chr.txt \
#    --region $chr \
#    --output $dir/$pfx.chr$chr.pgt.bcf
#done


#bcftools concat --no-version -o $dir/$pfx.pgt.bcf -Ob --write-index $dir/$pfx.chr{{1..22},X}.pgt.bcf

bcftools annotate --no-version -o $dir/$pfx.bcf -Ob --annotations $dir/$pfx.pgt.bcf --columns -FMT/GT --write-index $dir/$pfx.unphased.bcf
###########


cut -f2 call_rate.txt | paste calculated_gender.txt - > tsv.txt

#pfx="..." # output prefix
tsv=/mainfs/ddnb/Ahmed/call_mCA_genotype/tsv.txt # file with sample statistics (sample_id, computed_gender, call_rate)
lst="..." # file with list of samples to analyze for asymmetries (e.g. samples with 1p CN-LOH)
#cnp="..." # file with list of regions to genotype in BED format
#mhc_reg="..." # MHC region to skip
#kir_reg="..." # KIR region to skip





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
