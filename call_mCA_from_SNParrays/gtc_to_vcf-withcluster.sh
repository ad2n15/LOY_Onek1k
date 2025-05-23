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

vcf=$out_prefix.vcf

bcftools annotate --no-version -o $out_prefix.unphased.bcf -Ob --write-index $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR

'


#vcf=$out_prefix.vcf


#cat Onek1k-withclusterfile.tsv | cut -f1,21 | awk -F"\t" '$2<.97 {print $1}' | cut -f1 -d"."  > samples_xcl_list.txt; \

#  bcftools view --no-version -Ou -S ^samples_xcl_list.txt  $out_prefix.unphased.bcf | \
#  bcftools +/local/software/samtools/1.20/./libexec/bcftools/fill-tags.so --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
#  bcftools view --no-version -Ou -G | \
#  bcftools annotate --no-version -o $out_prefix.xcl.bcf -Ob --write-index \
#    -i 'FILTER!="." && FILTER!="PASS"  || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-.97' \
#    -x ^INFO/ExcHet,^INFO/F_MISSING
#rm samples_xcl_list.txt




#bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 $out_prefix.unphased.bcf $out_prefix.xcl.bcf | \
#  bcftools view --no-version -Ou --min-ac 0 --exclude-uncalled | \
#  bcftools annotate --no-version -Ou --remove ID,QUAL,INFO,^FMT/GT | \
#  bcftools +/local/software/samtools/1.20/./libexec/bcftools/scatter.so --no-version -Ob --output ./ --scatter $(echo {{1..22},X} | tr ' ' ',') --prefix $out_prefix.



$panel_pfx
$panel_sfx


for chr in {1..22} X; do
  bcftools index --force $out_prefix.chr$chr.bcf
  zcat $map | sed 's/^23/X/' | awk -v chr=$chr '$1==chr {print $2,$3,$4}' > genetic_map.chr$chr.txt
  phase_common \
    --thread $thr \
    --input $out_prefix.chr$chr.bcf \
    --reference $panel_pfx${chr}$panel_sfx.bcf \
    --map $dir/genetic_map.chr$chr.txt \
    --region $chr \
    --output $dir/$pfx.chr$chr.pgt.bcf
done
