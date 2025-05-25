module load samtools/1.20
export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/

vcf=Onek1k-withclusterfile.vcf
pfx=Onek1k-withclusterfile
dir=/mainfs/ddnb/Ahmed/call_mCA_genotype/

bcftools annotate --no-version -o $dir/$pfx.unphased.bcf -Ob --write-index $vcf \
  -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR
