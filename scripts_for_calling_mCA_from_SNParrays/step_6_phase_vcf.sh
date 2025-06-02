: '
This script is designed to phase VCF files generated from SNP array data using SHAPEIT5 and related tools in a high-performance computing (HPC) environment managed by SLURM. The script prepares input files, sets up environment modules, and processes each chromosome in parallel. Key steps include:

- Loading required modules and setting environment variables for plugin paths.
- Defining file paths for manifest, cluster, reference, and annotation files.
- Extracting sample call rate and gender information from a TSV file.
- Setting up output directories and file prefixes.
- Iterating over all autosomes and chromosome X to:
  - Fill in allele count (AN, AC) tags in BCF files using bcftools.
  - Index the resulting BCF files.
  - Prepare chromosome-specific genetic maps.
  - Run SHAPEIT5 phasing for each chromosome using the prepared inputs.

The script is intended for use in large-scale genotyping projects, such as the OneK1K project, and assumes the presence of pre-processed input files and reference panels.
'
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

cat Onek1k-withclusterfile.tsv | cut -f1,21 | sed 's/.gtc//g' | sed 's/gtc/sample_id/g' > call_rate.txt
cat Onek1k-withclusterfile.tsv | cut -f1,22 | sed 's/.gtc//g' | sed 's/M/1/g' | sed 's/F/2/g' | sed 's/gtc/sample_id/g' > calculated_gender.txt




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
for chr in {1..22} X; do
  bcftools +/local/software/samtools/1.20/./libexec/bcftools/fill-tags.so \
  $dir/$pfx.$chr.bcf -Ob -o $dir/$pfx.chr$chr.AN.bcf -- -t AN,AC

  bcftools index --force $dir/$pfx.chr$chr.AN.bcf
  zcat $map | sed 's/^23/X/' | awk -v chr=$chr '$1==chr {print $2,$3,$4}' > $dir/genetic_map.chr$chr.txt
  $phase_common \
    --thread $thr \
    --input $dir/$pfx.chr$chr.AN.bcf \
    --reference $panel_pfx${chr}$panel_sfx.bcf \
    --map $dir/genetic_map.chr$chr.txt \
    --region $chr \
    --output $dir/$pfx.chr$chr.pgt.bcf
done

