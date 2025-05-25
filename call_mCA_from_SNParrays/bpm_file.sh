module load samtools/1.20

export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/


path_to_idat_folder=$1
bcftools +gtc2vcf -i -g $path_to_idat_folder
