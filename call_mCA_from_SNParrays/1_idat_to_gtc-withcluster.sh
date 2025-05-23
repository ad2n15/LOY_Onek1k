#!/bin/bash
#SBATCH --job-name=cellranger_job
#SBATCH --nodes=1



module load samtools/1.20
export BCFTOOLS_PLUGINS=/mainfs/ddnb/Ahmed/images/MoChA/plugins/

cd /mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_idat


#for file in *.gz; do
#    zcat "$file" > "${file%.gz}"
#done


cd ..

mkdir all_idat_unzipped
#cp all_idat/*idat all_idat_unzipped/

path_to_idat_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_idat_unzipped
path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_gtc


bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.bpm
egt_cluster_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1_ClusterFile.egt
bcftools +idat2gtc \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --idats $path_to_idat_folder \
  --output $path_to_gtc_folder
