: '
This script is designed to convert Illumina IDAT files to GTC format using the bcftools +idat2gtc plugin. 
It is intended to be run as a SLURM batch job and includes the following steps:

- Loads required modules (samtools) and sets the BCFTOOLS_PLUGINS environment variable.
- Navigates to the directory containing compressed IDAT files.
- (Commented out) Optionally unzips .gz files and copies IDAT files to a new directory.
- Sets paths for input IDAT files, output GTC files, and required manifest and cluster files.
- Runs the bcftools +idat2gtc plugin with the specified parameters to perform the conversion.

This script is useful for preprocessing genotyping array data as part of a larger bioinformatics workflow.
'
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

#mkdir all_idat_unzipped
#cp all_idat/*idat all_idat_unzipped/

path_to_idat_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_idat_unzipped
path_to_gtc_folder=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/all_gtc_verify


bpm_manifest_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1.bpm
egt_cluster_file=/mainfs/ddnb/Ahmed/Data/OneK1K/genotyping/GSA-24v2-0_A1_ClusterFile.egt
bcftools +idat2gtc \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --idats $path_to_idat_folder \
  --output $path_to_gtc_folder
