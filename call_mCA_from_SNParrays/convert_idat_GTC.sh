bpm_manifest_file="..."
egt_cluster_file="..."
bcftools +idat2gtc \
  --bpm $bpm_manifest_file \
  --egt $egt_cluster_file \
  --idats $path_to_idat_folder \
  --output $path_to_gtc_folder
