

module load apptainer


export run_verocyto="apptainer exec -H $PWD \
--bind /mainfs/ddnb/Ahmed/Data/OneK1K:/mainfs/ddnb/Ahmed/Data/OneK1K \
--bind /mainfs/ddnb/Ahmed/LOY_project:/mainfs/ddnb/Ahmed/LOY_project \
--bind /mainfs/ddnb/Ahmed/LOY_project/.local/bin:/mainfs/ddnb/Ahmed/LOY_project/.local/bin \
--bind /mainfs/ddnb/Ahmed/images/ref_mask:/mainfs/ddnb/Ahmed/images/ref_mask \
--bind /mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes:/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes \
--unsquash \
/mainfs/ddnb/Ahmed/images/velocyto_0.17.sif \
velocyto"

filtered_barcodes=$1
outdir=$2
bamfile=$3
repeat_msk="/mainfs/ddnb/Ahmed/images/ref_mask/repeat_msk.gtf"
GTFFILE="/mainfs/ddnb/Ahmed/reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
cpu=$4

$run_verocyto run -b $filtered_barcodes -o $outdir -m $repeat_msk -@ $cpu $bamfile $GTFFILE
