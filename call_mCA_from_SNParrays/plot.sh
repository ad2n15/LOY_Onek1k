module load conda
source activate /mainfs/ddnb/Ahmed/conda/r_env/
module load samtools/1.20


dir=/mainfs/ddnb/Ahmed/call_mCA_genotype/results/
pfx=OneK1K
mocha_dir=/mainfs/ddnb/Ahmed/images/MoChA/
$mocha_dir/summary_plot.R \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.calls.tsv \
  --pdf $dir/${pfx}_summary.pdf

$mocha_dir/pileup_plot.R \
  --cytoband /mainfs/ddnb/Ahmed/images/MoChA/mocha.GRCh37/cytoBand.txt.gz \
  --stats $dir/$pfx.stats.tsv \
  --calls $dir/$pfx.calls.filtered.tsv \
  --pdf $dir/${pfx}_pileup.pdf



