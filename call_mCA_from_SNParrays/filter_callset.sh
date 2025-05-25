cd /mainfs/ddnb/Ahmed/call_mCA_genotype/results

pfx=OneK1K


awk -F "\t" 'NR==FNR && FNR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR==FNR && FNR>1 && ($(f["call_rate"])<.97 || $(f["baf_auto"])>.03) {xcl[$(f["sample_id"])]++}
  NR>FNR && FNR==1 {for (i=1; i<=NF; i++) g[$i] = i; print}
  NR>FNR && FNR>1 {gender=$(g["computed_gender"]); len=$(g["length"]); bdev=$(g["bdev"]);
  rel_cov=$(g["rel_cov"]); lod_baf_phase=$(g["lod_baf_phase"]); lod_baf_conc=$(g["lod_baf_conc"])}
  NR>FNR && FNR>1 && !($(g["sample_id"]) in xcl) && $(g["type"])!~"^CNP" &&
    ( $(g["chrom"])~"X" && gender=="M" || bdev<0.1 || $(g["n50_hets"])<2e5 || lod_baf_conc!="nan" && lod_baf_conc>10.0 ) &&
    ( $(g["bdev_se"])!="nan" || lod_baf_phase!="nan" && lod_baf_phase>10.0 ) &&
    ( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )' \
  $pfx.stats.tsv $pfx.calls.tsv > $pfx.calls.filtered.tsv
awk 'NR==FNR {x[$1"_"$3"_"$4"_"$5]++} NR>FNR && ($0~"^track" || $4"_"$1"_"$2"_"$3 in x)' \
  $pfx.calls.filtered.tsv $pfx.ucsc.bed > $pfx.ucsc.filtered.bed



awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR>1 && ($(f["call_rate"])<.97 || $(f["baf_auto"])>.03) {print $(f["sample_id"])}' $pfx.stats.tsv > $pfx.remove.lines



awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i} NR>1 && $(f["computed_gender"])=="M" && $(f["chrom"])~"X" &&
  $(f["length"])>2e6 && $(f["rel_cov"])<2.5 {print $(f["sample_id"])}' $pfx.calls.tsv > $pfx.Y_loss.lines


awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR>1 && ($(f["computed_gender"])=="F" || $(f["computed_gender"])=="K") && $(f["chrom"])~"X" &&
  $(f["length"])>1e8 && $(f["rel_cov"])<2.5 {print $(f["sample_id"])}' $pfx.calls.tsv > $pfx.X_loss.lines
awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR>1 && ($(f["computed_gender"])=="F" || $(f["computed_gender"])=="K") && $(f["chrom"])~"X" &&
  $(f["length"])>1e8 && $(f["bdev"])>.01 && $(f["rel_cov"])<2.5 {print $(f["sample_id"])}' $pfx.calls.tsv > $pfx.X_loss_high.lines
awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
  NR>1 && ($(f["computed_gender"])=="F" || $(f["computed_gender"])=="K") && $(f["chrom"])~"X" &&
  $(f["length"])>1e8 && $(f["bdev"])<=.01 && $(f["rel_cov"])<2.5 {print $(f["sample_id"])}' $pfx.calls.tsv > $pfx.X_loss_low.lines









awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i; print}
  NR>1 {len=$(f["length"]); bdev=$(f["bdev"]); rel_cov=$(f["rel_cov"])}
  NR>1 && $(f["type"])!~"^CNP" &&
  ( $(f["chrom"])~"X" && $(f["computed_gender"])=="M" || bdev<0.1 || $(f["n50_hets"])<2e5 ) &&
  ( $(f["bdev_se"])!="nan" || $(f["lod_baf_phase"])!="nan" && $(f["lod_baf_phase"]) > 10.0 ) &&
  ( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )' \
  $pfx.calls.tsv > $pfx.mca.calls.tsv






awk -F"\t" -v OFS="\t" 'NR==FNR && $3!="X" && $3!="chrX" {x[$1]++} NR>FNR && $1 in x {print $1}' $pfx.mca.calls.tsv $pfx.stats.tsv > $pfx.auto.lines
for chr in {1..12} {16..20}; do
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["p_arm"])=="T" && $(f["q_arm"])!="T" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y*y<x*x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}p_cnloh.lines
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["p_arm"])!="N" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y>x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}p_loss.lines
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["p_arm"])=="T" && $(f["q_arm"])=="T" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y<-x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}_gain.lines
done
for chr in {1..22}; do
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["p_arm"])!="T" && $(f["q_arm"])=="T" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y*y<x*x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}q_cnloh.lines
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["q_arm"])!="N" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y>x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}q_loss.lines
done
for chr in 13 14 15 21 22; do
  awk -F"\t" -v OFS="\t" -v chr=$chr 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
    NR>1 && ($(f["chrom"])==chr || $(f["chrom"])=="chr"chr) && $(f["p_arm"])=="C" && $(f["q_arm"])=="T" && $(f["rel_cov"])>1 {
    x=$(f["bdev"]); y=(1/($(f["rel_cov"])-1)-1)/2; if (y<-x) print $(f["sample_id"])}' \
    $pfx.mca.calls.tsv > $pfx.${chr}_gain.lines
done




cat $pfx.{{3,4,6,8,11,17,18}p_loss,{1,6,11,13,14,16,17,22}q_loss,13q_cnloh,{2,3,4,5,12,17,18,19}_gain}.lines | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]++} NR>FNR && $1 in x {print $1}' - $pfx.stats.tsv > $pfx.cll.lines
cat $pfx.{{5,12,20}q_loss,{9p,14q,22q}_cnloh,{1,8}_gain}.lines  | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]++} NR>FNR && $1 in x {print $1}' - $pfx.stats.tsv > $pfx.myeloid.lines
cat $pfx.{{1,8,10,17}p_loss,{1,6,7,10,11,13,14,15,22}q_loss,{1q,7q,9q,12q,13q,16p}_cnloh,{2,3,9,12,15,17,18,19,22}_gain}.lines  | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]++} NR>FNR && $1 in x {print $1}' - $pfx.stats.tsv > $pfx.lymphoid.lines




(echo -e "sample_id\tY_loss\tX_loss\tX_loss_high\tX_loss_low"
awk -F"\t" -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i} NR>1 {print $(f["sample_id"]),$(f["computed_gender"])}' $pfx.stats.tsv | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=1} NR>FNR {if ($2=="M") phe=0+x[$1]; else phe="NA"; print $0,phe}' $pfx.Y_loss.lines - | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=1} NR>FNR {if ($2=="F" || $2=="K") phe=0+x[$1]; else phe="NA"; print $0,phe}' $pfx.X_loss.lines - | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=1} NR>FNR {if ($2=="F" || $2=="K") phe=0+x[$1]; else phe="NA"; print $0,phe}' $pfx.X_loss_high.lines - | \
  awk -F"\t" -v OFS="\t" 'NR==FNR {x[$1]=1} NR>FNR {if ($2=="F" || $2=="K") phe=0+x[$1]; else phe="NA"; print $0,phe}' $pfx.X_loss_low.lines - | \
  cut -f1,3-) > $pfx.pheno.tsv




for type in auto cll myeloid lymphoid {{{1..12},{16..20}}p,{1..22}q}_{cnloh,loss} {1..22}_gain; do
  n=$(cat $pfx.$type.lines | wc -l);
  if [ "$n" -gt 0 ]; then
    awk -F"\t" -v OFS="\t" -v type=$type 'NR==FNR {x[$1]=1}
      NR>FNR {if (FNR==1) col=type; else col=0+x[$1]; print $0,col}' \
      $pfx.$type.lines $pfx.pheno.tsv | sponge $pfx.pheno.tsv
  fi
done



: '

bcftools annotate \
  --no-version \
  --output $dir/$pfx.imp.as.bcf \
  --output-type b \
  --columns FMT/AS \
  $dir/$pfx.imp.bcf \
  --annotations $dir/$pfx.as.bcf \
  --write-index




bcftools +extendFMT \
  --no-version -Ou \
  --format AS \
  --phase \
  --dist 500000 \
  --regions $reg \
  --samples $lst \
  $dir/$pfx.imp.as.bcf | \
bcftools +mochatools \
  --no-version \
  --output $dir/$pfx.bal.bcf \
  --output-type b \
  -- --summary AS \
  --test AS \
  --drop-genotypes \
  --write-index

fmt="%CHROM\t%POS\t%ID\t%AS{0}\t%AS{1}\t%binom_AS\n"
bcftools query \
  --include "binom_AS<1e-6" \
  --format "$fmt" \
  $dir/$pfx.as.bcf | \
  column -ts $'\t'

'
