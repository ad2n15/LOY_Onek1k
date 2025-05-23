cat /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls_with_barcodes2/* | grep -v "donor_id" > t1
wc -l t1

sort t1 | uniq | wc -l

cut -f1,2 t1 | sort | uniq | wc -l

awk '{line[$0]++} END {for (i in line) if (line[i] == 1) print i}' t1 > all_loy_calls_cranger_and_velo_with_barcodes2.csv
wc -l all_loy_calls_cranger_and_velo_with_barcodes2.csv

