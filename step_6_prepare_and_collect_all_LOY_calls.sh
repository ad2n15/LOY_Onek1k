# This bash script processes merged LOY (Loss Of Y chromosome) call files to prepare a final CSV file
# for downstream analysis in R. The steps include:
# 1. Concatenating all LOY call files while removing header lines containing "donor_id".
# 2. Counting the total number of LOY call entries after header removal.
# 3. Counting the number of unique LOY call entries.
# 4. Counting the number of unique combinations of the first two columns (e.g., donor and barcode).
# 5. Extracting LOY calls that appear only once (unique lines) and saving them to a CSV file.
# 6. Counting the number of unique LOY calls in the final CSV file.
# The resulting CSV file ("all_loy_calls_cranger_and_velo_with_barcodes2.csv") will be used for further analysis in R.

cat /mainfs/ddnb/Ahmed/Data/OneK1K/analysis/merged_LOY_calls_with_barcodes2/* | grep -v "donor_id" > t1
wc -l t1

sort t1 | uniq | wc -l

cut -f1,2 t1 | sort | uniq | wc -l

awk '{line[$0]++} END {for (i in line) if (line[i] == 1) print i}' t1 > all_loy_calls_cranger_and_velo_with_barcodes2.csv
wc -l all_loy_calls_cranger_and_velo_with_barcodes2.csv

