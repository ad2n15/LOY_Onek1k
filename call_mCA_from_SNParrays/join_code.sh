#!/bin/bash

# Sort the files
sort -k1,1 mCA_joined_improved-2.txt > mCA_sorted.txt
sort -k1,1 Onk1k_id_pool_age_sex-2-sorted.txt > Onk1k_sorted.txt

# Perform inner join
join -1 1 -2 1 mCA_sorted.txt Onk1k_sorted.txt > inner_join.txt

# Find unique records
comm -23 mCA_sorted.txt Onk1k_sorted.txt > unique_mCA.txt
comm -13 mCA_sorted.txt Onk1k_sorted.txt > unique_Onk1k.txt

# Combine all records
cat inner_join.txt unique_mCA.txt unique_Onk1k.txt > full_outer_join.txt

# Clean up
rm mCA_sorted.txt Onk1k_sorted.txt inner_join.txt unique_mCA.txt unique_Onk1k.txt

