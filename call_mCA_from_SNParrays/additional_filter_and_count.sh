cat mca_joined_numbers_2.txt | awk '{if ($1 !="155" && $1 !="155_2") {print }}'  | wc -l
