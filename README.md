# LOY_OneK1K

This repository contains the complete pipeline used to identify and analyze mosaic loss of chromosome Y (LOY) in immune cells using single-cell RNA-seq data from the OneK1K study.

---

## 🔬 Project Overview

LOY is the most frequent somatic event in aging men and has been associated with increased risk for multiple age-related diseases. In this study, we utilized single-cell transcriptomic data and SNP array genotyping to identify LOY at cellular resolution and explore its effects on immune cell phenotypes and gene expression programs.

Our integrated pipeline spans data preprocessing, LOY calling, and data harmonization across multiple platforms (Cell Ranger, Velocyto, SNP arrays), culminating in biological insights visualized through an interactive Shiny app.

---

## 📁 Repository Structure

| Folder/File | Description |
|-------------|-------------|
| `call_mCA_from_SNParrays/` | Bash scripts and annotations for calling mosaic chromosomal alterations (mCA) from SNP array data. |
| `mapping_and_counting/` | Bash scripts and annotations for initial mapping and gene counting steps. |
| `script/` | Miscellaneous annotated Bash scripts used throughout the pipeline. |
| `filter_expression_data_by_individual.R` | R script used to filter gene expression data by individual. |
| `step_1_count_and_demultiplex_cellranger.sh` | Demultiplexing and read count processing using Cell Ranger. |
| `step_2_count_velocyto.sh` | Running Velocyto to count spliced/unspliced reads. |
| `step_3a_call_LOY_cellranger.sh` | Identify LOY cells using Cell Ranger output. |
| `step_3b_call_LOY_velocyto.sh` | Identify LOY cells using Velocyto output. |
| `step_4_merge_cellranger_velocyto_LOY_calls.sh` | Merge LOY calls from Cell Ranger and Velocyto. |
| `step_5_process_barcodes_to_match_rds_file.sh` | Align LOY calls with metadata from Seurat RDS objects. |
| `step_6_prepare_and_collect_all_LOY_calls.sh` | Aggregate final LOY calls for analysis. |
| `README.md` | This file. |
| `LOY_Onek1k/` | Summary folder containing core results and scripts. |

---

## 📊 Interactive Results Portal

Explore LOY patterns by donor, cell type, and gene expression using our Shiny app:

👉 https://lagreen.shinyapps.io/1k1kloy/

---


