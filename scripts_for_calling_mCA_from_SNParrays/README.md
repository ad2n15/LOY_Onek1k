# mCA Calling Pipeline

This repository provides a step-by-step pipeline for calling mosaic chromosomal alterations (mCA) from SNP array data using the [MOCHA](https://github.com/freeseek/mocha) tool.

> **Note:** The provided scripts were originally developed and tested on the Iridis-5 HPC at the University of Southampton. You may need to modify paths, module loading commands, and job submission scripts to adapt the pipeline for use on other computing environments.

## Prerequisites

- [MOCHA](https://github.com/freeseek/mocha)
- PLINK
- bcftools
- R (for downstream analysis)
- SNP array data in appropriate format

## Pipeline Steps

The pipeline consists of the following scripts:

### Step 1: Convert IDAT to GTC
- `step_1_idat_to_gtc.sh`
- Converts raw IDAT files to GTC format.

### Step 2: Convert GTC to VCF
- `step_2_gtc_to_vcf.sh`
- Converts GTC files to VCF format for downstream processing.

### Step 3: Annotate VCF
- `step_3_annotate_vcf.sh`
- Adds necessary annotations to the VCF files.

### Step 4: Filter VCF
- `step_4_filter.sh`
- Applies quality control filters to the annotated VCF files.

### Step 5: Extract Genotypes
- `step_5_extract_genotypes.sh`
- Extracts genotype information required for phasing.

### Step 6: Phase VCF
- `step_6_phase_vcf.sh`
- Phases the filtered VCF files using a phasing tool (e.g., Eagle or SHAPEIT).

### Step 7: Call mCA
- `step_7_call_mca.sh`
- Runs MOCHA on the phased VCF files to detect mCA events.

## Usage

1. Clone this repository.
2. follow [MOCHA](https://github.com/freeseek/mocha) to install all dependencies 
3. Adjust any HPC-specific commands (such as module loading or job submission) as needed for your environment.
4- run the script after adjustment step by step


## References

- [MOCHA GitHub](https://github.com/freeseek/mocha)
- [PLINK](https://www.cog-genomics.org/plink/)

---

For questions or issues, please open an issue in this repository.