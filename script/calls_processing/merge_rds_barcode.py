import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Merge two files based on standardized barcodes.")
parser.add_argument('--file1', type=str, required=True, help="Path to the first input file.")
parser.add_argument('--file2', type=str, required=True, help="Path to the second input file.")
parser.add_argument('--output_file', type=str, required=True, help="Path to save the output merged file.")

# Parse arguments
args = parser.parse_args()

# Load files
file1 = pd.read_csv(args.file1, header=None, names=["donor_id", "barcode_file1"])
file2 = pd.read_csv(args.file2)

# Standardize barcode formats by removing suffixes
file1['barcode_standardized'] = file1['barcode_file1'].str.replace(r'-\d+$', '', regex=True)
file2['barcode_standardized'] = file2['barcode'].str.replace(r'-\d+$', '', regex=True)

# Merge files based on standardized barcodes, keeping original barcodes
merged_data = pd.merge(file1, file2, on="barcode_standardized", how="inner", suffixes=('_file1', '_file2'))

# Save to output file
merged_data.to_csv(args.output_file, index=False)
print(f"Merged data saved to {args.output_file}")

