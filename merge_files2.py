import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Merge two files based on original barcodes.")
parser.add_argument('--file1', type=str, required=True, help="Path to the first input file.")
parser.add_argument('--file2', type=str, required=True, help="Path to the second input file.")
parser.add_argument('--output_file', type=str, required=True, help="Path to save the output merged file.")

# Parse arguments
args = parser.parse_args()

# Load files
file1 = pd.read_csv(args.file1, header=None, names=["donor_id", "barcode_file1"])
file2 = pd.read_csv(args.file2)

# Merge files directly on barcodes, keeping original barcodes
merged_data = pd.merge(file1, file2, left_on="barcode_file1", right_on="barcode", how="inner", suffixes=('_file1', '_file2'))

# Save to output file
merged_data.to_csv(args.output_file, index=False)
print(f"Merged data saved to {args.output_file}")

