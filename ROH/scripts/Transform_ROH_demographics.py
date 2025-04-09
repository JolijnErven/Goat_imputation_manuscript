#!/usr/bin/env python

"""
This script processes an input file (output file from get_ROH_demographics.py) 
and seperates ROH sum and count and transposes data
"""

import argparse
import pandas as pd

# Argument Parser
parser = argparse.ArgumentParser(description='Create input for mapping for Dublin samples')
parser.add_argument("-i", "--input", required=True, help="Input txt file")
parser.add_argument("-o", "--out", required=True, help="Output file prefix")

def transform(input_lines, output_prefix):
    """Processes input data and generates transposed output files for ROH sum and count."""
    plot_data = []
    
    for line in input_lines:
        if not line.startswith("ID"):  # Skip header line
            values = line.strip().split('\t')
            plot_data.append([values[0]] + [float(x) for x in values[1:]])
    
    # Create DataFrame
    columns = ['Samples', 'roh count 0.5-1.0 Mb', '0.5-1.0 Mb', 'froh 0.5-1.0 Mb',
               'roh count 1.0-2.0 Mb', '1.0-2.0 Mb', 'froh 1.0-2.0 Mb',
               'roh count 2.0-4.0 Mb', '2.0-4.0 Mb', 'froh 2.0-4.0 Mb',
               'roh count 4.0-8.0 Mb', '4.0-8.0 Mb', 'froh 4.0-8.0 Mb',
               'roh count 8.0-16.0 Mb', '8.0-16.0 Mb', 'froh 8.0-16.0 Mb',
               'roh count >16.0 Mb', '>16.0 Mb', 'froh >16.0 Mb', 'sum_ROH', 'Froh']
    df = pd.DataFrame(plot_data, columns=columns)
    
    # Convert all numeric columns to proper types
    numeric_cols = columns[1:]  # Exclude 'Samples'
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)
    
    # Melt DataFrame for transformation
    trans_df = pd.melt(df, id_vars=['Samples'], value_vars=numeric_cols)
    
    # Filter specific categories
    all_count_ROH = trans_df[trans_df["variable"].str.contains('roh count')]
    all_sum_ROH = trans_df[trans_df["variable"].isin(['0.5-1.0 Mb','1.0-2.0 Mb','2.0-4.0 Mb','4.0-8.0 Mb','8.0-16.0 Mb','>16.0 Mb'])]
    
    # Save output
    all_count_ROH.to_csv(f"{output_prefix}_count_transpose.txt", sep='\t', index=False)
    all_sum_ROH.to_csv(f"{output_prefix}_sum_transpose.txt", sep='\t', index=False)

if __name__ == "__main__":
    args = parser.parse_args()
    
    with open(args.input, "r") as file:
        lines = file.readlines()
    
    transform(lines, args.out)
