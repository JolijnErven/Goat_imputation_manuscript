#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file, which is the roh_file output from extract_roh.py", type=str)
parser.add_argument("-o", "--output_prefix", help="Output file prefix", type=str)

def demographic_whole_genome(roh_file, prefix):
    """Calculate ROH demographics for the whole genome while preserving output format."""
    roh_categories = {
        "0.5-1.0 Mb": (0, 1000000),
        "1.0-2.0 Mb": (1000000, 2000000),
        "2.0-4.0 Mb": (2000000, 4000000),
        "4.0-8.0 Mb": (4000000, 8000000),
        "8.0-16.0 Mb": (8000000, 16000000),
        ">16.0 Mb": (16000000, float("inf"))
    }
    genome_size = 2466191000 # autosome size in KB - Constant for Froh calculation
    
    record, sample_data = [], {}
    current_sample = None
    
    for line in roh_file:
        line_parts = line.split()
        if line_parts[0] == "FID":
            continue
        
        sample, roh_size = line_parts[1], float(line_parts[5])
        
        if sample != current_sample:
            if current_sample is not None:
                summed_roh = sum(sample_data[cat][1] for cat in roh_categories)
                record.append([current_sample] + [val for cat in roh_categories for val in sample_data[cat]] + [summed_roh, summed_roh / genome_size])
            
            current_sample = sample
            sample_data = {cat: [0, 0] for cat in roh_categories}  # [count, sum]
        
        for cat, (low, high) in roh_categories.items():
            if low < roh_size <= high:
                sample_data[cat][0] += 1
                sample_data[cat][1] += roh_size
                break
    
    if current_sample:
        summed_roh = sum(sample_data[cat][1] for cat in roh_categories)
        record.append([current_sample] + [val for cat in roh_categories for val in sample_data[cat]] + [summed_roh, summed_roh / genome_size])
    
    with open(f"ROH_demo_plots_{prefix}.txt", 'w') as f:
        headers = [
            "ID",
            "roh count 0.5-1.0 Mb", "0.5-1.0 Mb", "froh 0.5-1.0 Mb",
            "roh count 1.0-2.0 Mb", "1.0-2.0 Mb", "froh 1.0-2.0 Mb",
            "roh count 2.0-4.0 Mb", "2.0-4.0 Mb", "froh 2.0-4.0 Mb",
            "roh count 4.0-8.0 Mb", "4.0-8.0 Mb", "froh 4.0-8.0 Mb",
            "roh count 8.0-16.0 Mb", "8.0-16.0 Mb", "froh 8.0-16.0 Mb",
            "roh count >16.0 Mb", ">16.0 Mb", "froh >16.0 Mb",
            "sum_ROH", "Froh"
        ]
        f.write("\t".join(headers) + "\n")
        for row in record:
            formatted_row = [str(row[i]) + ("\t" + str(row[i+1]) + "\t" + str(row[i+1]/genome_size)) for i in range(1, len(row)-2, 2)]
            f.write(str(row[0]) + "\t" + "\t".join(formatted_row) + "\t" + str(row[-2]) + "\t" + str(row[-1]) + "\n")
    
    return f

if __name__ == "__main__":
    args = parser.parse_args()
    with open(args.input) as f:
        roh_data = f.readlines()
    demographic_whole_genome(roh_data, args.output_prefix)
