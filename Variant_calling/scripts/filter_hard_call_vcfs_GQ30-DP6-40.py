#!/usr/bin/env python

# Import modules
import gzip
import sys

# This script filters individual gzipped VCF files based on Genotype Quality (GQ) and Read Depth (DP).
# Created by Kevin Daly
# For each variant, it checks:
#   - If GQ >= 30
#   - If DP is between 6 and 40 (inclusive)
# Variants that do not meet these criteria will have their genotypes replaced with "./." (missing data).
# Usage: python script.py <input_vcf.gz> | bgzip -c > output_filtered.vcf.gz
#
with gzip.open(sys.argv[1]) as FILE:
        for LINE in FILE:
                LINE=LINE.rstrip("\n")
                if LINE.startswith("#"):
                        print LINE
                else:
                        SPLINE = LINE.split()
                        FIELDS = SPLINE[8]
                        GENO = SPLINE[9]
                        if GENO == ".":
                                print LINE
                        else:
                                if "GQ" in FIELDS:
                                        GQ = int(GENO.split(":")[8])
                                        DP = int(GENO.split(":")[2])
                                else:
                                        GQ = 30
                                        DP = int(GENO.split(":")[1])
                                if GQ >= 30 and (DP >= 6) and (DP <= 40):
                                        print LINE
                                else:
                                        print LINE.replace("0/0","./.").replace("0/1","./.").replace("1/0","./.").replace("1/1","./.")