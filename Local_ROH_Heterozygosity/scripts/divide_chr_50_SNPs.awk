BEGIN {
    window_size = 50;
    heterozygote_count = 0;
    snp_count = 0;
    start_position = 0;
    end_position = 0;
}

#Process lines that do not start with #
!/^#/ {
    snp_count++;  # Increment the SNP count
    position = $2;  # Extract the position (2nd column)

    # If it is the first SNP, set the start position
    if (snp_count == 1) {
        start_position = position;
    }

    # Update the end position with the current SNP position
    end_position = position;

    # Process each genotype in the VCF file (Only works for individual VCF)
    for (i = 10; i <= NF; i++) {
        # Extract genotype from the FORMAT field (e.g., GT=0/1 or GT=1|0)
        split($i, geno, ":");

        # Check phased (|) and unphased (/)
        if (index(geno[1], "/") > 0) {
            split(geno[1], gt, "/");  # Unphased
        } else if (index(geno[1], "|") > 0) {
            split(geno[1], gt, "|");  # Phased
        } else {
            continue;  # Skip if no valid genotype
        }

        # Check if the genotype is heterozygous (e.g., 0/1, 1/0, 0|1, 1|0)
        if (gt[1] != gt[2]) {
            heterozygote_count++;
        }
    }

    # If 50 SNPs have been processed, output the count and reset
    if (snp_count % window_size == 0) {
        print start_position, end_position, heterozygote_count;
        heterozygote_count = 0;  # Reset for the next window
        start_position = position + 1;  # Set new start position for next window
    }
}

