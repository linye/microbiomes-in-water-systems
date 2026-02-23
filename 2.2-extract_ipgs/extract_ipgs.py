#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract nucleotide sequences corresponding to IPG clusters for 4-group thresholds (10%-60%).
Generates six FASTA files: IPG_4group_10%.fa ... IPG_4group_60%.fa
"""

from Bio import SeqIO
import os

# Input nucleotide FASTA file
extract_fa = "nucl_rep_extract.fa"

# Input cluster files (10%-60%)
cluster_files = [
    "IPG_4group_10%.txt",
    "IPG_4group_20%.txt",
    "IPG_4group_30%.txt",
    "IPG_4group_40%.txt",
    "IPG_4group_50%.txt",
    "IPG_4group_60%.txt"
]

# Output FASTA files
output_files = [
    "IPG_4group_10%.fa",
    "IPG_4group_20%.fa",
    "IPG_4group_30%.fa",
    "IPG_4group_40%.fa",
    "IPG_4group_50%.fa",
    "IPG_4group_60%.fa"
]

# Make sure output directory exists
output_dir = "./"
os.makedirs(output_dir, exist_ok=True)

# Loop over all cluster files
for cluster_file, output_fa in zip(cluster_files, output_files):
    # Read Representative IDs from cluster file
    representatives = set()
    with open(cluster_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            rep_name = line.split()[0]
            # Remove the last "_" and anything after
            processed_name = "_".join(rep_name.split("_")[:-1])
            representatives.add(processed_name)

    # Extract matching sequences from nucleotide FASTA
    with open(output_fa, 'w') as output_handle:
        for record in SeqIO.parse(extract_fa, "fasta"):
            sequence_name = "_".join(record.id.split("_")[:-1])
            if sequence_name in representatives:
                SeqIO.write(record, output_handle, "fasta")

    print(f"✅ Extracted sequences saved to {output_fa}")
