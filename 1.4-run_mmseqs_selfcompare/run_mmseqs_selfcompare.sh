#!/bin/bash

# This script performs self-comparison for each FASTA file using MMseqs2
# and exports pairwise alignments into TSV format.

for f in clusters30/*.fasta; do
    base=$(basename "$f" .fasta)

    echo "Processing $base..."

    # 1. Create MMseqs2 database
    mmseqs createdb "$f" "dbs/${base}_DB"

    # 2. Self-search (search-type 3 indicates nucleotide sequences)
    mmseqs search "dbs/${base}_DB" "dbs/${base}_DB" "results/${base}_res" tmp --search-type 3

    # 3. Convert alignment results to TSV format
    mmseqs convertalis "dbs/${base}_DB" "dbs/${base}_DB" "results/${base}_res" "results/${base}.tsv"
done

echo "All self-comparisons completed. Results are stored in the results/ directory."
