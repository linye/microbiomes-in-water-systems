#!/usr/bin/env python3
import os
from collections import defaultdict


def parse_fasta_lengths(fasta_file):
    """Parse FASTA file, extract sequence IDs and lengths"""
    seq_lengths = {}
    current_id = ""
    current_seq = ""

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    seq_lengths[current_id] = len(current_seq.replace('\n', '').replace(' ', ''))
                # Extract sequence ID (part before the first space)
                current_id = line[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line

        # Process the last sequence
        if current_id and current_seq:
            seq_lengths[current_id] = len(current_seq.replace('\n', '').replace(' ', ''))

    return seq_lengths


def count_orfs_per_plasmid(orf_fasta_file):
    """Count the number of ORFs per plasmid"""
    plasmid_orf_count = defaultdict(int)

    with open(orf_fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract ORF name, e.g., NZ_CP138514.1_1
                orf_name = line[1:].split()[0]
                # Extract plasmid ID (part before the last underscore)
                plasmid_id = '_'.join(orf_name.split('_')[:-1])
                plasmid_orf_count[plasmid_id] += 1

    return plasmid_orf_count


def count_args_per_plasmid(blast_result_file):
    """Count the number of ARGs (antimicrobial resistance genes) matched per plasmid"""
    plasmid_arg_count = defaultdict(int)
    processed_orfs = set()  # Used for deduplication, avoid counting multiple hits from the same ORF

    with open(blast_result_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            orf_name = parts[0]
            # Extract plasmid ID (part before the last underscore)
            plasmid_id = '_'.join(orf_name.split('_')[:-1])

            # Use ORF name for deduplication, ensure each ORF counted only once
            if orf_name not in processed_orfs:
                plasmid_arg_count[plasmid_id] += 1
                processed_orfs.add(orf_name)

    return plasmid_arg_count


def generate_summary_table(original_fasta, orf_fasta, blast_result, output_file):
    """Generate summary table"""
    print(f"Processing data: {original_fasta}")

    # Get plasmid lengths
    plasmid_lengths = parse_fasta_lengths(original_fasta)
    print(f"Found {len(plasmid_lengths)} plasmid sequences")

    # Get ORF counts per plasmid
    orf_counts = count_orfs_per_plasmid(orf_fasta)
    print(f"Counted ORFs for {len(orf_counts)} plasmids")

    # Get ARG counts per plasmid
    arg_counts = count_args_per_plasmid(blast_result)
    print(f"Counted ARGs for {len(arg_counts)} plasmids")

    # Generate summary table
    with open(output_file, 'w') as out:
        # Write header
        out.write("Plasmid_ID\tLength\tORF_Count\tARG_Count\n")

        # Write data
        for plasmid_id in plasmid_lengths:
            length = plasmid_lengths.get(plasmid_id, 0)
            orf_count = orf_counts.get(plasmid_id, 0)
            arg_count = arg_counts.get(plasmid_id, 0)

            out.write(f"{plasmid_id}\t{length}\t{orf_count}\t{arg_count}\n")

    print(f"Results saved to: {output_file}")
    print("-" * 50)


def main():
    # Target plasmid data file paths
    target_plasmid_fasta = "../plasmid_30_new.fasta"
    target_orf_fasta = "../plasmid_30_new_nucleotide_seq.fasta"
    target_blast_result = "plasmid_30_new_resfinder_cov70.txt"
    target_output = "plasmid_30_new_arg.txt"

    # Random plasmid data file paths
    random_plasmid_fasta = "../plasmid_random_1000.fasta"
    random_orf_fasta = "../plasmid_random_1000_nucleotide_seq.fasta"
    random_blast_result = "plasmid_random_1000_resfinder_cov70.txt"
    random_output = "plasmid_random_1000_arg.txt"

    # Check file existence
    files_to_check = [
        (target_plasmid_fasta, "Target plasmid FASTA"),
        (target_orf_fasta, "Target plasmid ORF FASTA"),
        (target_blast_result, "Target plasmid BLAST result"),
        (random_plasmid_fasta, "Random plasmid FASTA"),
        (random_orf_fasta, "Random plasmid ORF FASTA"),
        (random_blast_result, "Random plasmid BLAST result")
    ]

    for file_path, description in files_to_check:
        if not os.path.exists(file_path):
            print(f"Error: {description} does not exist - {file_path}")
            return

    # Process target plasmid data
    generate_summary_table(target_plasmid_fasta, target_orf_fasta, target_blast_result, target_output)

    # Process random plasmid data
    generate_summary_table(random_plasmid_fasta, random_orf_fasta, random_blast_result, random_output)

    print("All processing complete!")


if __name__ == "__main__":
    main()