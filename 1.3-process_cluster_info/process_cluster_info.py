"""
This script parses cluster memberships and metadata to generate
group-level sample counts and study counts for each representative cluster.
"""

import os
from collections import defaultdict
import time

# Input files
CLUSTER_FILE = "protein_clusters.txt"
GROUP_FILE = "env.txt"
STUDY_FILE = "study.txt"

# Output files
OUTPUT_SAMPLE_FILE = "cluster_env_sample_number_matrix.txt"
OUTPUT_STUDY_FILE = "cluster_env_study_number_matrix.txt"


def process_group_file():
    """Parse the group file and return sample-to-group mapping and sorted groups."""
    sample_to_group = {}
    groups = set()

    print("Parsing group file...")
    with open(GROUP_FILE, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            sample, group = parts[0], parts[1]
            sample_to_group[sample] = group
            groups.add(group)

            if line_num % 1000000 == 0:
                print(f"Processed {line_num} lines of group data")

    sorted_groups = sorted(groups)
    print(f"Group parsing finished. Total groups: {len(sorted_groups)}")
    return sample_to_group, sorted_groups


def process_study_file():
    """Parse the study file and return sample-to-study mapping."""
    sample_to_study = {}

    print("Parsing study file...")
    with open(STUDY_FILE, 'r') as f:
        # Skip header
        next(f)

        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            sample = parts[0]
            study_id = parts[2]
            sample_to_study[sample] = study_id

            if line_num % 100000 == 0:
                print(f"Processed {line_num} lines of study data")

    print(f"Study parsing finished. Total samples: {len(sample_to_study)}")
    return sample_to_study


def process_clusters(sample_to_group, sample_to_study, groups):
    """Generate matrices of sample counts and study counts per group."""
    print("Processing cluster file and generating matrices...")
    start_time = time.time()
    processed_count = 0

    with open(OUTPUT_SAMPLE_FILE, 'w') as out_sample_f, \
         open(OUTPUT_STUDY_FILE, 'w') as out_study_f:

        # Write headers
        header = "Representative\t" + "\t".join(groups) + "\n"
        out_sample_f.write(header)
        out_study_f.write(header)

        with open(CLUSTER_FILE, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                representative = parts[0]

                # Include representative itself
                all_sequences = [representative]
                if len(parts) > 1:
                    all_sequences.extend(parts[1].split(','))

                # Extract sample IDs
                samples = [seq.split('_')[0] for seq in all_sequences]

                # Initialize counters
                group_sample_counts = defaultdict(int)
                group_study_counts = defaultdict(set)

                # Counting
                for sample in samples:
                    group = sample_to_group.get(sample)
                    if group:
                        group_sample_counts[group] += 1

                        study_id = sample_to_study.get(sample)
                        if study_id:
                            group_study_counts[group].add(study_id)

                # Prepare output rows
                sample_counts = [str(group_sample_counts.get(group, 0)) for group in groups]
                study_counts = [str(len(group_study_counts.get(group, set()))) for group in groups]

                # Write
                out_sample_f.write(f"{representative}\t" + "\t".join(sample_counts) + "\n")
                out_study_f.write(f"{representative}\t" + "\t".join(study_counts) + "\n")

                processed_count += 1
                if line_num % 100000 == 0:
                    elapsed = time.time() - start_time
                    rate = processed_count / elapsed if elapsed > 0 else 0
                    print(f"Processed {line_num} lines, {rate:.2f} clusters/s, elapsed {elapsed:.2f}s")

    print(f"Cluster processing finished. Total clusters: {processed_count}")
    return processed_count


def main():
    # Step 1: group mapping
    sample_to_group, groups = process_group_file()

    # Step 2: study mapping
    sample_to_study = process_study_file()

    # Step 3: cluster processing
    process_clusters(sample_to_group, sample_to_study, groups)

    print("Results written to:")
    print(f"1. Sample count matrix: {OUTPUT_SAMPLE_FILE}")
    print(f"2. Study count matrix: {OUTPUT_STUDY_FILE}")


if __name__ == "__main__":
    main()
