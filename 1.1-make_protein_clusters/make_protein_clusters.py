"""
This script processes a large protein clustering table and writes cluster memberships
by grouping sequences in a memory-efficient chunk-wise manner.
"""

import os
from collections import defaultdict

# Configuration
INPUT_FILE = "protein_cluster.tsv"
OUTPUT_FILE = "protein_clusters.txt"
CHUNK_SIZE = 1000000  # Process every 1 million lines

def process_clusters():
    """Process the clustering file in chunks."""
    # Initialize cluster dictionary
    cluster_map = defaultdict(set)

    # Open input file
    with open(INPUT_FILE, 'r') as f:
        chunk_count = 0
        total_count = 0

        for line in f:
            items = line.strip().split('\t')
            if len(items) < 2:
                continue

            a, b = items[0], items[1]
            cluster_map[a].add(b)

            # Counters
            total_count += 1
            chunk_count += 1

            # Periodically flush to disk to release memory
            if chunk_count >= CHUNK_SIZE:
                write_clusters(cluster_map, OUTPUT_FILE)
                cluster_map = defaultdict(set)  # Reset dictionary
                chunk_count = 0
                print(f'Processed {total_count} lines')

        # Handle the remaining data
        if cluster_map:
            write_clusters(cluster_map, OUTPUT_FILE)

    print(f'Finished processing! Total lines: {total_count}')

def write_clusters(cluster_map, output_file):
    """Write clustering results to file."""
    mode = 'a' if os.path.exists(output_file) else 'w'
    with open(output_file, mode) as out:
        for protein_id, members in cluster_map.items():
            out.write(f"{protein_id}\t{','.join(members)}\n")

if __name__ == "__main__":
    process_clusters()
    print(f'Clustering results saved to {OUTPUT_FILE}')
