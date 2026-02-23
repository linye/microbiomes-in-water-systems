"""
This script summarizes environment co-occurrence patterns of clusters based on
presence and reports both exclusive and non-exclusive combination counts.
"""

import itertools

# Input / output
INPUT_FILE = "cluster_env_sample_number_matrix.txt"
OUTPUT_FILE = "env_combination_results.txt"

# Environment order must match the column order in the matrix
GROUPS = ['AS', 'DW', 'GW', 'MW', 'NW', 'WW']

# Threshold for strict presence
STRICT_THRESHOLD = 2


def main():
    strict_dict = {}

    print("Parsing matrix file...")
    with open(INPUT_FILE, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(GROUPS) + 1:
                continue

            # Keep the same column selection logic as in the original script
            counts = list(map(int, parts[1:4] + parts[5:8]))

            present_groups = set()
            for idx, count in enumerate(counts):
                if count >= STRICT_THRESHOLD:
                    present_groups.add(GROUPS[idx])

            key = frozenset(present_groups)
            strict_dict[key] = strict_dict.get(key, 0) + 1

    print("Calculating all possible combinations...")
    all_combinations = []
    n = len(GROUPS)
    for r in range(1, n + 1):
        for combo in itertools.combinations(GROUPS, r):
            all_combinations.append(frozenset(combo))

    # Inclusive (normal) counts
    normal_count = {}
    for combo in all_combinations:
        total = 0
        for key in strict_dict:
            if combo.issubset(key):
                total += strict_dict[key]
        normal_count[combo] = total

    # Prepare output
    results = []
    for combo in all_combinations:
        exclusive_val = strict_dict.get(combo, 0)
        non-exclusive_val = normal_count.get(combo, 0)
        group_str = ','.join(sorted(combo))
        results.append((group_str, exclusive_val, non-exclusive_val))

    results.sort(key=lambda x: x[1], reverse=True)

    with open(OUTPUT_FILE, 'w') as out_file:
        out_file.write("group\texclusive\tnon-exclusive\n")
        for res in results:
            out_file.write(f"{res[0]}\t{res[1]}\t{res[2]}\n")

    print(f"Results written to {OUTPUT_FILE}")


if __name__ == '__main__':
    main()
