#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter protein clusters that meet threshold criteria for each environment
and generate multiple output files. 
This script generates six files corresponding to 10%-60% thresholds in a single run.
"""

import os

# Input file path
input_file = "cluster_env_sample_number_matrix.txt"

# Output directory (can be customized)
output_dir = "./"
os.makedirs(output_dir, exist_ok=True)

# List of criteria for each file (from 10% to 60%)
criteria_dict_list = [
    {'file_suffix': '10%', 'criteria': {'DW': 44,  'GW': 23,  'NW': 79,  'WW': 59}},
    {'file_suffix': '20%', 'criteria': {'DW': 88,  'GW': 46,  'NW': 158, 'WW': 118}},
    {'file_suffix': '30%', 'criteria': {'DW': 132, 'GW': 69,  'NW': 237, 'WW': 177}},
    {'file_suffix': '40%', 'criteria': {'DW': 176, 'GW': 92,  'NW': 316, 'WW': 236}},
    {'file_suffix': '50%', 'criteria': {'DW': 220, 'GW': 115, 'NW': 395, 'WW': 295}},
    {'file_suffix': '60%', 'criteria': {'DW': 264, 'GW': 135, 'NW': 468, 'WW': 342}},
]

# Read header from input file
with open(input_file, 'r') as f:
    header = f.readline().strip().split('\t')

# Get column indices for each group (assumes header contains group names)
group_indices = {group: header.index(group) for group in ['DW', 'GW', 'NW', 'WW']}

# Read all lines into memory (skip header)
lines = [line for line in open(input_file, 'r')][1:]

# Generate output files for each threshold
for item in criteria_dict_list:
    suffix = item['file_suffix']
    criteria = item['criteria']
    output_file = os.path.join(output_dir, f"IPG_4group_{suffix}.txt")

    with open(output_file, 'w') as f_out:
        f_out.write('\t'.join(header) + '\n')  # Write header

        for line in lines:
            columns = line.strip().split('\t')
            meets_criteria = all(
                int(columns[group_indices[group]]) >= threshold
                for group, threshold in criteria.items()
            )
            if meets_criteria:
                f_out.write(line)

    print(f"✅ File generated: {output_file}, number of rows meeting criteria: {sum(1 for line in open(output_file)) - 1}")
