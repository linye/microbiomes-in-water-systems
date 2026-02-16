#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script calculates the average proportion of sequences in each group above varying thresholds and 
generates a TSV file with average proportions and a line plot for visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input and output file paths
input_file = "IPG_4group_10%.txt"
output_file = "cluster.txt"

# Total number of samples per group
group_totals = {
    "AS": 607,
    "DW": 411,
    "GW": 214,
    "MW": 778,
    "NW": 651,
    "WW": 526,
}

# Column names for reading the input file
columns = ["Sequence", "AS", "DW", "GW", "MW", "NW", "WW"]
data = pd.read_csv(input_file, sep="\t", names=columns, skiprows=1)

# Percentage thresholds to evaluate
percentages = range(10, 61)  # 10% to 60%

# Store results
results = []

# Main loop over thresholds
for pct in percentages:
    threshold = pct / 100.0

    # Select sequences meeting threshold criteria
    selected_sequences = data[
        # Uncomment these lines if AS or MW groups should be included
        #(data["AS"] >= group_totals["AS"] * threshold) &
        (data["DW"] >= group_totals["DW"] * threshold) &
        (data["GW"] >= group_totals["GW"] * threshold) &
        #(data["MW"] >= group_totals["MW"] * threshold) &
        (data["NW"] >= group_totals["NW"] * threshold) &
        (data["WW"] >= group_totals["WW"] * threshold)
    ]

    if not selected_sequences.empty:
        # Compute average proportion for each group
        avg_ratios = {
            group: selected_sequences[group].mean() / group_totals[group] for group in group_totals
        }
        results.append([pct] + [avg_ratios[group] for group in group_totals])
    else:
        results.append([pct] + [0] * len(group_totals))

# Save results to TSV
df_results = pd.DataFrame(results, columns=["Percentage"] + list(group_totals.keys()))
df_results.to_csv(output_file, sep="\t", index=False)

# Plot average proportion curves
plt.figure(figsize=(10, 6))
for group in group_totals:
    plt.plot(df_results["Percentage"].to_numpy(),
             df_results[group].to_numpy(),
             label=group,
             marker="o")

plt.xlabel("Percentage Threshold (%)")
plt.ylabel("Average Sample Proportion")
plt.title("Average Sample Proportion by Percentage Threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save plot in PNG and PDF formats
plt.savefig("cluster_plot.png")
plt.savefig("cluster_plot.pdf")
plt.show()

print(f"✅ Results saved to {output_file} and plot saved as cluster_plot.png/pdf")
