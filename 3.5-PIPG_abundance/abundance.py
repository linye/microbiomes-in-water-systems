import pandas as pd
import os
import glob
import re
from collections import defaultdict

# Configuration file paths
LENGTHS_FILE = "plasmid_lengths.txt"
COVERAGE_DIR = "coverage_counts"
OUTPUT_PREFIX = "plasmid_abundance"

# 1. Read plasmid length information
lengths_df = pd.read_csv(LENGTHS_FILE, sep='\t', header=None,
                         names=['plasmid_id', 'length'])
print(f"Loaded plasmid length information: {len(lengths_df)} plasmids")

# 2. Collect all sample files
sample_files = glob.glob(os.path.join(COVERAGE_DIR, "*_plasmid_coverage_counts.txt"))
if not sample_files:
    print(f"Error: No coverage count files found in directory {COVERAGE_DIR}")
    exit(1)

print(f"Found {len(sample_files)} sample files")

# 3. Create base dataframe
base_df = lengths_df.copy()

# 4. Initialize result dataframes
counts_df = lengths_df.copy()
rpkm_df = lengths_df.copy()
tpm_df = lengths_df.copy()
coverage_pct_df = lengths_df.copy()

# 5. Process each sample file
sample_stats = defaultdict(dict)  # Store statistics for each sample

for file_path in sample_files:
    # Extract complete sample name (e.g., NW_SRR11088437)
    filename = os.path.basename(file_path)

    # Use regex to extract the full sample ID
    # Pattern: everything before the first "_plasmid"
    match = re.match(r'(.+?)_plasmid', filename)
    if match:
        sample_name = match.group(1)
    else:
        # If regex fails, use filename (without extension) as fallback
        sample_name = os.path.splitext(filename)[0]
        print(f"Warning: Could not extract sample ID from filename {filename}, using {sample_name} as fallback")

    print(f"Processing sample: {sample_name}")

    # Read sample data
    sample_df = pd.read_csv(file_path, sep='\t', header=None,
                            names=['plasmid_id', 'count'])

    # Merge with base dataframe
    merged_df = pd.merge(base_df, sample_df, on='plasmid_id', how='left')
    merged_df['count'] = merged_df['count'].fillna(0)  # Handle missing values

    # Calculate total reads
    total_reads = merged_df['count'].sum()
    sample_stats[sample_name]['total_reads'] = total_reads

    # Calculate RPKM (avoid division by zero)
    if total_reads > 0:
        merged_df['RPKM'] = (merged_df['count'] * 1e9) / (total_reads * merged_df['length'])
    else:
        merged_df['RPKM'] = 0
        print(f"Warning: Sample {sample_name} total reads is 0, cannot compute RPKM")

    # Calculate RPK (for TPM)
    merged_df['RPK'] = merged_df['count'] / merged_df['length'].clip(lower=1) * 1000

    # Calculate total RPK (avoid division by zero)
    total_rpk = merged_df['RPK'].sum()
    if total_rpk > 0:
        merged_df['TPM'] = (merged_df['RPK'] / total_rpk) * 1e6
    else:
        merged_df['TPM'] = 0
        print(f"Warning: Sample {sample_name} total RPK is 0, cannot compute TPM")

    # Calculate coverage percentage
    if total_reads > 0:
        merged_df['coverage_pct'] = (merged_df['count'] / total_reads) * 100
    else:
        merged_df['coverage_pct'] = 0

    # Add to respective metric dataframes
    counts_df[sample_name] = merged_df['count']
    rpkm_df[sample_name] = merged_df['RPKM']
    tpm_df[sample_name] = merged_df['TPM']
    coverage_pct_df[sample_name] = merged_df['coverage_pct']

    # Record statistics
    sample_stats[sample_name]['rpkm_mean'] = merged_df['RPKM'].mean()
    sample_stats[sample_name]['tpm_mean'] = merged_df['TPM'].mean()
    sample_stats[sample_name]['max_rpkm'] = merged_df['RPKM'].max()
    sample_stats[sample_name]['max_tpm'] = merged_df['TPM'].max()

# 6. Save results to different files
# Raw counts
counts_df.to_csv(f"{OUTPUT_PREFIX}_counts.csv", index=False)
print(f"Raw counts saved to: {OUTPUT_PREFIX}_counts.csv")

# RPKM values
rpkm_df.to_csv(f"{OUTPUT_PREFIX}_RPKM.csv", index=False)
print(f"RPKM values saved to: {OUTPUT_PREFIX}_RPKM.csv")

# TPM values
tpm_df.to_csv(f"{OUTPUT_PREFIX}_TPM.csv", index=False)
print(f"TPM values saved to: {OUTPUT_PREFIX}_TPM.csv")

# Coverage percentage
coverage_pct_df.to_csv(f"{OUTPUT_PREFIX}_coverage_pct.csv", index=False)
print(f"Coverage percentages saved to: {OUTPUT_PREFIX}_coverage_pct.csv")

# 7. Save sample summary statistics
stats_df = pd.DataFrame.from_dict(sample_stats, orient='index')
stats_df['sample'] = stats_df.index
stats_df = stats_df[['sample', 'total_reads', 'rpkm_mean', 'tpm_mean', 'max_rpkm', 'max_tpm']]
stats_df.to_csv("sample_summary_stats.csv", index=False)
print("Sample summary statistics saved to sample_summary_stats.csv")

# 8. Print final summary
print("\nProcessing complete! Summary:")
print(f"Number of samples: {len(sample_files)}")
print(f"Number of plasmids: {len(lengths_df)}")
print(f"Total reads range: {stats_df['total_reads'].min():,} - {stats_df['total_reads'].max():,}")
print(f"Mean RPKM range: {stats_df['rpkm_mean'].min():.2f} - {stats_df['rpkm_mean'].max():.2f}")
print(f"Mean TPM range: {stats_df['tpm_mean'].min():.2f} - {stats_df['tpm_mean'].max():.2f}")

# 9. Show example of first 5 rows
print("\nFirst 5 rows of RPKM file:")
print(rpkm_df.head().to_string(index=False))