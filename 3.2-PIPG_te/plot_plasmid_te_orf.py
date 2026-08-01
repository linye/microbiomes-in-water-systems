#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import os

# Set backend for non-interactive environment
plt.switch_backend('Agg')
sns.set_style("whitegrid")


def load_data():
    """Load TE data files"""
    # Read the two TE statistics files
    shared_df = pd.read_csv('plasmid_30_new_te.txt', sep='\t')
    random_df = pd.read_csv('plasmid_random_1000_te.txt',
                            sep='\t')

    # Add type labels
    shared_df['type'] = 'shared-plasmid'
    random_df['type'] = 'random-plasmid'

    # Classify random plasmids by size
    random_df['size_category'] = random_df['Length'].apply(
        lambda x: 'small-plasmid' if x < 100000 else 'large-plasmid'
    )

    # Calculate TE/ORF ratio (avoid division by zero)
    shared_df['TE_ORF_ratio'] = shared_df['TE_Count'] / shared_df['ORF_Count'].replace(0, np.nan)
    random_df['TE_ORF_ratio'] = random_df['TE_Count'] / random_df['ORF_Count'].replace(0, np.nan)

    return shared_df, random_df


def plot_te_orf_ratio_boxplot(shared_df, random_df):
    """Draw boxplot of TE/ORF ratio (log scale)"""
    # Prepare data
    small_plasmids = random_df[random_df['size_category'] == 'small-plasmid']
    large_plasmids = random_df[random_df['size_category'] == 'large-plasmid']

    plot_data = []

    for _, row in small_plasmids.iterrows():
        if not pd.isna(row['TE_ORF_ratio']) and row['TE_ORF_ratio'] > 0:
            plot_data.append({'type': 'small-plasmid', 'ratio': row['TE_ORF_ratio']})

    for _, row in large_plasmids.iterrows():
        if not pd.isna(row['TE_ORF_ratio']) and row['TE_ORF_ratio'] > 0:
            plot_data.append({'type': 'large-plasmid', 'ratio': row['TE_ORF_ratio']})

    for _, row in shared_df.iterrows():
        if not pd.isna(row['TE_ORF_ratio']) and row['TE_ORF_ratio'] > 0:
            plot_data.append({'type': 'shared-plasmid', 'ratio': row['TE_ORF_ratio']})

    plot_df = pd.DataFrame(plot_data)

    if len(plot_df) == 0:
        print("Warning: No valid TE/ORF ratio data to plot")
        return

    plt.figure(figsize=(10, 8))
    ax = sns.boxplot(data=plot_df, x='type', y='ratio',
                     order=['small-plasmid', 'large-plasmid', 'shared-plasmid'])
    ax.set_yscale('log')

    plt.xlabel('Plasmid Type')
    plt.ylabel('TE/ORF Ratio (log scale)')
    plt.title('TE/ORF Ratio Comparison by Plasmid Type (Log Scale)')
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig('te_orf_ratio_boxplot_log.png', dpi=300, bbox_inches='tight')
    plt.savefig('te_orf_ratio_boxplot_log.pdf', bbox_inches='tight')
    plt.close()

    # Print statistics
    print("\nTE/ORF ratio statistics:")
    for plasmid_type in ['small-plasmid', 'large-plasmid', 'shared-plasmid']:
        data = plot_df[plot_df['type'] == plasmid_type]['ratio']
        if len(data) > 0:
            print(f"{plasmid_type}: median={data.median():.4f}, mean={data.mean():.4f}, n={len(data)}")


def additional_analysis(shared_df, random_df):
    """Perform additional data analysis (printing only, no plots)"""

    # 1. TE prevalence analysis
    print("\n=== TE Prevalence Analysis ===")
    for name, df, desc in [('shared', shared_df, 'Shared plasmids'),
                           ('random', random_df, 'Random plasmids'),
                           ('small', random_df[random_df['size_category'] == 'small-plasmid'], 'Small plasmids'),
                           ('large', random_df[random_df['size_category'] == 'large-plasmid'], 'Large plasmids')]:
        te_positive = len(df[df['TE_Count'] > 0])
        total = len(df)
        if total > 0:
            print(f"{desc}: {te_positive}/{total} ({te_positive / total * 100:.1f}%) carry TEs")

    # 2. Correlation between length and TE count
    print("\n=== Correlation between Length and TE Count ===")
    for name, df in [('Shared', shared_df), ('Random', random_df)]:
        if len(df) > 1:
            te_positive_df = df[df['TE_Count'] > 0]
            if len(te_positive_df) > 1:
                corr = te_positive_df['Length'].corr(te_positive_df['TE_Count'])
                print(f"{name} plasmids - correlation (TE-positive only): {corr:.3f}")
            corr_all = df['Length'].corr(df['TE_Count'])
            print(f"{name} plasmids - correlation (all): {corr_all:.3f}")


def statistical_tests(shared_df, random_df):
    """Perform statistical tests"""
    print("\n=== Statistical Test Results ===")

    small_plasmids = random_df[random_df['size_category'] == 'small-plasmid']
    large_plasmids = random_df[random_df['size_category'] == 'large-plasmid']

    # 1. Length difference test
    if len(shared_df) > 1 and len(random_df) > 1:
        u_stat, p_value = stats.mannwhitneyu(shared_df['Length'], random_df['Length'])
        print(f"Shared vs Random plasmids length difference: p-value = {p_value:.4f}")

    # 2. TE/ORF ratio difference test
    groups = {
        'small-plasmid': small_plasmids['TE_ORF_ratio'].dropna(),
        'large-plasmid': large_plasmids['TE_ORF_ratio'].dropna(),
        'shared-plasmid': shared_df['TE_ORF_ratio'].dropna()
    }

    valid_groups = {}
    for k, v in groups.items():
        positive_values = v[v > 0]
        if len(positive_values) > 0:
            valid_groups[k] = positive_values

    if len(valid_groups) >= 2:
        h_stat, p_value = stats.kruskal(*valid_groups.values())
        print(f"TE/ORF ratio difference among plasmid types: p-value = {p_value:.4f}")

        if p_value < 0.05:
            print("\nPost-hoc pairwise comparisons:")
            from itertools import combinations
            for group1, group2 in combinations(valid_groups.keys(), 2):
                u_stat, p_val = stats.mannwhitneyu(valid_groups[group1], valid_groups[group2])
                print(f"{group1} vs {group2}: p-value = {p_val:.4f}")


def generate_summary_report(shared_df, random_df):
    """Generate summary report"""
    print("\n" + "=" * 60)
    print("                    TE Data Analysis Summary Report")
    print("=" * 60)

    print(f"\n1. Dataset Overview:")
    print(f"   - Shared plasmids: {len(shared_df)} plasmids")
    print(f"   - Random plasmids: {len(random_df)} plasmids")
    print(f"   - Small plasmids (<100kb): {len(random_df[random_df['size_category'] == 'small-plasmid'])}")
    print(f"   - Large plasmids (≥100kb): {len(random_df[random_df['size_category'] == 'large-plasmid'])}")

    print(f"\n2. Basic Statistics:")
    print(f"   Shared plasmids mean length: {shared_df['Length'].mean():.0f} ± {shared_df['Length'].std():.0f} bp")
    print(f"   Random plasmids mean length: {random_df['Length'].mean():.0f} ± {random_df['Length'].std():.0f} bp")

    print(f"\n3. TE carriage:")
    shared_te_positive = len(shared_df[shared_df['TE_Count'] > 0])
    random_te_positive = len(random_df[random_df['TE_Count'] > 0])
    print(
        f"   Shared plasmids TE-positive: {shared_te_positive}/{len(shared_df)} ({shared_te_positive / len(shared_df) * 100:.1f}%)")
    print(
        f"   Random plasmids TE-positive: {random_te_positive}/{len(random_df)} ({random_te_positive / len(random_df) * 100:.1f}%)")

    print(f"\n4. TE/ORF ratio summary:")
    small_plasmids = random_df[random_df['size_category'] == 'small-plasmid']
    large_plasmids = random_df[random_df['size_category'] == 'large-plasmid']

    for name, df in [('Small plasmids', small_plasmids),
                     ('Large plasmids', large_plasmids),
                     ('Shared plasmids', shared_df)]:
        ratios = df['TE_ORF_ratio'].dropna()
        positive_ratios = ratios[ratios > 0]
        if len(positive_ratios) > 0:
            print(
                f"   {name}: median={np.median(positive_ratios):.4f}, range={np.min(positive_ratios):.4f}-{np.max(positive_ratios):.4f}")


def main():
    """Main function"""
    # Check input files exist
    te_target = 'plasmid_30_new_te.txt'
    te_random = 'plasmid_random_1000_te.txt'
    if not os.path.exists(te_target) or not os.path.exists(te_random):
        print("Error: Please ensure the TE statistics files exist.")
        print("Required files:")
        print(f"1. {te_target}")
        print(f"2. {te_random}")
        return

    print("Starting TE data analysis...")

    # Load data
    shared_df, random_df = load_data()

    # Generate boxplot (log scale)
    plot_te_orf_ratio_boxplot(shared_df, random_df)

    # Additional analysis (prints only)
    additional_analysis(shared_df, random_df)

    # Statistical tests
    statistical_tests(shared_df, random_df)

    # Summary report
    generate_summary_report(shared_df, random_df)

    print("\nTE analysis complete! Generated plots:")
    print("- te_orf_ratio_boxplot_log.png/pdf: Boxplot of TE/ORF ratio (log scale)")


if __name__ == "__main__":
    main()