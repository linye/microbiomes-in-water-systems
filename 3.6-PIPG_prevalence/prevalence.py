import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict

def create_result_folder():
    """Create results folder if it doesn't exist"""
    if not os.path.exists('result'):
        os.makedirs('result')
        print("Created result folder")
    else:
        print("Result folder already exists")

def load_plasmid_lengths():
    """
    Load plasmid length information and classify by size (≥100kb = Large, <100kb = Small)
    """
    print("Reading plasmid length file...")
    length_df = pd.read_csv('length.txt', sep='\t')

    # Classify plasmids by size
    length_df['Size_Category'] = length_df['Length'].apply(
        lambda x: 'Large-plasmid' if x >= 100000 else 'Small-plasmid'
    )

    print(f"Loaded {len(length_df)} plasmids")
    print(f"Size distribution: {length_df['Size_Category'].value_counts().to_dict()}")

    return length_df

def process_data(length_df):
    """
    Process coverage data: filter plasmids with coverage > 60% in each sample.
    Returns a DataFrame with columns: Sample, Group, Plasmid, Size_Category
    """
    # Read coverage file (format: first row is sample IDs, second row is 'length', then data)
    print("Reading coverage file...")
    coverage_df = pd.read_csv('coverage.txt', sep='\t', index_col=0)

    # Skip the 'length' row (second row)
    coverage_data = coverage_df.iloc[1:]
    coverage_data = coverage_data.astype(float)

    print("Filtering plasmids with coverage > 60%...")
    filtered_records = []

    for sample in coverage_data.index:
        sample_coverage = coverage_data.loc[sample]
        # Find plasmids with coverage > 0.6
        high_cov_plasmids = sample_coverage[sample_coverage > 0.6].index

        if len(high_cov_plasmids) > 0:
            # Extract group from sample name (first part before '_')
            group = sample.split('_')[0]
            for plasmid in high_cov_plasmids:
                # Get size category from length_df
                size_cat_series = length_df[length_df['Plasmid'] == plasmid]['Size_Category']
                size_cat = size_cat_series.iloc[0] if not size_cat_series.empty else 'Unknown'
                filtered_records.append({
                    'Sample': sample,
                    'Group': group,
                    'Plasmid': plasmid,
                    'Size_Category': size_cat
                })

    result_df = pd.DataFrame(filtered_records)

    print(f"Filtered data statistics:")
    print(f"Total data points (plasmid-sample pairs): {len(result_df)}")
    print(f"Number of samples: {result_df['Sample'].nunique()}")
    print(f"Group distribution: {result_df['Group'].value_counts().to_dict()}")
    print(f"Size category distribution: {result_df['Size_Category'].value_counts().to_dict()}")

    return result_df

def calculate_plasmid_prevalence(data_df):
    """
    Calculate the prevalence (percentage of samples) of each plasmid within each group,
    stratified by plasmid size category.
    """
    print("Calculating plasmid prevalence with size categories...")

    prevalence_records = []
    groups = sorted(data_df['Group'].unique())
    plasmids = sorted(data_df['Plasmid'].unique())
    size_categories = sorted(data_df['Size_Category'].unique())

    for group in groups:
        group_samples = data_df[data_df['Group'] == group]['Sample'].unique()
        n_samples_in_group = len(group_samples)

        for plasmid in plasmids:
            # Get plasmid size category
            plasmid_size = data_df[data_df['Plasmid'] == plasmid]['Size_Category'].iloc[0]

            # Count samples in which this plasmid appears in this group
            plasmid_samples = data_df[(data_df['Group'] == group) & (data_df['Plasmid'] == plasmid)]['Sample'].unique()
            n_samples_with_plasmid = len(plasmid_samples)

            prevalence_pct = (n_samples_with_plasmid / n_samples_in_group) * 100 if n_samples_in_group > 0 else 0

            prevalence_records.append({
                'Group': group,
                'Plasmid': plasmid,
                'Size_Category': plasmid_size,
                'Sample_Count': n_samples_with_plasmid,
                'Prevalence_Pct': prevalence_pct,
                'Total_Samples': n_samples_in_group
            })

    prevalence_df = pd.DataFrame(prevalence_records)

    print(f"Prevalence data calculated for {len(plasmids)} plasmids across {len(groups)} groups")
    print(f"Size categories: {size_categories}")

    return prevalence_df

def plot_plasmid_prevalence_histograms(prevalence_df):
    """
    Create stacked histograms (one row per group) showing:
    - Left: distribution of sample counts per plasmid (by size)
    - Right: distribution of prevalence percentages (by size)
    """
    print("Creating plasmid prevalence histograms with size categories...")

    groups = sorted(prevalence_df['Group'].unique())
    n_groups = len(groups)
    size_categories = sorted(prevalence_df['Size_Category'].unique())

    # Colors for size categories
    size_colors = {
        'Small-plasmid': 'lightblue',
        'Large-plasmid': 'lightcoral',
        'Unknown': 'lightgray'
    }

    fig, axes = plt.subplots(n_groups, 2, figsize=(14, 3 * n_groups))
    if n_groups == 1:
        axes = np.array([axes])

    # Bin edges for histograms
    sample_count_bins = np.arange(0, 11)          # 0–10 integers
    prevalence_pct_bins = np.linspace(0, 100, 11) # 0–100% in steps of 10

    for i, group in enumerate(groups):
        group_data = prevalence_df[prevalence_df['Group'] == group]

        # ----- Left: sample count histogram -----
        sample_count_data = []
        sample_count_labels = []
        for size_cat in size_categories:
            values = group_data[group_data['Size_Category'] == size_cat]['Sample_Count']
            if not values.empty:
                sample_count_data.append(values)
                sample_count_labels.append(size_cat)

        if sample_count_data:
            axes[i, 0].hist(sample_count_data, bins=sample_count_bins,
                            color=[size_colors.get(cat, 'gray') for cat in sample_count_labels],
                            label=sample_count_labels, alpha=0.8, edgecolor='black', stacked=True)

        axes[i, 0].set_title(f'{group} - Sample Count by Size', fontsize=12)
        axes[i, 0].set_xlabel('Number of Samples with Plasmid')
        axes[i, 0].set_ylabel('Number of Plasmids')
        axes[i, 0].set_xlim(-0.5, 10.5)
        axes[i, 0].grid(True, alpha=0.3)
        axes[i, 0].legend()

        # Statistics text
        total_plasmids = len(group_data)
        small_count = len(group_data[group_data['Size_Category'] == 'Small-plasmid'])
        large_count = len(group_data[group_data['Size_Category'] == 'Large-plasmid'])
        avg_count = group_data['Sample_Count'].mean()
        stats_text = f'Total: {total_plasmids}\nSmall: {small_count}\nLarge: {large_count}\nAvg: {avg_count:.1f}'
        axes[i, 0].text(0.02, 0.98, stats_text, transform=axes[i, 0].transAxes,
                        verticalalignment='top', fontsize=8,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # ----- Right: prevalence percentage histogram -----
        prev_data = []
        prev_labels = []
        for size_cat in size_categories:
            values = group_data[group_data['Size_Category'] == size_cat]['Prevalence_Pct']
            if not values.empty:
                prev_data.append(values)
                prev_labels.append(size_cat)

        if prev_data:
            axes[i, 1].hist(prev_data, bins=prevalence_pct_bins,
                            color=[size_colors.get(cat, 'gray') for cat in prev_labels],
                            label=prev_labels, alpha=0.8, edgecolor='black', stacked=True)

        axes[i, 1].set_title(f'{group} - Prevalence (%) by Size', fontsize=12)
        axes[i, 1].set_xlabel('Prevalence (%)')
        axes[i, 1].set_ylabel('Number of Plasmids')
        axes[i, 1].set_xlim(0, 100)
        axes[i, 1].grid(True, alpha=0.3)
        axes[i, 1].legend()

        avg_prev = group_data['Prevalence_Pct'].mean()
        common = len(group_data[group_data['Prevalence_Pct'] > 50])
        rare = len(group_data[group_data['Prevalence_Pct'] < 10])
        stats_text = f'Avg: {avg_prev:.1f}%\nCommon: {common}\nRare: {rare}'
        axes[i, 1].text(0.02, 0.98, stats_text, transform=axes[i, 1].transAxes,
                        verticalalignment='top', fontsize=8,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    plt.tight_layout()
    plt.savefig('result/plasmid_prevalence_histograms.png', dpi=300, bbox_inches='tight')
    plt.savefig('result/plasmid_prevalence_histograms.pdf', bbox_inches='tight')
    plt.show()
    return fig

def plot_prevalence_summary(prevalence_df):
    """
    Create summary plots:
    - Average prevalence per group and size category
    - Number of common (>50%) and rare (<10%) plasmids per group and size
    - Overall prevalence distribution by size
    """
    print("Creating plasmid prevalence summary plots with size categories...")

    size_colors = {
        'Small-plasmid': 'lightblue',
        'Large-plasmid': 'lightcoral',
        'Unknown': 'lightgray'
    }

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # 1. Average prevalence per group and size
    group_size_avg = prevalence_df.groupby(['Group', 'Size_Category'])['Prevalence_Pct'].mean().unstack()
    group_size_avg.plot(kind='bar', ax=axes[0, 0],
                        color=[size_colors.get(c, 'gray') for c in group_size_avg.columns])
    axes[0, 0].set_title('Average Plasmid Prevalence by Group and Size', fontsize=14)
    axes[0, 0].set_xlabel('Group')
    axes[0, 0].set_ylabel('Average Prevalence (%)')
    axes[0, 0].tick_params(axis='x', rotation=45)
    axes[0, 0].legend(title='Size Category')

    # 2. Number of common plasmids (>50%) per group and size
    common = prevalence_df[prevalence_df['Prevalence_Pct'] > 50]
    common_counts = common.groupby(['Group', 'Size_Category']).size().unstack(fill_value=0)
    common_counts.plot(kind='bar', ax=axes[0, 1],
                       color=[size_colors.get(c, 'gray') for c in common_counts.columns])
    axes[0, 1].set_title('Number of Common Plasmids (>50% prevalence) by Group and Size', fontsize=14)
    axes[0, 1].set_xlabel('Group')
    axes[0, 1].set_ylabel('Number of Plasmids')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].legend(title='Size Category')

    # 3. Number of rare plasmids (<10%) per group and size
    rare = prevalence_df[prevalence_df['Prevalence_Pct'] < 10]
    rare_counts = rare.groupby(['Group', 'Size_Category']).size().unstack(fill_value=0)
    rare_counts.plot(kind='bar', ax=axes[1, 0],
                     color=[size_colors.get(c, 'gray') for c in rare_counts.columns])
    axes[1, 0].set_title('Number of Rare Plasmids (<10% prevalence) by Group and Size', fontsize=14)
    axes[1, 0].set_xlabel('Group')
    axes[1, 0].set_ylabel('Number of Plasmids')
    axes[1, 0].tick_params(axis='x', rotation=45)
    axes[1, 0].legend(title='Size Category')

    # 4. Overall prevalence distribution by size
    for size_cat in sorted(prevalence_df['Size_Category'].unique()):
        data = prevalence_df[prevalence_df['Size_Category'] == size_cat]['Prevalence_Pct']
        axes[1, 1].hist(data, bins=30, alpha=0.7,
                        color=size_colors.get(size_cat, 'gray'),
                        label=size_cat, edgecolor='black')

    axes[1, 1].set_title('Distribution of Plasmid Prevalence by Size', fontsize=14)
    axes[1, 1].set_xlabel('Prevalence (%)')
    axes[1, 1].set_ylabel('Number of Plasmid-Group Pairs')
    axes[1, 1].axvline(x=50, color='red', linestyle='--', label='Common (50%)')
    axes[1, 1].axvline(x=10, color='blue', linestyle='--', label='Rare (10%)')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig('result/plasmid_prevalence_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig('result/plasmid_prevalence_summary.pdf', bbox_inches='tight')
    plt.show()

def save_data(data_df, prevalence_df, length_df):
    """
    Save filtered data and prevalence results.
    """
    # Save the filtered plasmid-sample pairs (with coverage > 60%)
    data_df.to_csv('result/filtered_plasmid_samples.csv', index=False)

    # Save prevalence data
    prevalence_df.to_csv('result/plasmid_prevalence_data.csv', index=False)

    # Save plasmid length information with size categories
    length_df.to_csv('result/plasmid_length_info.csv', index=False)

    print("\nData saved:")
    print("- result/filtered_plasmid_samples.csv: Filtered plasmid-sample pairs (coverage > 60%)")
    print("- result/plasmid_prevalence_data.csv: Prevalence statistics per group and plasmid")
    print("- result/plasmid_length_info.csv: Plasmid lengths and size categories")

def main():
    """Main function"""
    print("Starting plasmid prevalence analysis...")

    # Check required input files
    required_files = ['coverage.txt', 'length.txt']
    for f in required_files:
        if not os.path.exists(f):
            print(f"Error: {f} not found.")
            return

    try:
        create_result_folder()

        # Load plasmid lengths and size categories
        length_df = load_plasmid_lengths()

        # Process coverage data (filter > 60%)
        filtered_data = process_data(length_df)

        if len(filtered_data) == 0:
            print("No plasmids with coverage > 60% found in any sample.")
            return

        # Calculate prevalence
        prevalence_df = calculate_plasmid_prevalence(filtered_data)

        # Generate prevalence histograms
        plot_plasmid_prevalence_histograms(prevalence_df)

        # Generate summary plots
        plot_prevalence_summary(prevalence_df)

        # Save results
        save_data(filtered_data, prevalence_df, length_df)

        print("\nProcessing completed!")
        print("Generated plots in 'result' folder:")
        print("- plasmid_prevalence_histograms.png/pdf: Stacked histograms by group and size")
        print("- plasmid_prevalence_summary.png/pdf: Summary bar plots and distribution")

    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()