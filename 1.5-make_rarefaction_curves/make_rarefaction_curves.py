"""
This script computes rarefaction curves of protein clusters by progressively
adding samples. It generates a global curve and environment-specific curves.
"""

from collections import defaultdict

# Input files
CLUSTER_FILE = "protein_cluster.tsv"
ENV_FILE = "env.txt"

# Output prefix
GLOBAL_OUT = "rarecurve.txt"


def load_environment_info():
    """Load sample-to-environment mapping."""
    sample_to_env = {}
    env_samples = defaultdict(set)

    with open(ENV_FILE) as f:
        next(f)  # skip header if present
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            sample, env = parts[0], parts[1]
            sample_to_env[sample] = env
            env_samples[env].add(sample)

    print(f"Loaded {len(sample_to_env)} samples across {len(env_samples)} environments")
    return sample_to_env, env_samples


def load_cluster_data():
    """Load pairwise cluster relationships."""
    data = []
    with open(CLUSTER_FILE) as f:
        for line in f:
            data.append(line.strip())
    print(f"Loaded {len(data)} cluster links")
    return data


def compute_curve(data, samples_subset, step, outfile):
    """Compute rarefaction curve for a subset of samples."""
    samples_subset = list(samples_subset)
    max_n = len(samples_subset)

    with open(outfile, 'w') as out:
        for x in range(step, max_n + 1, step):
            selected = set(samples_subset[:x])

            observed = set()
            for line in data:
                items = line.split('\t')
                if len(items) < 2:
                    continue

                name1 = items[0].split('_')[0]
                name2 = items[1].split('_')[0]

                if name1 in selected and name2 in selected:
                    observed.add(items[0])

            out.write(f"{x}\t{len(observed)}\n")

    print(f"Curve written to {outfile}")


def main():
    sample_to_env, env_samples = load_environment_info()
    data = load_cluster_data()

    # ---- Global rarefaction ----
    all_samples = set(sample_to_env.keys())
    compute_curve(data, all_samples, step=50, outfile=GLOBAL_OUT)

    # ---- Environment-specific rarefaction ----
    for env, samples in env_samples.items():
        outfile = f"rarecurve-{env}.txt"
        compute_curve(data, samples, step=20, outfile=outfile)


if __name__ == "__main__":
    main()
