#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script computes Jaccard similarity, directional gene flow, UpSet plots,
and net directional flow networks among aquatic environments based on cluster occurrence data.
"""

import itertools
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from upsetplot import UpSet, from_memberships


# ===============================
# Utilities
# ===============================

def save_figure(basepath, dpi=300):
    """Save current figure as PNG and PDF"""
    plt.savefig(f"{basepath}.png", dpi=dpi)
    plt.savefig(f"{basepath}.pdf")


# ===============================
# Data loading and processing
# ===============================

def load_data(file):
    df = pd.read_csv(file, sep='\t')
    df['group_set'] = df['group'].apply(lambda x: frozenset(x.split(',')))
    envs = sorted(set(itertools.chain.from_iterable(df['group_set'])))
    return df, envs


def calc_env_total(df, envs):
    total = {e: 0 for e in envs}
    for _, r in df.iterrows():
        for e in r['group_set']:
            total[e] += r['exclusive']
    return total


def calc_shared_only_df(df):
    return df[df['group_set'].apply(lambda x: len(x) >= 2)].copy()


def calc_jaccard(df, envs, env_total):
    mat = pd.DataFrame(0.0, index=envs, columns=envs)
    for e1, e2 in itertools.combinations(envs, 2):
        shared = df.loc[
            df['group_set'].apply(lambda x: {e1, e2}.issubset(x)),
            'exclusive'
        ].sum()
        union = env_total[e1] + env_total[e2] - shared
        val = shared / union if union > 0 else 0
        mat.loc[e1, e2] = mat.loc[e2, e1] = val
    np.fill_diagonal(mat.values, 1)
    return mat


def calc_direction_gene_flow(df, envs, env_total):
    mat = pd.DataFrame(0.0, index=envs, columns=envs)
    for src, tgt in itertools.product(envs, envs):
        if src == tgt:
            continue
        shared = df.loc[
            df['group_set'].apply(lambda x: {src, tgt}.issubset(x)),
            'exclusive'
        ].sum()
        mat.loc[src, tgt] = shared / env_total[tgt] if env_total[tgt] > 0 else 0
    return mat


def calc_net_flow(direction_mat):
    net_flow = {}
    for env in direction_mat.index:
        outflow = direction_mat.loc[env, :].sum()
        inflow = direction_mat.loc[:, env].sum()
        net_flow[env] = outflow - inflow
    return net_flow


# ===============================
# Plotting functions
# ===============================

def draw_heatmap(mat, title, outfile):
    plt.figure(figsize=(8, 6))
    sns.heatmap(mat, cmap='viridis', square=True,
                cbar_kws={'label': 'Similarity'})
    plt.title(title)
    plt.tight_layout()
    save_figure(outfile)
    plt.close()


def draw_direction_gene_flow_network(mat, threshold, outfile):
    G = nx.DiGraph()
    for src, tgt in itertools.product(mat.index, mat.columns):
        val = mat.loc[src, tgt]
        if val >= threshold:
            G.add_edge(src, tgt, weight=val)

    pos = nx.spring_layout(G, seed=42, weight='weight')
    weights = [G[u][v]['weight'] * 15 for u, v in G.edges()]

    plt.figure(figsize=(7, 7))
    nx.draw_networkx_nodes(G, pos, node_size=1800,
                           node_color='lightsteelblue',
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos, font_size=11, font_weight='bold')
    nx.draw_networkx_edges(G, pos, width=weights,
                           arrowstyle='->', arrowsize=18)
    plt.title('Directional gene flow network')
    plt.axis('off')
    plt.tight_layout()
    save_figure(outfile)
    plt.close()


def draw_upset(df, outfile):
    memberships = []
    counts = []

    for _, row in df.iterrows():
        memberships.append(list(row['group_set']))
        counts.append(row['exclusive'])

    data = from_memberships(memberships, data=counts)
    upset = UpSet(data, subset_size='sum',
                  show_counts=True, sort_by='cardinality')

    plt.figure(figsize=(10, 6))
    upset.plot()
    plt.suptitle('Distribution of clusters across aquatic environments')
    plt.tight_layout()
    save_figure(outfile)
    plt.close()


def draw_net_flow_bar(net_flow, outfile):
    envs = list(net_flow.keys())
    values = [net_flow[e] for e in envs]

    plt.figure(figsize=(8, 4))
    bars = plt.bar(envs, values)
    plt.axhline(0, linewidth=1)

    for bar, val in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2,
                 val, f'{val:.2f}',
                 ha='center',
                 va='bottom' if val >= 0 else 'top')

    plt.ylabel('Net directional gene flow')
    plt.title('Net directional flow among aquatic environments')
    plt.tight_layout()
    save_figure(outfile)
    plt.close()


def get_dominant_direction_edges(direction_mat, threshold):
    edges = []
    envs = direction_mat.index.tolist()
    for i, src in enumerate(envs):
        for tgt in envs[i+1:]:
            val_st = direction_mat.loc[src, tgt]
            val_ts = direction_mat.loc[tgt, src]
            if max(val_st, val_ts) < threshold:
                continue
            edges.append((src, tgt, val_st) if val_st >= val_ts else (tgt, src, val_ts))
    return edges


def draw_global_direction_network(direction_mat, env_total, net_flow, threshold, outfile):
    G = nx.DiGraph()
    for env in direction_mat.index:
        G.add_node(env)

    edges = get_dominant_direction_edges(direction_mat, threshold)
    for src, tgt, val in edges:
        G.add_edge(src, tgt, weight=val)

    pos = nx.spring_layout(G, seed=42, weight='weight', k=2.2, iterations=200)

    # adjust positions for better visualization
    pos['NW'] += np.array([-0.50, -0.25])
    pos['MW'] += np.array([-0.70, 0])
    pos['GW'] += np.array([0.15, 0])

    richness = np.array(list(env_total.values()))
    min_r, max_r = richness.min(), richness.max()

    def norm(x):
        return (x - min_r) / (max_r - min_r) if max_r > min_r else 0.5

    node_sizes = [800 + norm(env_total[n]) * 2200 for n in G.nodes()]
    flows = np.array([net_flow[n] for n in G.nodes()])
    vmax = np.max(np.abs(flows))
    cmap = plt.cm.RdBu_r
    node_colors = cmap(plt.Normalize(-vmax, vmax)(flows))

    plt.figure(figsize=(8, 8))
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes,
                           node_color=node_colors,
                           edgecolors='black')
    nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')

    widths = [np.sqrt(d['weight']) * 12 for _, _, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, width=widths,
                           arrows=True, arrowstyle='-|>',
                           arrowsize=22,
                           connectionstyle='arc3,rad=0.0')

    labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels,
                                 bbox=dict(facecolor='white', alpha=0.7))

    sm = plt.cm.ScalarMappable(cmap=cmap,
                              norm=plt.Normalize(-vmax, vmax))
    sm.set_array([])
    plt.colorbar(sm, shrink=0.75,
                 label='Net directional gene flow (source → sink)')

    plt.title('Global directional gene flow among environments')
    plt.axis('off')
    plt.tight_layout()
    save_figure(outfile)
    plt.close()


# ===============================
# Main
# ===============================

def main():
    # Load occurrence data
    df, envs = load_data('env_combination_results.txt')
    env_total = calc_env_total(df, envs)

    # Jaccard similarity (all clusters)
    jac_all = calc_jaccard(df, envs, env_total)
    jac_all.to_csv('jaccard_all_matrix.tsv', sep='\t', float_format='%.4f')
    draw_heatmap(jac_all, 'Jaccard similarity (all clusters)', 'jaccard_all_heatmap')

    # Jaccard similarity (shared clusters only)
    df_shared = calc_shared_only_df(df)
    env_total_shared = calc_env_total(df_shared, envs)
    jac_shared = calc_jaccard(df_shared, envs, env_total_shared)
    jac_shared.to_csv('jaccard_shared_matrix.tsv', sep='\t', float_format='%.4f')
    draw_heatmap(jac_shared, 'Jaccard similarity (shared clusters only)', 'jaccard_shared_heatmap')

    # Directional gene flow
    direction_mat = calc_direction_gene_flow(df, envs, env_total)
    direction_mat.to_csv('directional_gene_flow_matrix.tsv', sep='\t', float_format='%.4f')
    draw_heatmap(direction_mat, 'Directional gene flow (source → target)', 'directional_gene_flow_heatmap')

    draw_direction_gene_flow_network(direction_mat, threshold=0.01, outfile='directional_gene_flow_network')

    # UpSet plots
    draw_upset(df, 'upset_all_clusters')
    df[['group', 'exclusive']].to_csv('upset_all_clusters.tsv', sep='\t', index=False)

    draw_upset(df_shared, 'upset_shared_clusters')
    df_shared[['group', 'exclusive']].to_csv('upset_shared_clusters.tsv', sep='\t', index=False)

    # Net directional flow
    net_flow = calc_net_flow(direction_mat)
    pd.DataFrame.from_dict(net_flow, orient='index',
                           columns=['net_directional_flow']) \
        .to_csv('net_flow_values.tsv', sep='\t', float_format='%.4f')

    draw_net_flow_bar(net_flow, 'net_flow_barplot')

    # Dominant directional edges
    edges = get_dominant_direction_edges(direction_mat, threshold=0.02)
    pd.DataFrame(edges, columns=['source', 'target', 'direction_gene_flow']) \
        .to_csv('global_flow_edges.tsv', sep='\t', index=False, float_format='%.4f')

    # Node attributes
    pd.DataFrame({
        'environment': envs,
        'cluster_richness': [env_total[e] for e in envs],
        'net_flow': [net_flow[e] for e in envs]
    }).to_csv('global_flow_nodes.tsv', sep='\t', index=False, float_format='%.4f')

    # Global directional gene flow network
    draw_global_direction_network(direction_mat, env_total, net_flow,
                                  threshold=0.02,
                                  outfile='global_directional_gene_flow_network')

    print('Analysis finished. Figures and data saved.')


if __name__ == '__main__':
    main()
