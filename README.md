# Analysis of Microbiomes in Water Systems

This repository contains code and data for analyzing microbiomes in water systems.

## 1. Gene Catalogue Construction

### 1.1 `make_protein_clusters.py`

Process MMseqs-generated data.  
This script processes a large protein clustering table and writes cluster memberships by grouping sequences in a memory-efficient chunk-wise manner.

**Input**
- `protein_cluster.tsv`

**Output**
- `protein_clusters.txt`

---

### 1.2 `extract_rep_nucl_sequence.py`

This script extracts representative sequence IDs from an MMseqs-generated representative protein FASTA file and retrieves the corresponding nucleotide sequences from a reference FASTA file.

**Input**
- `protein_rep_seq.fasta`
- `nucl.fa`

**Output**
- `nucl_rep_extract.fa`

---

### 1.3 `process_cluster_info.py`

This script parses cluster memberships and metadata to generate environmental-level sample counts and study counts for each cluster.

**Input**
- `protein_clusters.txt`
- `env.txt`
- `study.txt`

**Output**
- `cluster_env_sample_number_matrix.txt`
- `cluster_env_study_number_matrix.txt`

---

### 1.4 `run_mmseqs_selfcompare.sh`

This script performs self-comparison for each FASTA file using MMseqs2 and exports pairwise alignments into TSV format.

**Input**
- `clusters30/*.fasta`

**Output**
- `results/*.tsv`

---

### 1.5 `make_rarefaction_curves.py`

This script computes rarefaction curves of protein clusters by progressively adding samples.  
It generates a global curve and environment-specific curves.

**Input**
- `protein_cluster.tsv`
- `env.txt`

**Output**
- `rarecurve.txt`
- `rarecurve-AS.txt`
- `rarecurve-DW.txt`
- `rarecurve-GW.txt`
- `rarecurve-MW.txt`
- `rarecurve-NW.txt`
- `rarecurve-WW.txt`

---

### 1.6 `summarize_env_combination.py`

This script summarizes environment co-occurrence patterns of clusters based on study-level presence and reports both strict and inclusive combination counts.

**Input**
- `cluster_env_study_number_matrix.txt`

**Output**
- `env_combination_results.txt`

---

### 1.7 `jaccard_anlysis.py`

This script computes Jaccard similarity, directional gene flow, UpSet plots, and net directional flow networks among aquatic environments based on cluster occurrence data.

**Input**
- `env_combination_results.txt`

**Output**
- `jaccard_all_matrix.tsv`
- `jaccard_all_heatmap.png / .pdf`
- `jaccard_shared_matrix.tsv`
- `jaccard_shared_heatmap.png / .pdf`
- `directional_gene_flow_matrix.tsv`
- `directional_gene_flow_heatmap.png / .pdf`
- `directional_gene_flow_network.png / .pdf`
- `upset_all_clusters.tsv`
- `upset_all_clusters.png / .pdf`
- `upset_shared_clusters.tsv`
- `upset_shared_clusters.png / .pdf`
- `net_flow_values.tsv`
- `net_flow_barplot.png / .pdf`
- `global_flow_edges.tsv`
- `global_flow_nodes.tsv`
- `global_directional_gene_flow_network.png / .pdf`

---

## 2. IPG Identification

### 2.1 `filter_protein_clusters.py`

This script filters protein clusters that meet threshold criteria for each environment and generates six output files corresponding to 10%–60% thresholds.

**Input**
- `cluster_env_sample_number_matrix.txt`

**Output**
- `IPG_4group_10%.txt`
- `IPG_4group_20%.txt`
- `IPG_4group_30%.txt`
- `IPG_4group_40%.txt`
- `IPG_4group_50%.txt`
- `IPG_4group_60%.txt`

---

### 2.2 `extract_ipg.py`

This script extracts nucleotide sequences corresponding to IPG clusters for different thresholds (10%–60%).

**Input**
- `nucl_rep_extract.fa`
- `IPG_4group_10%.txt`
- `IPG_4group_20%.txt`
- `IPG_4group_30%.txt`
- `IPG_4group_40%.txt`
- `IPG_4group_50%.txt`
- `IPG_4group_60%.txt`

**Output**
- `IPG_4group_10%.fa`
- `IPG_4group_20%.fa`
- `IPG_4group_30%.fa`
- `IPG_4group_40%.fa`
- `IPG_4group_50%.fa`
- `IPG_4group_60%.fa`

---

### 2.3 `calculate_average_percent.py`

This script calculates the average proportion of sequences in each group above varying thresholds and generates a TSV file with average proportions and a line plot for visualization.

**Input**
- `IPG_4group_10%.txt`

**Output**
- `cluster.txt`
- `cluster_plot.png / .pdf`
