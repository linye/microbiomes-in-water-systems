"""
This script extracts representative sequence IDs from an MMseqs protein FASTA file
and retrieves the corresponding nucleotide sequences from a reference FASTA file.
"""

from Bio import SeqIO

# Input and output files
PROTEIN_REP_SEQ = "protein_rep_seq.fasta"
NUCL_FILE = "nucl.fa"
OUTPUT_FILE = "nucl_rep_extract.fa"

# Step 1: Extract target IDs from protein_rep_seq.fasta
target_ids = set()
with open(PROTEIN_REP_SEQ, 'r') as fasta_file:
    for line in fasta_file:
        if line.startswith('>'):
            # Keep the part before the last underscore
            trimmed_id = line[1:].strip().rsplit('_', 1)[0]
            target_ids.add(trimmed_id)

print(f"Successfully extracted {len(target_ids)} target IDs")

# Step 2: Filter sequences from nucl.fa
matched_count = 0
with open(OUTPUT_FILE, 'w') as out_f:
    for record in SeqIO.parse(NUCL_FILE, "fasta"):
        # Keep the part before the last underscore in the header
        nucl_id = record.description.rsplit('_', 1)[0]

        if nucl_id in target_ids:
            matched_count += 1
            SeqIO.write(record, out_f, "fasta")

print(f"Extraction completed, total matched sequences: {matched_count}")
print(f"Results saved to: {OUTPUT_FILE}")
