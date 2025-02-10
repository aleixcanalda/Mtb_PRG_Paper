from Bio import AlignIO
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the alignment from .aln file
alignment = AlignIO.read("Manila/dnaa_mycobacterium.aln", "clustal")

# Number of sequences in the alignment
num_sequences = len(alignment)
identity_matrix = np.zeros((num_sequences, num_sequences))

# Function to calculate pairwise identity
def calculate_identity(seq1, seq2):
    matches = 0
    total = 0
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':  # Ignore gaps
            total += 1
            if a == b:
                matches += 1
    return matches / total if total > 0 else 0

# Calculate pairwise identity for all sequence pairs
for i in range(num_sequences):
    for j in range(i, num_sequences):
        identity = calculate_identity(alignment[i].seq, alignment[j].seq)
        identity_matrix[i][j] = identity
        identity_matrix[j][i] = identity  # Symmetric matrix

# Plot the heatmap
sns.heatmap(
    identity_matrix,
    annot=True,
    cmap="YlGnBu",
    xticklabels=[record.id for record in alignment],
    yticklabels=[record.id for record in alignment]
)
plt.title("Pairwise Sequence Identity of Essential and Non-essential Genes")
plt.xlabel("Species")
plt.ylabel("Species")
plt.tight_layout()
plt.show()
