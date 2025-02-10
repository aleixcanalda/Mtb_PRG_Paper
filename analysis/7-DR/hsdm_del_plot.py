from Bio.Seq import Seq
import matplotlib.pyplot as plt

# Define the sequences
reference_seq = "atgccgcccaggaagaagcaggcgccgcaggcgccgtcgacgatgaaggagctcaaagacacgctctggaaggccgccgacaagctgcgcgggtcgctgtcggccagccaatacaaggacgtgatcctcggcctggtgttccttaagtacgtgtccgacgcgtatgacgaacggcgcgaggcaatccgtgccgagttggcggccgaaggaatggaggagtctcagatagaagacctgatcgacgatcccgagcagtaccagggttacggcgtattcgtcgtgccggtgagtgcgcgctggaagttcttggcagagaacacaaaaggcaagccagccgttggtggtgagccggcgaagaacatcggtcagctgatcgacgaggcgatggacgcggtaatgaaggccaatccaacactcggtgggacgctgccgaggctgtataacaaggacaacatcgaccagcgccggctcggtgagctgatcgacctatttaacagtgcgcgcttcagccggcagggcgagcaccgcgcccgggatctgatgggtgaggtctacgaatacttcctcggcaatttcgctcgcgcggaagggaagcggggtggcgagttctttaccccgcccagcgtggtcaaggtgatcgtggaggtgctggagccgtcgagtgggcgggtgtatgacccgtgctgcggttccggaggcatgtttgtgcagaccgagaagttcatctacgaacacgacggcgatccgaaggatgtctcgatctatggccaggaaagcattgaggagacctggcggatggcgaagatgaacctcgccatccacggcatcgacaacaaggggctcggcgcccgatggagtgataccttcgcccgcgaccagcacccggacgtgcagatggactacgtgatggccaatctgccgttcaacatcaaagactgggcccgcaacgaggaagacccacgctggcgcttcggtgttccgcccgccaataacgccaactacgcatggattcagcacatcctgtacaagttggcgccgggaggtcgggcgggcgtggtgatggccaacgggtcgatgtcgtcgaactccaacggcgagggggatattcgcgcgcagatcgtggaggcggatttggtttcctgcatggtcgcgttacccacccagctgttccgcagcaccggaatcccggtgtgcctgtggtttttcgccaaagacaaggcggcaggtaagcaagggtctatcgaccggtgcgggcaggtgctgttcatcgacgctcgtgaactgggcgacctagtggaccgggccgagcgggcgctgaccaacgaggagatcgtccgcatcggggataccttccacgcgtggcgcgggtcgaagtcggctgccgtcaaagggattatgtacgaggatgttccggggttctgtaagtcggcgacgttggcggagatcaaggcgaccgactatgcgctcacgccggggcggtatgtgggtacgcccgcggtcgaggacgacggagagccgatcgacgagaagatggcccggttgtcgaaggcgttgctggaggcgttcgatgagtcggcgcggctggaaagggtcgtgcgggagcaactggggcggcttagatga"
reference_seq=reference_seq.upper()
deletion_seq = reference_seq[:189] + reference_seq[291:]

# Translate both sequences
ref_protein = Seq(reference_seq).translate()
del_protein = Seq(deletion_seq).translate()

# Analyze frameshift
deletion_length = len(reference_seq) - len(deletion_seq)
is_frameshift = deletion_length % 3 != 0

# Find stop codons
def find_stop_positions(seq):
    protein = seq.translate()
    stops = [i for i, aa in enumerate(protein) if aa == "*"]
    return stops

ref_stops = find_stop_positions(Seq(reference_seq))
del_stops = find_stop_positions(Seq(deletion_seq))

# Determine if there's a new stop codon
new_stop_positions = [pos for pos in del_stops if pos not in ref_stops]

# Plot sequences with stops
fig, ax = plt.subplots(figsize=(12, 3))

# Plot reference sequence
ax.plot(range(len(reference_seq)), [1] * len(reference_seq), 'o-', label="Reference")
for stop in ref_stops:
    ax.text(stop * 3, 1.1, "*", color="blue", ha="center", fontsize=12)

# Plot deleted sequence
ax.plot(range(len(deletion_seq)), [0] * len(deletion_seq), 'o-', label="With Deletion")
for stop in del_stops:
    color = "red" if stop in new_stop_positions else "blue"
    ax.text(stop * 3, -0.1, "*", color=color, ha="center", fontsize=12)

# Highlight deletion
# Highlight deletion region in the reference sequence
ax.axhline(1, xmin=248  / len(reference_seq), xmax=334/ len(reference_seq),color="red", linewidth=3, label="Deleted Region")

# Add dashed lines to connect deletion
ax.plot([188, 188], [1, 0], 'k--', linewidth=1)  # From start of deletion
ax.plot([292, 188], [1, 0], 'k--', linewidth=1)  # From end of deletion

# Formatting
ax.set_yticks([0, 1])
ax.set_yticklabels(["hsdM 102bp DEL", "Reference"])
ax.set_xlabel("Nucleotide Position")
ax.legend()
plt.title("hsdM DEL")

# Print analysis
print("Reference Protein:", ref_protein)
print("Deleted Protein  :", del_protein)
print("Frameshift?" if is_frameshift else "No Frameshift")
if new_stop_positions:
    print(f"New premature stop codons introduced at amino acid positions: {new_stop_positions}")
else:
    print("No new premature stop codons introduced.")

plt.show()
