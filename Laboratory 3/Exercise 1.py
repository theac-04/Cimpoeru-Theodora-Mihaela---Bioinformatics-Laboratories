# 1. design a software application which is able to read a FASTA file of 1 genome from either a virus or a bacteria.
# measure the frequencies of the codons from the sequence of the genome.
# in order to measure the frequency, use sliding a window of 3 positions to slide over the genome & counting the
# frequency of each codon seen in the sliding window.
# the application must show top 10 highest frequency codons. also your application must show, for each top 10 codons,
# the associated aminoacid.
import io
import numpy as np
import matplotlib.pyplot as plt
import contextlib

FASTA_PATH = r"C:\Users\THEA\Desktop\poli 2023-2027\facultate\anul 4 (2026 - 2027)\sem 1\Bioinformatics\labs\lab 3\GCF_002118885.1_ViralProj385635_genomic.fna"


def read_fasta(path):
    """Read a FASTA file. Returns (header, sequence) with the sequence as RNA letters (A, C, G, U)."""
    header = None
    parts = io.StringIO()
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
            else:
                parts.write(line)
    return header, parts.getvalue().upper().replace("T", "U")  # the codon table uses U instead of T


# read the genome from the FASTA file
header, s = read_fasta(FASTA_PATH)

aminoacids = {
    ("Alanine", "A"): ["GCU", "GCC", "GCA", "GCG"],
    ("Cysteine", "C"): ["UGU", "UGC"],
    ("Aspartic acid", "D"): ["GAU", "GAC"],
    ("Glutamic acid", "E"): ["GAA", "GAG"],
    ("Phenylalanine", "F"): ["UUU", "UUC"],
    ("Glycine", "G"): ["GGU", "GGC", "GGA", "GGG"],
    ("Histidine", "H"): ["CAU", "CAC"],
    ("Isoleucine", "I"): ["AUU", "AUC", "AUA"],
    ("Lysine", "K"): ["AAA", "AAG"],
    ("Leucine", "L"): ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    ("Methionine", "M"): ["AUG"],
    ("Asparagine", "N"): ["AAU", "AAC"],
    ("Proline", "P"): ["CCU", "CCC", "CCA", "CCG"],
    ("Glutamine", "Q"): ["CAA", "CAG"],
    ("Arginine", "R"): ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    ("Serine", "S"): ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    ("Threonine", "T"): ["ACU", "ACC", "ACA", "ACG"],
    ("Valine", "V"): ["GUU", "GUC", "GUA", "GUG"],
    ("Tryptophan", "W"): ["UGG"],
    ("Tyrosine", "Y"): ["UAU", "UAC"],
    ("Stop", "*"): ["UAA", "UAG", "UGA"],
}

# codon -> (amino acid name, one-letter code)
codon_to_aa = {codon: aa for aa, codons in aminoacids.items() for codon in codons}

# fixed ordering of the 64 codons, so the counts can live in a numpy array
all_codons = sorted(codon_to_aa)
index_of = {codon: i for i, codon in enumerate(all_codons)}
counts = np.zeros(len(all_codons), dtype=np.int64)

# slide a window of 3 positions, one base at a time, over the whole genome
skipped = 0
for i in range(len(s) - 2):
    idx = index_of.get(s[i:i + 3])
    if idx is None:  # window contains something other than A/C/G/U (e.g. N)
        skipped += 1
    else:
        counts[idx] += 1

total = int(counts.sum())

# top 10 codons (stable sort: ties keep alphabetical order)
top_idx = np.argsort(-counts, kind="stable")[:10]

print(f"File   : {FASTA_PATH}")
print(f"Record : {header}")
print(f"Length : {len(s):,} nt | windows counted: {total:,} | windows skipped (non-ACGU): {skipped:,}\n")
print(f"{'Rank':<5}{'Codon':<7}{'Count':>10}{'Freq %':>9}   Amino acid")
print("-" * 55)
labels, values = [], []
for rank, i in enumerate(top_idx, start=1):
    codon = all_codons[i]
    name, letter = codon_to_aa[codon]
    pct = 100 * counts[i] / total
    print(f"{rank:<5}{codon:<7}{counts[i]:>10,}{pct:>8.2f}%   {name} ({letter})")
    labels.append(f"{codon}\n{name}\n({letter})")
    values.append(pct)

# bar chart of the top 10
plt.figure(figsize=(12, 6))
plt.bar(range(10), values, color="steelblue")
plt.xticks(range(10), labels, fontsize=8)
plt.ylabel("Frequency (% of all windows)")
plt.title(f"Top 10 codons by sliding window\n{header}", fontsize=10)
plt.tight_layout()
with contextlib.suppress(Exception):  # don't crash if the IDE has no GUI backend
    plt.show()