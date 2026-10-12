s = "TCGATGACCTGCAGCTAC"
s = s.replace('T', 'U')

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

# Invert dictionary to map each codon directly to (name, abvr)
GENETIC_CODE = {
    codon: (name, abvr)
    for (name, abvr), codons in aminoacids.items()
    for codon in codons
}

# Split into 3-letter triplets
arr = [s[i:i+3] for i in range(0, len(s), 3)]

names = []
abvrs = []
translating = False  # Track whether translation has started

for codon in arr:
    if codon in GENETIC_CODE:
        name, abvr = GENETIC_CODE[codon]

        # 1. Start translating when Methionine (AUG) is found
        if not translating and abvr == "M":
            translating = True

        # 2. If actively translating, check for stop or append amino acid
        if translating:
            if abvr == "*":
                break  # Stop immediately on the first stop codon
            names.append(name)
            abvrs.append(abvr)

print("Full names:", names)
print("Protein sequence:", "".join(abvrs))