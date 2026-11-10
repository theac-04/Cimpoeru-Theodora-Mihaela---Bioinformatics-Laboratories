import sys
import os
from pathlib import Path

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(filepath: str) -> str:
    seq_parts = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()


def detect_transposons_sliding_window(
        seq: str,
        ir_len: int = 8,
        min_dist: int = 100,
        max_dist: int = None,
        tsd_len: int = 5
):
    found = set()
    L = len(seq)
    if max_dist is None:
        max_dist = L

    for s1 in range(tsd_len, L - min_dist - ir_len):
        s1_pattern = seq[s1: s1 + ir_len]
        if 'N' in s1_pattern:
            continue
        target_s2 = revcomp(s1_pattern)
        s2_search_start = s1 + min_dist
        s2_search_end = min(L - ir_len - tsd_len, s1 + max_dist)
        s2 = seq.find(target_s2, s2_search_start, s2_search_end)

        while s2 != -1:
            end_s2 = s2 + ir_len - 1
            if end_s2 + 1 + tsd_len <= L:
                left_tsd = seq[s1 - tsd_len: s1]
                right_tsd = seq[end_s2 + 1: end_s2 + 1 + tsd_len]

                if left_tsd == right_tsd and 'N' not in left_tsd:
                    ext = ir_len
                    while (
                            s1 + ext < s2 - ext and
                            seq[s1 + ext] == seq[end_s2 - ext].translate(COMP)
                    ):
                        ext += 1
                    found.add((s1, end_s2, ext))
            s2 = seq.find(target_s2, s2 + 1, s2_search_end)

    return sorted(found)


if __name__ == "__main__":
    REPO_ROOT = Path(__file__).resolve().parent

    print(f"Loading FASTA file: {REPO_ROOT / "Laboratory 6" / "GCF_003926195.1_ASM392619v1_genomic.fna"}...")
    genome_seq = read_fasta(str(REPO_ROOT / "Laboratory 6" / "GCF_003926195.1_ASM392619v1_genomic.fna"))
    seq_len = len(genome_seq)
    print(f"Genome loaded successfully! Total length: {seq_len:,} bp")
    print("Searching genome using s1 and s2 sliding windows...")
    candidates = detect_transposons_sliding_window