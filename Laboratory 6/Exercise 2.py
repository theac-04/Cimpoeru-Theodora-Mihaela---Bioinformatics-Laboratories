"""
Exercise 2: detect the positions (start, end) of transposable elements in a DNA sequence.
"""
import sys

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def detect(seq, min_ir=8, tsd_len=5, min_len=24, max_len=300):
    found = set()
    L = len(seq)
    for i in range(tsd_len, L - min_len):
        rc = revcomp(seq[i:i + min_ir])
        j = seq.find(rc, i + min_len - min_ir)
        while j != -1 and j - i + min_ir <= max_len:
            end = j + min_ir - 1
            if end + tsd_len < L and seq[i - tsd_len:i] == seq[end + 1:end + 1 + tsd_len]:
                n = min_ir
                while i + n < end - n and seq[i + n] == seq[end - n].translate(COMP):
                    n += 1
                found.add((i, end, n))
            j = seq.find(rc, j + 1)
    return sorted(found)


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "sequence.txt"
    with open(path) as f:
        seq = "".join(f.read().split()).upper()
    print(f"Sequence length: {len(seq)} b")

    results = detect(seq)
    print(f"\nDetected {len(results)} transposable element(s) (1-based):")
    for i, (s, e, n) in enumerate(results, 1):
        print(f"  E{i}: start={s + 1:4d}  end={e + 1:4d}  length={e - s + 1:3d}  IR={n}")

    try:
        with open("truth.txt") as f:
            truth = sorted(tuple(int(x) for x in line.split()) for line in f if line.strip())
        got = sorted((s + 1, e + 1) for s, e, _ in results)
        print("\nMatches ground truth:", got == truth)
    except FileNotFoundError:
        pass