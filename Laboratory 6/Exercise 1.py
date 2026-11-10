"""
Exercise 1: build an artificial DNA sequence (200-400 b) simulating 4 transposons.

"""

import random

COMP = str.maketrans("ACGT", "TGCA")
IR_LEN = 10
TSD_LEN = 5


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_transposon(body_len, rng):
    ir = rand_dna(IR_LEN, rng)
    return ir + rand_dna(body_len, rng) + revcomp(ir)


def insert(seq, elements, tr, pos):
    tsd = seq[pos:pos + TSD_LEN]
    new_seq = seq[:pos] + tsd + tr + seq[pos:]
    shift = len(tsd) + len(tr)
    start = pos + len(tsd)
    updated = []
    for s, e in elements:
        if s >= pos + TSD_LEN:
            s += shift
            e += shift
        elif e >= pos:
            e += shift
        updated.append((s, e))
    updated.append((start, start + len(tr) - 1))
    return new_seq, updated


def build_sequence(seed=7):
    rng = random.Random(seed)
    seq = rand_dna(90, rng)
    elems = []
    seq, elems = insert(seq, elems, make_transposon(40, rng), 30)
    seq, elems = insert(seq, elems, make_transposon(30, rng), len(seq) - 20)
    s, _ = elems[0]
    seq, elems = insert(seq, elems, make_transposon(18, rng), s + IR_LEN + 15)
    s, _ = elems[1]
    seq, elems = insert(seq, elems, make_transposon(14, rng), s + IR_LEN + 10)
    return seq, sorted(elems)


def show(seq, elems):
    marks = [" "] * len(seq)
    for idx, (s, e) in enumerate(elems, 1):
        for p in range(s, e + 1):
            marks[p] = str(idx) if marks[p] == " " else "*"   # * = overlap
    for k in range(0, len(seq), 60):
        print(f"{k + 1:4d} {seq[k:k + 60]}")
        print(f"     {''.join(marks[k:k + 60])}")


if __name__ == "__main__":
    seq, truth = build_sequence()
    print(f"Sequence length: {len(seq)} b\n")
    show(seq, truth)
    print("\nTrue positions (1-based start, end):")
    for i, (s, e) in enumerate(truth, 1):
        print(f"  T{i}: {s + 1:4d} - {e + 1:4d}  (len {e - s + 1})")

    with open("sequence.txt", "w") as f:
        f.write(seq + "\n")
    with open("truth.txt", "w") as f:
        for s, e in truth:
            f.write(f"{s + 1} {e + 1}\n")
    print("\nSaved sequence.txt and truth.txt")