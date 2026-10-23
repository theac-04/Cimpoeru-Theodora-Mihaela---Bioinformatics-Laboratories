import random

MIN_OVERLAP = 10

def split_genome(genome, a, b, min_overlap=MIN_OVERLAP):
    if a <= min_overlap:
        raise ValueError("a must be > min_overlap so fragments can advance")
    frags, pos = [], 0
    while True:
        L = random.randint(a, b)
        frags.append(genome[pos:pos + L])
        if pos + L >= len(genome):
            break
        overlap = random.randint(min_overlap, L - 1)
        pos += L - overlap
    return frags

def overlap(x, y, min_len):
    for k in range(min(len(x), len(y)), min_len - 1, -1):
        if x.endswith(y[:k]):
            return k
    return 0

def reconstruct(frags, min_overlap=MIN_OVERLAP):
    frags = list(set(frags))
    frags = [f for f in frags if not any(f != g and f in g for g in frags)]
    while len(frags) > 1:
        best, bi, bj = 0, None, None
        for i, x in enumerate(frags):
            for j, y in enumerate(frags):
                if i != j:
                    k = overlap(x, y, min_overlap)
                    if k > best:
                        best, bi, bj = k, i, j
        if best < min_overlap:
            break  # cannot merge any further
        merged = frags[bi] + frags[bj][best:]
        frags = [f for n, f in enumerate(frags) if n not in (bi, bj)] + [merged]
    return max(frags, key=len)

if __name__ == "__main__":
    genome = "".join(random.choice("ACGT") for _ in range(200))
    a, b = 25, 40

    frags = split_genome(genome, a, b)
    random.shuffle(frags)
    result = reconstruct(frags)

    print(f"Interval [a, b]: [{a}, {b}]; min overlap: {MIN_OVERLAP}")
    print(f"Fragments ({len(frags)}), lengths: {[len(f) for f in frags]}\n")
    print("Original:", genome)
    print("Reconstructed:", result)
    print("\nMatch:", result == genome)