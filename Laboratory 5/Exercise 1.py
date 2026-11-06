import random
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



def fetch_ncbi_sequence(accession="NM_000518", email="your.email@example.com"):
    from Bio import Entrez, SeqIO
    Entrez.email = email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as h:
        return str(SeqIO.read(h, "fasta").seq).upper()


def get_sequence(accession="NM_000546", fallback_seed=1):
    try:
        seq = fetch_ncbi_sequence(accession)
        src = f"NCBI {accession}"
    except Exception as e:  # offline -> random stand-in sequence
        rnd = random.Random(fallback_seed)
        seq = "".join(rnd.choice("ACGT") for _ in range(2000))
        src = f"random fallback (NCBI unavailable: {type(e).__name__})"
    if not 1000 <= len(seq) <= 3000:
        raise ValueError(f"Sequence length {len(seq)} outside 1000-3000")
    print(f"Sequence source: {src}, length = {len(seq)} nt")
    return seq


def take_samples(seq, n=10, lo=100, hi=3000):
    samples = []
    for _ in range(n):
        length = random.randint(lo, min(hi, len(seq)))
        start = random.randint(0, len(seq) - length)
        samples.append(seq[start:start + length])
    return samples


def migration_distance(length, gel_len=100.0, max_len=3000, min_len=100):
    import math
    lo, hi = math.log10(min_len), math.log10(max_len)
    frac = (math.log10(length) - lo) / (hi - lo)
    return gel_len * (0.95 - 0.80 * frac)


def simulate_gel(samples, ladder=(100, 250, 500, 1000, 1500, 2000, 3000)):
    fig, ax = plt.subplots(figsize=(9, 6), facecolor="#202020")
    ax.set_facecolor("#1a3a5a")
    lanes = [("Ladder", list(ladder))] + [(f"S{i+1}", [len(s)]) for i, s in enumerate(samples)]
    for x, (name, sizes) in enumerate(lanes):
        ax.add_patch(plt.Rectangle((x - 0.35, 99), 0.7, 3, color="black"))   # well
        for size in sizes:
            y = migration_distance(size)
            ax.add_patch(plt.Rectangle((x - 0.3, y - 0.8), 0.6, 1.6,
                color="#9fe8ff" if name != "Ladder" else "white", alpha=.9))
            if name == "Ladder":
                ax.text(x - 0.45, y, f"{size}", ha="right", va="center", color="white", fontsize=7)
    ax.set_xlim(-1.2, len(lanes) - 0.4)
    ax.set_ylim(100, 0)
    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels([n for n, _ in lanes], color="white")
    ax.set_ylabel("Migration distance (mm)", color="white")
    ax.tick_params(colors="white")
    ax.set_title("Simulated agarose gel electrophoresis", color="white")
    fig.savefig("gel.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print("Gel image saved to gel.png")


def exercise1(seq):
    samples = take_samples(seq)
    print("\n--- Exercise 1 ---")
    print(f"{'Sample':<8}{'Length (nt)':>12}{'Migration (mm)':>16}")
    for i, s in enumerate(samples, 1):
        print(f"S{i:<7}{len(s):>12}{migration_distance(len(s)):>16.1f}")
    simulate_gel(samples)


def find_repeats(seq, unit_min=3, unit_max=5, total_min=6, total_max=10):
    found = []
    for u in range(unit_min, unit_max + 1):
        pattern = re.compile(rf"(?=((.{{{u}}})\2+))")
        i = 0
        while i <= len(seq) - 2 * u:
            m = pattern.match(seq, i)
            if m:
                run = m.group(1)
                unit = m.group(2)
                primitive = not any(u % d == 0 and unit == unit[:d] * (u // d) for d in range(1, u))
                if primitive and total_min <= len(run) <= total_max:
                    found.append((i, unit, len(run) // u, run))
                i += max(len(run) - u + 1, 1)
            else:
                i += 1
    return sorted(found)


def exercise2(seq):
    print("\n--- Exercise 2 ---")
    reps = find_repeats(seq)
    print(f"{'Pos':>6}  {'Unit':<6}{'Copies':>7}  Repeat")
    for pos, unit, k, run in reps:
        print(f"{pos+1:>6}  {unit:<6}{k:>7}  {run}")
    print(f"Total repeats found: {len(reps)}")


if __name__ == "__main__":
    random.seed(42)
    sequence = get_sequence("NM_000546")
    exercise1(sequence)
    exercise2(sequence)