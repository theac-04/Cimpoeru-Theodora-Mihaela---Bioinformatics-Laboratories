import math
from collections import Counter
import matplotlib.pyplot as plt

from ex_1 import SEQUENCE, ENZYMES, find_cut_positions


def nucleotide_frequencies(seq):
    counts = Counter(seq)
    return {b: counts[b] / len(seq) for b in "ACGT"}


def motif_probability(motif, freq):
    p = 1.0
    for base in motif:
        p *= freq[base]
    return p


def poisson_at_least(k, lam):
    return 1 - sum(math.exp(-lam) * lam**i / math.factorial(i) for i in range(k))


def plot_results(seq, freq, results):
    L = len(seq)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    max_j = 8
    width = 0.8 / len(results)
    for i, r in enumerate(results):
        xs = [j + (i - (len(results) - 1) / 2) * width for j in range(max_j + 1)]
        ys = [math.exp(-r["lam"]) * r["lam"] ** j / math.factorial(j)
              for j in range(max_j + 1)]
        bars = ax1.bar(xs, ys, width=width, label=f"{r['name']} (λ={r['lam']:.2f})")
        ax1.plot(r["observed"] + (i - (len(results) - 1) / 2) * width,
                 ys[r["observed"]] if r["observed"] <= max_j else 0,
                 "k*", markersize=12)
    ax1.plot([], [], "k*", markersize=12, label="observed count")
    ax1.set_xlabel("Number of sites in the sequence")
    ax1.set_ylabel("Probability")
    ax1.set_title(f"Poisson model, L = {L} bp")
    ax1.set_xticks(range(max_j + 1))
    ax1.legend()

    lengths = list(range(100, 20001, 100))
    for r in results:
        probs = [1 - math.exp(-(n - len(r["site"]) + 1) * r["p"]) for n in lengths]
        ax2.plot(lengths, probs, label=r["name"])
        ax2.plot(L, 1 - math.exp(-r["lam"]), "ko")
    ax2.axvline(L, color="gray", linestyle="--", label=f"this sequence ({L} bp)")
    ax2.axhline(0.5, color="lightgray", linestyle=":")
    ax2.set_xlabel("Sequence length (bp)")
    ax2.set_ylabel("P(at least one site)")
    ax2.set_title("Chance of at least one site vs length")
    ax2.legend()

    plt.tight_layout()
    plt.savefig("restriction_probability.png", dpi=150)
    plt.show()


def main():
    seq = SEQUENCE.upper()
    L = len(seq)
    freq = nucleotide_frequencies(seq)
    uniform = {b: 0.25 for b in "ACGT"}

    print(f"Length: {L} bp")
    print("Nucleotide frequencies:", {b: round(f, 4) for b, f in freq.items()}, "\n")
    results = []

    for enz in ENZYMES:
        motif = enz["site"]
        n_pos = L - len(motif) + 1                 # possible start positions
        p_uniform = motif_probability(motif, uniform)
        p = motif_probability(motif, freq)
        lam = n_pos * p                            # expected number of sites
        observed = len(find_cut_positions(seq, enz))
        results.append({"name": enz["name"], "site": motif, "p": p,
                        "lam": lam, "observed": observed})

        print(f"{enz['name']} ({motif})")
        print(f"  P(site at one position), uniform  : {p_uniform:.6f}  (1 in {1/p_uniform:.0f})")
        print(f"  P(site at one position), observed : {p:.6f}  (1 in {1/p:.0f})")
        print(f"  expected sites in {L} bp           : {lam:.3f}")
        print(f"  P(no site)  = e^-lambda            : {math.exp(-lam):.4f}")
        print(f"  P(>=1 site) = 1 - e^-lambda        : {1 - math.exp(-lam):.4f}")
        print(f"  observed sites                     : {observed}")
        print(f"  P(>= {observed} sites by chance)        : {poisson_at_least(observed, lam):.5f}\n")

    plot_results(seq, freq, results)


if __name__ == "__main__":
    main()