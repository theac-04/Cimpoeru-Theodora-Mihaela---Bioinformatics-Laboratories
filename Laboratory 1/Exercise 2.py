import matplotlib.pyplot as plt

def count_combinations(sequence: str, length: int) -> dict:
    counts = {}

    for i in range(len(sequence) - length + 1):
        combo = sequence[i : i + length]
        if combo in counts:
            counts[combo] += 1
        else:
            counts[combo] = 1

    return dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))


def plot_side_by_side(freq_2: dict, freq_3: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].bar(list(freq_2.keys()), list(freq_2.values()), color="#2b5c8f")
    axes[0].set_title("2-Letter Combinations (Descending)")
    axes[0].set_xlabel("Dinucleotides")
    axes[0].set_ylabel("Frequency")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].grid(axis="y", linestyle="--", alpha=0.7)

    axes[1].bar(list(freq_3.keys()), list(freq_3.values()), color="#d95f02")
    axes[1].set_title("3-Letter Combinations (Descending)")
    axes[1].set_xlabel("Trinucleotides")
    axes[1].set_ylabel("Frequency")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].grid(axis="y", linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    sequence = "ATGCGATCGATCGATCGATC"

    freq_2 = count_combinations(sequence, 2)
    freq_3 = count_combinations(sequence, 3)

    print("2-Letter counts:", freq_2)
    print("3-Letter counts:", freq_3)

    plot_side_by_side(freq_2, freq_3)