"""
GWAS Visualization: Manhattan Plot & QQ Plot
=============================================
Generates publication-standard GWAS figures from synthetic summary statistics.
Demonstrates genomics data visualization skills central to population genetics work.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

SEED = 42
N_SNPS = 50_000
RESULTS_DIR = "results"

GENOME_WIDE_THRESH = 5e-8
SUGGESTIVE_THRESH = 1e-5

CHR_SIZES_MB = {
    1: 249, 2: 243, 3: 198, 4: 191, 5: 181, 6: 171, 7: 159, 8: 146,
    9: 141, 10: 136, 11: 135, 12: 133, 13: 115, 14: 107, 15: 102,
    16: 91,  17: 81,  18: 78,  19: 59,  20: 63,  21: 47,  22: 51,
}

rng = np.random.default_rng(SEED)


def simulate_gwas_summary(n_snps, signal_loci=None):
    """
    Simulate GWAS summary statistics with realistic null p-value distribution
    and injected association signals at specified loci.
    """
    chrs = list(CHR_SIZES_MB.keys())
    weights = np.array([CHR_SIZES_MB[c] for c in chrs], dtype=float)
    weights /= weights.sum()
    chr_assignments = rng.choice(chrs, size=n_snps, p=weights)

    bp_positions = np.array([
        rng.integers(1, CHR_SIZES_MB[c] * 1_000_000)
        for c in chr_assignments
    ])

    pvalues = rng.uniform(0, 1, size=n_snps)
    pvalues = pvalues ** 1.01   # mild inflation (λ ~ 1.02)

    snp_ids = [f"rs{rng.integers(1_000_000, 9_999_999)}" for _ in range(n_snps)]

    df = pd.DataFrame({
        "CHR": chr_assignments,
        "BP": bp_positions,
        "SNP": snp_ids,
        "P": pvalues,
    })

    if signal_loci:
        for chrom, center_mb, n_signals, min_logp, max_logp in signal_loci:
            center_bp = center_mb * 1_000_000
            window = 500_000
            mask = (df["CHR"] == chrom) & (df["BP"].between(center_bp - window, center_bp + window))
            n_in_window = mask.sum()
            if n_in_window == 0:
                continue
            logp_vals = rng.uniform(min_logp, max_logp, size=n_in_window)
            logp_vals = np.sort(logp_vals)[::-1]
            logp_vals[n_signals:] = rng.uniform(0, 3, size=n_in_window - n_signals)
            df.loc[mask, "P"] = 10 ** (-logp_vals)

    return df.sort_values(["CHR", "BP"]).reset_index(drop=True)


def compute_cumulative_positions(df):
    """Add genome-wide cumulative base-pair positions for Manhattan plot x-axis."""
    chrom_offsets = {}
    offset = 0
    for chrom in sorted(df["CHR"].unique()):
        chrom_offsets[chrom] = offset
        offset += CHR_SIZES_MB.get(chrom, 50) * 1_000_000 + 5_000_000

    df = df.copy()
    df["cumBP"] = df.apply(lambda r: r["BP"] + chrom_offsets[r["CHR"]], axis=1)
    df["neg_log10_p"] = -np.log10(df["P"].clip(lower=1e-35))
    return df, chrom_offsets


def plot_manhattan(df, chrom_offsets, output_path):
    palette = ["#2166ac", "#b2182b"]
    chrs_sorted = sorted(df["CHR"].unique())

    fig, ax = plt.subplots(figsize=(14, 5))

    xtick_positions = []
    xtick_labels = []

    for i, chrom in enumerate(chrs_sorted):
        sub = df[df["CHR"] == chrom]
        color = palette[i % 2]
        ax.scatter(sub["cumBP"], sub["neg_log10_p"],
                   c=color, s=3, alpha=0.6, linewidths=0, rasterized=True)

        mid = chrom_offsets[chrom] + CHR_SIZES_MB[chrom] * 500_000
        xtick_positions.append(mid)
        xtick_labels.append(str(chrom))

    ax.axhline(-np.log10(GENOME_WIDE_THRESH), color="red", lw=1.2,
               linestyle="--", label=f"Genome-wide (p={GENOME_WIDE_THRESH:.0e})")
    ax.axhline(-np.log10(SUGGESTIVE_THRESH), color="orange", lw=1,
               linestyle=":", label=f"Suggestive (p={SUGGESTIVE_THRESH:.0e})")

    sig = df[df["P"] < GENOME_WIDE_THRESH]
    if not sig.empty:
        peaks = sig.groupby("CHR").apply(lambda x: x.loc[x["neg_log10_p"].idxmax()])
        for _, row in peaks.iterrows():
            ax.annotate(row["SNP"], xy=(row["cumBP"], row["neg_log10_p"]),
                        fontsize=6.5, ha="center", va="bottom",
                        xytext=(0, 4), textcoords="offset points", color="darkred")

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=7)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("−log₁₀(p)", fontsize=12)
    ax.set_title("GWAS Manhattan Plot", fontsize=14, fontweight="bold")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_qq(df, output_path):
    """QQ plot with genomic inflation factor λ."""
    observed = np.sort(-np.log10(df["P"].clip(lower=1e-35)))[::-1]
    n = len(observed)
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))

    median_chi2_obs = np.median((-np.log10(df["P"])) * 2 * np.log(10))
    lambda_gc = median_chi2_obs / 0.4549  # chi2(1) median

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(expected, observed, s=4, alpha=0.5, color="#2166ac", rasterized=True)
    max_val = max(observed.max(), expected.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], "r--", lw=1.2, label="Expected (null)")
    ax.set_xlabel("Expected −log₁₀(p)", fontsize=12)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=12)
    ax.set_title(f"QQ Plot  (λ = {lambda_gc:.3f})", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    signal_loci = [
        (3,  95,  5,  8.5, 12.0),
        (7,  120, 3,  7.5, 10.0),
        (11, 65,  4,  8.0, 11.5),
        (17, 44,  2,  6.5,  9.0),
    ]

    print(f"Simulating GWAS summary statistics ({N_SNPS:,} SNPs)...")
    gwas = simulate_gwas_summary(N_SNPS, signal_loci=signal_loci)
    n_sig = (gwas["P"] < GENOME_WIDE_THRESH).sum()
    print(f"  Genome-wide significant SNPs: {n_sig}")

    gwas_plot, chrom_offsets = compute_cumulative_positions(gwas)

    print("\nPlotting Manhattan plot...")
    plot_manhattan(gwas_plot, chrom_offsets, os.path.join(RESULTS_DIR, "manhattan_plot.png"))

    print("Plotting QQ plot...")
    plot_qq(gwas, os.path.join(RESULTS_DIR, "qq_plot.png"))

    print("\nDone.")


if __name__ == "__main__":
    main()
