"""
RNA-seq Differential Expression Analysis
=========================================
Synthetic count data pipeline: normalization → statistical testing → volcano plot.
Demonstrates core DE analysis concepts used in real-world genomics pipelines.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

SEED = 42
N_GENES = 500
N_SAMPLES_PER_GROUP = 10
N_DE_GENES = 60          # genes with true differential expression
RESULTS_DIR = "results"

rng = np.random.default_rng(SEED)


def simulate_counts(n_genes, n_samples, de_gene_indices, de_fold_changes):
    """Simulate negative-binomial-like count data for two conditions."""
    base_means = rng.lognormal(mean=5, sigma=1.5, size=n_genes).clip(1, None)
    dispersion = rng.uniform(0.05, 0.4, size=n_genes)

    def nb_counts(means, n):
        r = 1.0 / dispersion          # (n_genes,)
        p = r / (r + means)           # (n_genes,) element-wise
        counts = rng.negative_binomial(
            np.tile(r[:, None], (1, n)),
            np.tile(p[:, None], (1, n)),
        )
        return counts.astype(float)

    control = nb_counts(base_means, n_samples)

    treatment_means = base_means.copy()
    treatment_means[de_gene_indices] *= de_fold_changes

    treatment = nb_counts(treatment_means, n_samples)
    return control, treatment, base_means


def cpm_normalize(counts):
    """Counts per million normalization."""
    lib_sizes = counts.sum(axis=0)
    return counts / lib_sizes * 1e6


def run_de_test(control_norm, treatment_norm):
    """Welch t-test on log-transformed CPM values."""
    eps = 0.5
    log_ctrl = np.log2(control_norm + eps)
    log_trt = np.log2(treatment_norm + eps)

    log2fc = log_trt.mean(axis=1) - log_ctrl.mean(axis=1)

    pvalues = np.array([
        stats.ttest_ind(log_trt[i], log_ctrl[i], equal_var=False).pvalue
        for i in range(log_ctrl.shape[0])
    ])

    _, padj, _, _ = multipletests(pvalues, method="fdr_bh")
    return log2fc, pvalues, padj


def plot_volcano(results_df, output_path):
    """Publication-ready volcano plot with significance thresholds."""
    fig, ax = plt.subplots(figsize=(8, 6))

    neg_log10_padj = -np.log10(results_df["padj"].clip(lower=1e-30))

    fc_thresh = 1.0
    pval_thresh = 0.05

    up = (results_df["log2FC"] >= fc_thresh) & (results_df["padj"] < pval_thresh)
    down = (results_df["log2FC"] <= -fc_thresh) & (results_df["padj"] < pval_thresh)
    ns = ~(up | down)

    ax.scatter(results_df.loc[ns, "log2FC"], neg_log10_padj[ns],
               color="#aaaaaa", alpha=0.5, s=12, label="Not significant")
    ax.scatter(results_df.loc[up, "log2FC"], neg_log10_padj[up],
               color="#d62728", alpha=0.8, s=18, label=f"Up-regulated (n={up.sum()})")
    ax.scatter(results_df.loc[down, "log2FC"], neg_log10_padj[down],
               color="#1f77b4", alpha=0.8, s=18, label=f"Down-regulated (n={down.sum()})")

    ax.axhline(-np.log10(pval_thresh), color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(fc_thresh, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(-fc_thresh, color="black", linestyle="--", linewidth=0.8, alpha=0.6)

    top_genes = results_df[up | down].nsmallest(8, "padj")
    for _, row in top_genes.iterrows():
        ax.annotate(row["gene"], xy=(row["log2FC"], -np.log10(row["padj"])),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points")

    ax.set_xlabel("Log₂ Fold Change (Treatment / Control)", fontsize=12)
    ax.set_ylabel("−log₁₀ Adjusted P-value (BH FDR)", fontsize=12)
    ax.set_title("Differential Expression: Volcano Plot", fontsize=14, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    de_indices = rng.choice(N_GENES, size=N_DE_GENES, replace=False)
    magnitudes = rng.uniform(1.5, 4.0, size=N_DE_GENES)
    directions = rng.choice([-1, 1], size=N_DE_GENES)
    # up: multiply, down: divide — keeps means strictly positive
    fold_changes = np.where(directions == 1, magnitudes, 1.0 / magnitudes)

    print("Simulating RNA-seq count data...")
    control, treatment, _ = simulate_counts(N_GENES, N_SAMPLES_PER_GROUP, de_indices, fold_changes)

    print("Normalizing (CPM)...")
    ctrl_norm = cpm_normalize(control)
    trt_norm = cpm_normalize(treatment)

    print("Running differential expression tests...")
    log2fc, pvalues, padj = run_de_test(ctrl_norm, trt_norm)

    gene_names = [f"GENE_{i:04d}" for i in range(N_GENES)]
    results = pd.DataFrame({
        "gene": gene_names,
        "log2FC": log2fc,
        "pvalue": pvalues,
        "padj": padj,
        "mean_ctrl_cpm": ctrl_norm.mean(axis=1),
        "mean_trt_cpm": trt_norm.mean(axis=1),
    }).sort_values("padj")

    sig = (results["padj"] < 0.05) & (results["log2FC"].abs() >= 1.0)
    print(f"\nResults: {sig.sum()} significant DE genes (FDR < 0.05, |log2FC| >= 1)")
    print(results.head(10).to_string(index=False))

    csv_path = os.path.join(RESULTS_DIR, "de_results.csv")
    results.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    print("\nGenerating volcano plot...")
    plot_volcano(results, os.path.join(RESULTS_DIR, "volcano_plot.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
