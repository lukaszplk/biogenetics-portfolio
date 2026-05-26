"""
Protein sequence feature extraction and analysis.

Demonstrates:
- Parsing FASTA files with Biopython
- Computing physicochemical descriptors (MW, pI, GRAVY, instability index)
- Amino acid composition and dipeptide frequency
- Clustering proteins by sequence features
- Visualising descriptor distributions
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
FASTA_FILE = os.path.join(DATA_DIR, "proteins.fasta")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


# ── Parsing ────────────────────────────────────────────────────────────────────

def load_proteins(fasta_path: str) -> dict[str, str]:
    """Return {id: sequence} from a FASTA file, stripping non-standard residues."""
    proteins = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        seq_clean = "".join(aa for aa in seq if aa in AMINO_ACIDS)
        if seq_clean:
            proteins[record.id] = seq_clean
    return proteins


# ── Physicochemical descriptors ────────────────────────────────────────────────

def compute_descriptors(seq: str) -> dict:
    """Compute ProtParam descriptors for a single protein sequence."""
    analysis = ProteinAnalysis(seq)
    return {
        "length": len(seq),
        "molecular_weight": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(),
        "gravy": analysis.gravy(),
        "instability_index": analysis.instability_index(),
        "aromaticity": analysis.aromaticity(),
        "helix_fraction": analysis.secondary_structure_fraction()[0],
        "sheet_fraction": analysis.secondary_structure_fraction()[1],
        "coil_fraction": analysis.secondary_structure_fraction()[2],
    }


def build_descriptor_table(proteins: dict[str, str]) -> pd.DataFrame:
    rows = []
    for prot_id, seq in proteins.items():
        desc = compute_descriptors(seq)
        desc["protein"] = prot_id.split("|")[0]
        rows.append(desc)
    return pd.DataFrame(rows).set_index("protein")


# ── Amino acid composition ─────────────────────────────────────────────────────

def aa_composition(seq: str) -> dict[str, float]:
    """Return relative frequency of each standard amino acid."""
    counts = Counter(seq)
    total = len(seq)
    return {aa: counts.get(aa, 0) / total for aa in AMINO_ACIDS}


def build_composition_matrix(proteins: dict[str, str]) -> pd.DataFrame:
    rows = {
        prot_id.split("|")[0]: aa_composition(seq)
        for prot_id, seq in proteins.items()
    }
    return pd.DataFrame(rows).T


# ── Clustering ─────────────────────────────────────────────────────────────────

def cluster_proteins(feature_matrix: pd.DataFrame, n_clusters: int = 2) -> pd.Series:
    scaler = StandardScaler()
    X = scaler.fit_transform(feature_matrix)
    labels = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(X)
    return pd.Series(labels, index=feature_matrix.index, name="cluster")


# ── Visualisation ──────────────────────────────────────────────────────────────

def plot_descriptor_heatmap(desc_df: pd.DataFrame, out_path: str) -> None:
    scaled = pd.DataFrame(
        StandardScaler().fit_transform(desc_df),
        index=desc_df.index,
        columns=desc_df.columns,
    )
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.heatmap(
        scaled,
        annot=True,
        fmt=".2f",
        cmap="RdYlBu_r",
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title("Physicochemical Descriptors (z-scored)", fontsize=13, pad=12)
    ax.set_xlabel("Descriptor")
    ax.set_ylabel("Protein")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def plot_aa_composition(comp_df: pd.DataFrame, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(14, 5))
    comp_df.T.plot(kind="bar", ax=ax, width=0.75)
    ax.set_title("Amino Acid Composition by Protein", fontsize=13, pad=12)
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("Relative Frequency")
    ax.legend(title="Protein", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def plot_pca(feature_matrix: pd.DataFrame, clusters: pd.Series, out_path: str) -> None:
    scaler = StandardScaler()
    X = scaler.fit_transform(feature_matrix)
    pcs = PCA(n_components=2).fit_transform(X)
    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=feature_matrix.index)
    pca_df["cluster"] = clusters

    fig, ax = plt.subplots(figsize=(7, 5))
    for cluster_id, group in pca_df.groupby("cluster"):
        ax.scatter(group["PC1"], group["PC2"], label=f"Cluster {cluster_id}", s=80)
    for name, row in pca_df.iterrows():
        ax.annotate(name, (row["PC1"], row["PC2"]), fontsize=8,
                    xytext=(5, 5), textcoords="offset points")
    ax.set_title("PCA of Protein Sequence Features", fontsize=13)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def plot_scatter_mw_pi(desc_df: pd.DataFrame, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    scatter = ax.scatter(
        desc_df["molecular_weight"] / 1000,
        desc_df["isoelectric_point"],
        c=desc_df["gravy"],
        cmap="coolwarm",
        s=100,
        edgecolors="k",
        linewidths=0.5,
    )
    for name, row in desc_df.iterrows():
        ax.annotate(name, (row["molecular_weight"] / 1000, row["isoelectric_point"]),
                    fontsize=8, xytext=(5, 5), textcoords="offset points")
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("GRAVY score")
    ax.set_title("Molecular Weight vs Isoelectric Point", fontsize=13)
    ax.set_xlabel("Molecular Weight (kDa)")
    ax.set_ylabel("Isoelectric Point (pI)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("Loading proteins...")
    proteins = load_proteins(FASTA_FILE)
    print(f"  {len(proteins)} proteins loaded: {list(proteins.keys())}")

    print("\nComputing physicochemical descriptors...")
    desc_df = build_descriptor_table(proteins)
    print(desc_df[["length", "molecular_weight", "isoelectric_point",
                    "gravy", "instability_index"]].to_string())

    print("\nComputing amino acid composition...")
    comp_df = build_composition_matrix(proteins)

    print("\nClustering proteins by combined features...")
    combined = pd.concat([desc_df, comp_df], axis=1)
    clusters = cluster_proteins(combined, n_clusters=min(2, len(proteins)))
    print(clusters.to_string())

    print("\nGenerating plots...")
    plot_descriptor_heatmap(
        desc_df, os.path.join(OUTPUT_DIR, "descriptor_heatmap.png"))
    plot_aa_composition(
        comp_df, os.path.join(OUTPUT_DIR, "aa_composition.png"))
    plot_pca(
        combined, clusters, os.path.join(OUTPUT_DIR, "pca_clusters.png"))
    plot_scatter_mw_pi(
        desc_df, os.path.join(OUTPUT_DIR, "mw_vs_pi.png"))

    desc_df.to_csv(os.path.join(OUTPUT_DIR, "descriptors.csv"))
    comp_df.to_csv(os.path.join(OUTPUT_DIR, "aa_composition.csv"))
    print("\nAll outputs saved to results/")


if __name__ == "__main__":
    main()
