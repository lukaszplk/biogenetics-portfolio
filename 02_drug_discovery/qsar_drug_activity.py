"""
QSAR Drug Activity Prediction
================================
ML pipeline for predicting drug bioactivity from molecular fingerprint features.
Demonstrates core cheminformatics / computational drug discovery workflow.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, classification_report
import matplotlib.pyplot as plt
import os

SEED = 42
N_MOLECULES = 1200
N_FEATURES = 256      # Morgan fingerprint bit length (ECFP4 analogue)
ACTIVE_RATE = 0.20    # 20% active compounds — realistic class imbalance
RESULTS_DIR = "results"

rng = np.random.default_rng(SEED)


def simulate_fingerprints(n_mols, n_bits, active_rate):
    """
    Simulate binary molecular fingerprint data with a hidden pharmacophore.
    A subset of bits is correlated with activity to mimic a real SAR signal.
    """
    n_active = int(n_mols * active_rate)
    n_inactive = n_mols - n_active

    pharmacophore_bits = rng.choice(n_bits, size=30, replace=False)

    active_fp = rng.binomial(1, 0.15, size=(n_active, n_bits)).astype(float)
    active_fp[:, pharmacophore_bits] = rng.binomial(1, 0.75, size=(n_active, len(pharmacophore_bits)))

    inactive_fp = rng.binomial(1, 0.12, size=(n_inactive, n_bits)).astype(float)

    X = np.vstack([active_fp, inactive_fp])
    y = np.array([1] * n_active + [0] * n_inactive)

    shuffle = rng.permutation(n_mols)
    return X[shuffle], y[shuffle], pharmacophore_bits


def build_models():
    return {
        "Logistic Regression": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(C=0.1, class_weight="balanced",
                                       max_iter=500, random_state=SEED))
        ]),
        "Random Forest": RandomForestClassifier(
            n_estimators=300, max_depth=8, class_weight="balanced",
            random_state=SEED, n_jobs=-1
        ),
        "Gradient Boosting": GradientBoostingClassifier(
            n_estimators=200, learning_rate=0.05, max_depth=4,
            subsample=0.8, random_state=SEED
        ),
    }


def evaluate_models(models, X, y):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    results = {}
    proba_dict = {}

    for name, model in models.items():
        proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
        auc = roc_auc_score(y, proba)
        pr_auc = average_precision_score(y, proba)
        proba_dict[name] = proba
        results[name] = {"ROC-AUC": auc, "PR-AUC": pr_auc}
        print(f"  {name:25s}  ROC-AUC={auc:.3f}  PR-AUC={pr_auc:.3f}")

    return results, proba_dict


def plot_roc_curves(y, proba_dict, output_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["#1f77b4", "#d62728", "#2ca02c"]

    for (name, proba), color in zip(proba_dict.items(), colors):
        fpr, tpr, _ = roc_curve(y, proba)
        auc = roc_auc_score(y, proba)
        ax.plot(fpr, tpr, color=color, lw=2, label=f"{name} (AUC={auc:.3f})")

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="Random classifier")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("QSAR Model Comparison — ROC Curves (5-fold CV)", fontsize=13, fontweight="bold")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_feature_importance(model, X, pharmacophore_bits, output_path):
    """Show which fingerprint bits drive predictions (Random Forest)."""
    importances = model.feature_importances_
    top_idx = np.argsort(importances)[::-1][:30]

    colors = ["#d62728" if i in pharmacophore_bits else "#aaaaaa" for i in top_idx]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(range(30), importances[top_idx], color=colors, edgecolor="none")
    ax.set_xticks(range(30))
    ax.set_xticklabels([f"Bit {i}" for i in top_idx], rotation=90, fontsize=7)
    ax.set_ylabel("Feature Importance (Gini)", fontsize=11)
    ax.set_title("Top 30 Molecular Fingerprint Bits — Feature Importance", fontsize=13, fontweight="bold")

    from matplotlib.patches import Patch
    legend = [Patch(color="#d62728", label="True pharmacophore bit"),
              Patch(color="#aaaaaa", label="Non-pharmacophore bit")]
    ax.legend(handles=legend, frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("Simulating molecular fingerprint dataset...")
    X, y = simulate_fingerprints(N_MOLECULES, N_FEATURES, ACTIVE_RATE)[:2]
    _, _, pharmacophore_bits = simulate_fingerprints(N_MOLECULES, N_FEATURES, ACTIVE_RATE)
    print(f"  {N_MOLECULES} molecules | {y.sum()} actives ({y.mean()*100:.1f}%) | {N_FEATURES} bits")

    print("\nTraining and evaluating models (5-fold CV)...")
    models = build_models()
    results, proba_dict = evaluate_models(models, X, y)

    results_df = pd.DataFrame(results).T.reset_index().rename(columns={"index": "Model"})
    csv_path = os.path.join(RESULTS_DIR, "model_comparison.csv")
    results_df.to_csv(csv_path, index=False)
    print(f"\n  Saved: {csv_path}")

    print("\nPlotting ROC curves...")
    plot_roc_curves(y, proba_dict, os.path.join(RESULTS_DIR, "roc_curves.png"))

    print("\nFitting Random Forest for feature importance...")
    rf = models["Random Forest"]
    rf.fit(X, y)
    plot_feature_importance(rf, X, set(pharmacophore_bits),
                            os.path.join(RESULTS_DIR, "feature_importance.png"))
    print("\nDone.")


if __name__ == "__main__":
    main()
