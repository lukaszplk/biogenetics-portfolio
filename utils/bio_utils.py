"""
Shared utility functions used across portfolio projects.
"""

import numpy as np
import matplotlib.pyplot as plt


def set_publication_style():
    """Apply consistent matplotlib style for publication-quality figures."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "figure.dpi": 100,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def benjamini_hochberg(pvalues):
    """
    BH FDR correction.
    Returns adjusted p-values (q-values).
    """
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked_pvals = pvalues[order]
    adjusted = np.minimum(1.0, ranked_pvals * n / (np.arange(1, n + 1)))
    # enforce monotonicity
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result
