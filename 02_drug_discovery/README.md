# ML Drug Activity Prediction (QSAR)

Quantitative Structure-Activity Relationship (QSAR) modeling pipeline for predicting
drug bioactivity from molecular descriptors. A core workflow in computational drug discovery.

## What It Does

1. Generates a synthetic dataset of molecules with Morgan-fingerprint-like binary features
2. Trains and compares multiple ML classifiers (Logistic Regression, Random Forest, XGBoost-style Gradient Boosting)
3. Evaluates models with ROC-AUC, PR-AUC, and classification report
4. Plots ROC curves and feature importance

## Output

- `results/roc_curves.png` — ROC comparison across models
- `results/feature_importance.png` — top predictive molecular features
- `results/model_comparison.csv` — AUC scores for all models

## Run

```bash
python qsar_drug_activity.py
```

## Key Concepts

- Binary molecular fingerprints (Morgan/ECFP analogues)
- Class imbalance handling (SMOTE / class weighting)
- Model comparison and ROC-AUC evaluation
- Feature importance for interpretability in lead optimization
