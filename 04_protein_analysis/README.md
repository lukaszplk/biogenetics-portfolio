# Protein Sequence Feature Extraction

Physicochemical analysis and clustering of human protein sequences using Biopython.

## What it does

- Parses FASTA files and strips non-standard residues
- Computes ProtParam descriptors: molecular weight, isoelectric point, GRAVY score, instability index, aromaticity, secondary structure fractions
- Builds amino acid composition profiles
- Clusters proteins by combined sequence features (KMeans + PCA)
- Generates four plots: descriptor heatmap, AA composition bar chart, PCA scatter, MW vs pI scatter

## Proteins analysed

| Protein | UniProt | Function |
|---|---|---|
| BRCA2 | NP_000050.2 | DNA repair, tumour suppressor |
| TP53 | NP_000537.3 | Cell cycle regulation, apoptosis |
| EGFR | NP_005219.2 | Cell signalling receptor |
| Insulin | NP_000198.1 | Glucose metabolism hormone |
| β-Actin | NP_001092.1 | Cytoskeletal structural protein |

## Usage

```bash
python protein_analysis.py
```

Outputs saved to `results/`:
- `descriptor_heatmap.png`
- `aa_composition.png`
- `pca_clusters.png`
- `mw_vs_pi.png`
- `descriptors.csv`
- `aa_composition.csv`

## Tools

`Biopython` · `pandas` · `scikit-learn` · `matplotlib` · `seaborn`
