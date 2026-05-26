# Population Genetics — GWAS Visualization

Tools for visualizing Genome-Wide Association Study (GWAS) summary statistics.
Produces the two standard publication figures: Manhattan plot and QQ plot.

## What It Does

1. Generates synthetic GWAS summary statistics (24 chromosomes, 50k SNPs)
2. Injects simulated trait-associated loci on several chromosomes
3. Plots a Manhattan plot with genome-wide and suggestive significance thresholds
4. Plots a QQ plot with genomic inflation factor λ (lambda)

## Output

- `results/manhattan_plot.png` — chromosome-level association landscape
- `results/qq_plot.png` — p-value distribution quality control plot

## Run

```bash
python gwas_visualization.py
```

## Key Concepts

- GWAS summary statistics (CHR, BP, SNP, P)
- Genome-wide significance threshold (5×10⁻⁸)
- Suggestive significance threshold (1×10⁻⁵)
- Genomic inflation factor λ for QC
- LD clumping concept (peak SNP labeling)
