# scrna-breast-ddr

**Single-cell resolution of DNA damage response (DDR) gene programs in breast cancer.**

This project applies single-cell RNA-seq to ask *which cells* in the breast tumor
microenvironment carry an active DDR program — extending prior **pan-cancer bulk**
RNA-seq work on DDR-associated lncRNAs (GlncDDR, *Bioinformatics Advances* 2026) to
**single-cell resolution**.

---

## Status

| Stage | Status |
|-------|--------|
| Standard scRNA-seq pipeline (established on PBMC3k) | ✅ Complete — `notebooks/00_explore.ipynb` |
| Breast-cancer dataset + DDR signature analysis | 🚧 In progress |

The PBMC3k notebook establishes and validates the end-to-end workflow (QC →
normalization → feature selection → dimensionality reduction → clustering →
marker-based annotation) that is then applied to the breast-cancer DDR analysis.

---

## Scientific context

DDR is the target of an entire drug class — PARP inhibitors — and breast cancers with
homologous recombination deficiency (HRD; e.g. BRCA1/2-mutant) depend on it. Resolving
*which* malignant cells carry an active DDR program is therefore a therapeutically
relevant, stratification-oriented question that bulk RNA-seq (which averages over all
cell types) cannot answer.

---

## Pipeline (PBMC3k walkthrough)

`notebooks/00_explore.ipynb` runs the standard Scanpy workflow end-to-end:

1. **QC** — per-cell metrics, mitochondrial fraction, distribution-based filtering
2. **Normalization** — counts-per-10k (CP10K) + log1p
3. **Feature selection** — 2,000 highly variable genes
4. **Dimensionality reduction** — scaling → PCA (30 PCs) → neighbor graph
5. **Clustering & visualization** — Leiden clustering, UMAP
6. **Annotation** — canonical markers (CST3 / NKG7 / MS4A1 → monocytes / NK / B cells)

---

## Repository structure

```
scrna-breast-ddr/
├── notebooks/
│   └── 00_explore.ipynb     # standard scRNA-seq pipeline (PBMC3k)
├── data/                    # datasets (gitignored)
├── src/                     # reusable functions
├── figures/                 # exported figures
├── environment.yml          # conda environment
└── README.md
```

---

## Setup & run

```bash
conda env create -f environment.yml
conda activate scrna
python -m ipykernel install --user --name scrna --display-name "Python (scrna)"
jupyter lab
```

Then open `notebooks/00_explore.ipynb` and run top to bottom
(Kernel → Restart & Run All).

---

## Author

**Snehal Shah, PhD** — Computational Biologist
- LinkedIn: https://linkedin.com/in/connectsnehalshah
- Related publication (GlncDDR): https://doi.org/10.1093/bioadv/vbag119
