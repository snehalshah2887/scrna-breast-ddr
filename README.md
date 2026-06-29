# scrna-breast-ddr

Single-cell resolution of DNA-damage-response (DDR) gene programs in breast cancer.

This project applies single-cell RNA-seq to determine which cells in the breast tumour
microenvironment express an active DDR program, extending prior pan-cancer bulk RNA-seq work on
DDR-associated lncRNAs (GlncDDR, *Bioinformatics Advances* 2026) to single-cell resolution.

---

## Status

| Stage | Status |
|-------|--------|
| Standard scRNA-seq pipeline (established on PBMC3k) | Complete (`notebooks/00_explore.ipynb`) |
| Breast-cancer dataset and DDR analysis | Complete (`notebooks/01_breast_ddr.ipynb`) |
| Primary processing (reads to count matrix, Cell Ranger) | In progress (`notebooks/02_primary_processing.ipynb`) |

The PBMC3k notebook establishes and validates the downstream workflow (quality control,
normalisation, feature selection, dimensionality reduction, clustering, and marker-based
annotation), which is then applied to the breast-cancer DDR analysis in `01_breast_ddr.ipynb`.

The breast-cancer analysis in `01_breast_ddr.ipynb` is a complete single-cell analysis, run from raw
counts through to the biological result. The CELLxGENE atlas it uses was supplied already aligned and
quantified by the original authors, so the one upstream step not performed there is primary
processing: the conversion of sequencing reads into a count matrix. That step is documented separately
in `02_primary_processing.ipynb` (in progress), which uses Cell Ranger on 10x Genomics FASTQ files to
produce a count matrix that then enters the same downstream pipeline.

---

## Scientific context

DDR is the target of an entire drug class, the PARP inhibitors, and breast cancers with
homologous-recombination deficiency (HRD; for example BRCA1/2-mutant tumours) depend on it.
Resolving which malignant cells express an active DDR program is therefore a therapeutically
relevant, stratification-oriented question that bulk RNA-seq, which averages over all cell types,
cannot address.

---

## Analysis (`01_breast_ddr.ipynb`)

The breast-cancer notebook runs the full workflow on the Wu et al. (2021) atlas:

1. Load and inspect the dataset, and rebuild the working object from raw counts.
2. Quality control (light confirmatory thresholds, since the atlas is already pre-filtered).
3. Normalisation (counts per 10,000 and log1p).
4. Feature selection (2,000 batch-aware highly variable genes).
5. Scaling, PCA (50 components), and Harmony integration across the 26 patients.
6. Leiden clustering and UMAP.
7. Independent annotation by a three-method consensus (canonical markers, PanglaoDB
   over-representation analysis, and CellTypist label transfer), validated against the published
   labels.
8. DDR analysis (Section 9): a DDR-activity score from 491 DDR genes, its distribution across cell
   types and molecular subtypes, the detectability and DDR-association of the predicted lncRNAs, and
   their cell-type localisation.

Main result: the GlncDDR-predicted lncRNAs are expressed predominantly in malignant epithelial
cells, the compartment with the highest DDR activity, which yields a prioritised candidate list for
follow-up. A direct, cell-intrinsic DDR-regulatory role is not established by these data.

---

## Pipeline reference (`00_explore.ipynb`)

The PBMC3k notebook establishes and validates the standard Scanpy workflow (quality control,
normalisation, feature selection, dimensionality reduction, clustering, and marker-based
annotation) that is applied to the breast-cancer analysis.

---

## Repository structure

```
scrna-breast-ddr/
├── notebooks/
│   ├── 00_explore.ipynb              standard scRNA-seq pipeline (PBMC3k)
│   ├── 01_breast_ddr.ipynb           breast-cancer DDR analysis (CELLxGENE atlas)
│   └── 02_primary_processing.ipynb   primary processing: sequencing reads to count matrix (in progress)
├── data/
│   ├── glncddr_ddr_proteins.csv      491 DDR protein-coding genes (tracked)
│   ├── glncddr_lncrna_clean.csv      6,841 predicted DDR-lncRNAs, biotype-cleaned (tracked)
│   ├── glncddr_lncrna_DDRcorr.csv    per-lncRNA DDR correlation output (tracked)
│   └── (large inputs are gitignored; see Data)
├── src/
│   └── build_gene_lists.py           regenerate the GlncDDR-derived gene lists
├── figures/                          key figures exported from 01_breast_ddr.ipynb
├── environment.yml                   conda environment (pinned)
├── LICENSE                           MIT License
└── README.md
```

---

## Data

Two categories of data are used.

Tracked in the repository (small derived gene lists, inputs to Section 9):

- `data/glncddr_ddr_proteins.csv`: 491 DDR protein-coding genes (Pearl 2015, Knijnenburg 2018, Weir
  2022), extracted from the GlncDDR training positives.
- `data/glncddr_lncrna_clean.csv`: 6,841 GlncDDR-predicted lncRNAs (positive across all three models
  at probability at least 0.6), filtered to lncRNA biotype using Ensembl annotation.
- `data/glncddr_lncrna_DDRcorr.csv`: per-lncRNA correlation with the DDR score, generated by Section
  9.5.

The two input lists (`glncddr_ddr_proteins.csv` and `glncddr_lncrna_clean.csv`) can be regenerated from
the published GlncDDR model outputs with `python src/build_gene_lists.py`.

Not tracked (large files; download or generate locally):

- `data/raw/breast_atlas.h5ad`: Wu et al. (2021) breast atlas (GSE176078, approximately 100,000
  cells). Download from CZ CELLxGENE:

```bash
mkdir -p data/raw
curl -L -o data/raw/breast_atlas.h5ad \
  "https://datasets.cellxgene.cziscience.com/22a27631-aecf-463b-86c6-a8334a2f2cf2.h5ad"
```

- `data/human_cell_marker.txt`: CellMarker 2.0 human marker table, used for marker provenance in
  Section 8. Download from the CellMarker 2.0 database.

Source atlas: [CZ CELLxGENE, A single-cell and spatially resolved atlas of human breast cancers](https://cellxgene.cziscience.com/collections/dea97145-f712-431c-a223-6b5f565f362a).
GlncDDR predictions: [github.com/BioDataLearning/GlncDDR](https://github.com/BioDataLearning/GlncDDR).

---

## Setup and run

```bash
conda env create -f environment.yml
conda activate scrna
python -m ipykernel install --user --name scrna --display-name "Python (scrna)"
jupyter lab
```

Open a notebook and run it top to bottom (Kernel, then Restart and Run All).

---

## License

Released under the MIT License. See `LICENSE`.

---

## Author

Snehal Shah, PhD, Computational Biologist
- LinkedIn: https://linkedin.com/in/connectsnehalshah
- Related publication (GlncDDR): https://doi.org/10.1093/bioadv/vbag119
