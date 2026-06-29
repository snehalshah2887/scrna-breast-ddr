"""
build_gene_lists.py

Reproduce the two GlncDDR-derived gene lists used by notebooks/01_breast_ddr.ipynb:

  data/glncddr_ddr_proteins.csv   491 DDR protein-coding training positives
  data/glncddr_lncrna_clean.csv   predicted DDR-lncRNAs, biotype-cleaned to lncRNA only

Both are derived from the published GlncDDR model outputs
(https://github.com/BioDataLearning/GlncDDR). This script documents the exact
derivation so the committed CSVs are reproducible.

Usage:
    python src/build_gene_lists.py

Requirements: pandas, openpyxl (Excel), and network access (GlncDDR repo + Ensembl REST).
"""

import json
import os
import time
import urllib.request

import pandas as pd

GLNCDDR = "https://raw.githubusercontent.com/BioDataLearning/GlncDDR/main"
OUT = os.path.join(os.path.dirname(__file__), "..", "data")
PROBA_THRESHOLD = 0.6  # a transcript must be positive in all three models at this probability

# Ensembl biotypes that count as a long non-coding RNA
LNCRNA_BIOTYPES = {
    "lncRNA", "lincRNA", "antisense", "antisense_RNA", "processed_transcript",
    "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA",
    "bidirectional_promoter_lncRNA", "macro_lncRNA", "non_coding",
}


def _download(url, path):
    if not os.path.exists(path):
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, path)
    return path


def build_ddr_proteins(workdir):
    """491 DDR protein-coding positives = Gene_type == 1 in the protein training embeddings."""
    path = _download(f"{GLNCDDR}/embeddings/train_emb_len100.csv",
                     os.path.join(workdir, "train_emb_len100.csv"))
    te = pd.read_csv(path)
    genes = sorted(te.loc[te["Gene_type"] == 1, "Genes"].dropna().unique())
    df = pd.DataFrame({"Genes": genes})
    df.to_csv(os.path.join(OUT, "glncddr_ddr_proteins.csv"), index=False)
    print(f"DDR proteins: {len(df)} genes -> data/glncddr_ddr_proteins.csv")
    return df


def ensembl_biotypes(ensembl_ids, workdir):
    """Look up the Ensembl biotype for each gene via the Ensembl REST API (batched, cached)."""
    cache_path = os.path.join(workdir, "biotype_cache.json")
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    todo = [g for g in ensembl_ids if g not in cache]
    server = "https://rest.ensembl.org/lookup/id"
    for i in range(0, len(todo), 800):
        batch = todo[i:i + 800]
        req = urllib.request.Request(
            server, data=json.dumps({"ids": batch}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        r = None
        for _ in range(5):
            try:
                r = json.loads(urllib.request.urlopen(req, timeout=90).read()); break
            except Exception:
                time.sleep(3)
        if r:
            for k, v in r.items():
                cache[k] = (v or {}).get("biotype")
            for g in batch:
                cache.setdefault(g, None)
        json.dump(cache, open(cache_path, "w"))
        print(f"  biotypes {min(i + 800, len(todo))}/{len(todo)}")
    return cache


def build_lncrnas(workdir):
    """Predicted lncRNAs = positive in all three models at proba >= 0.6, filtered to lncRNA biotype."""
    path = _download(f"{GLNCDDR}/output_dir/lncrna_prediction.xlsx",
                     os.path.join(workdir, "lncrna_prediction.xlsx"))
    sheets = pd.read_excel(path, sheet_name=None)  # one sheet per model (LR, RF, SVM)

    # consensus: kept only if positive at >= threshold in every model
    sets = [set(df.loc[df["predict_proba"] >= PROBA_THRESHOLD, "Ensembl"]) for df in sheets.values()]
    consensus = set.intersection(*sets)
    base = next(iter(sheets.values()))[["Ensembl", "Genes"]]
    lst = base[base["Ensembl"].isin(consensus)].drop_duplicates("Ensembl").reset_index(drop=True)
    print(f"3-model consensus (proba >= {PROBA_THRESHOLD}): {len(lst)} transcripts")

    # biotype clean: keep only true lncRNAs (drops pseudogenes, miRNA, snRNA, etc.)
    bt = ensembl_biotypes(lst["Ensembl"].tolist(), workdir)
    lst["biotype"] = lst["Ensembl"].map(bt)
    lst = lst[lst["biotype"].isin(LNCRNA_BIOTYPES)].copy()

    # attach the per-model probabilities for provenance
    for name, df in sheets.items():
        lst[name] = lst["Ensembl"].map(dict(zip(df["Ensembl"], df["predict_proba"])))
    prob_cols = list(sheets.keys())
    lst["min_proba"] = lst[prob_cols].min(axis=1)
    lst["mean_proba"] = lst[prob_cols].mean(axis=1)

    lst.to_csv(os.path.join(OUT, "glncddr_lncrna_clean.csv"), index=False)
    print(f"lncRNAs (biotype-cleaned): {len(lst)} -> data/glncddr_lncrna_clean.csv")
    return lst


def main():
    workdir = os.path.join(os.path.dirname(__file__), "..", "data", "_build")
    os.makedirs(workdir, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    build_ddr_proteins(workdir)
    build_lncrnas(workdir)
    print("done.")


if __name__ == "__main__":
    main()
