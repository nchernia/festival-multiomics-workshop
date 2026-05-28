#!/usr/bin/env python3
"""Refresh the AlphaGenome cache in the workshop .h5mu for the CD8A enhancer.

Computes the featured CD8A peak the same way the notebook does (paired
correlation over a ±250 kb TSS window), then caches AlphaGenome CHIP-TF
predictions (RUNX/ETS T-cell factors) for it. No GLUE retrain.

Run:
    ALPHA_GENOME_API_KEY=<key> python preprocessing/refresh_alphagenome_cd8a.py \
        --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu
"""
import argparse, os, re, sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats
import muon as mu

GENE = "MS4A1"
WINDOW = 250_000
TARGET_TFS = ["EBF1","PAX5","SPIB","TCF3","POU2F2","SPI1","IRF8"]
# CHIP-TF biosamples to request: K562 + GM12878 (lymphoblastoid) have the
# broadest ENCODE TF ChIP coverage; add CD8/CD4 T-cell terms in case present.
BIOSAMPLES = ["EFO:0002784", "EFO:0002067", "CL:0000236"]  # GM12878, K562, B cell


def dense(ad, name):
    col = ad[:, name].X
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()


def featured_peak(rna, atac, gene, tss_lookup, window=WINDOW):
    chrom, tss = tss_lookup[gene]
    pc = atac.var_names.str.extract(r"(?P<chrom>chr\w+)[:\-](?P<s>\d+)[:\-](?P<e>\d+)")
    mids = ((pc["s"].astype(int) + pc["e"].astype(int)) // 2).values
    in_win = (pc["chrom"].values == chrom) & (mids >= tss - window) & (mids <= tss + window)
    names = atac.var_names[in_win]
    g = dense(rna, gene)
    best, best_rho = None, -2
    for nm in names:
        a = dense(atac, nm)
        if a.std() == 0:
            continue
        # require DISTAL (>10 kb from TSS) so we feature an enhancer, not the promoter
        mid = (int(nm.split(":")[1].split("-")[0]) + int(nm.split("-")[1])) // 2
        if abs(mid - tss) < 10_000:
            continue
        r = stats.spearmanr(a, g)[0]
        if r > best_rho:
            best_rho, best = r, nm
    return best, best_rho


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5mu", default="workshop_data/pbmc_10k_multiome_workshop.h5mu")
    args = ap.parse_args()

    key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not key:
        sys.exit("Set ALPHA_GENOME_API_KEY (run: ALPHA_GENOME_API_KEY=... python ...)")

    from alphagenome.data import genome
    from alphagenome.models import dna_client

    mdata = mu.read(args.h5mu)
    rna, atac = mdata.mod["rna"], mdata.mod["atac"]
    _tss = mdata.uns["gene_tss"]
    TSS = {g: (c, int(t)) for g, c, t in zip(_tss["gene"], _tss["chrom"], _tss["tss"])}

    # Correlation needs log-normalized accessibility. rna.X is already
    # log-normalized in the .h5mu; atac.X is raw counts. Normalize a COPY of
    # atac so we never mutate (and re-save) the stored matrices.
    atac_n = atac.copy()
    sc.pp.normalize_total(atac_n, target_sum=1e4); sc.pp.log1p(atac_n)

    peak, rho = featured_peak(rna, atac_n, GENE, TSS)
    print(f"Featured {GENE} peak: {peak}  (cell-level rho={rho:.3f})")
    m = re.match(r"(chr\w+)[:\-](\d+)[:\-](\d+)", peak)
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

    model = dna_client.create(key)
    peak_iv = genome.Interval(chromosome=chrom, start=start, end=end)
    context_iv = peak_iv.resize(dna_client.SEQUENCE_LENGTH_1MB)

    print(f"predict_interval(CHIP_TF) at {context_iv} ...")
    out = model.predict_interval(
        interval=context_iv,
        requested_outputs=[dna_client.OutputType.CHIP_TF],
        ontology_terms=BIOSAMPLES,
    )
    tf_pred = out.chip_tf
    tf_md = tf_pred.metadata

    # Which of our target TFs actually have tracks?
    keep_idx, kept = [], []
    str_cols = tf_md.select_dtypes(include="object")
    for tf in TARGET_TFS:
        mask = str_cols.apply(lambda c: c.str.upper().str.contains(tf, na=False)).any(axis=1)
        hits = np.where(mask.values)[0]
        if len(hits):
            keep_idx.append(int(hits[0])); kept.append(tf)
    keep_idx = sorted(set(keep_idx))
    print(f"Found {len(kept)} target TF tracks: {kept}  (of {len(tf_md)} total CHIP-TF tracks)")
    if not kept:
        print("None of the target TFs found. First 30 available TF track names:")
        print(list(tf_md.iloc[:30, 0]))
        sys.exit("Adjust TARGET_TFS / BIOSAMPLES and re-run.")

    tf_values = tf_pred.values[:, keep_idx].astype(np.float32)
    tf_md_keep = tf_md.iloc[keep_idx].reset_index(drop=True)
    tf_md_dict = {c: np.asarray(tf_md_keep[c].astype(str).values, dtype=object)
                  for c in tf_md_keep.columns}

    cache = {
        "featured_peak": peak,
        "featured_gene": GENE,
        "context_chrom": chrom,
        "context_start": int(context_iv.start),
        "context_end":   int(context_iv.end),
        "peak_start":    int(start),
        "peak_end":      int(end),
        "chip_tf_values":   tf_values,
        "chip_tf_metadata": tf_md_dict,
    }
    mdata.uns["alphagenome_cache"] = cache
    mdata.write(args.h5mu)
    print(f"Wrote refreshed cache to {args.h5mu}")
    print(f"  featured_peak = {peak}")
    print(f"  TFs cached    = {kept}")


if __name__ == "__main__":
    main()
