#!/usr/bin/env python3
"""Refresh the AlphaGenome cache in the workshop .h5mu for the MS4A1 enhancer.

Computes the featured MS4A1 peak the same way the notebook does (paired
correlation over a ±250 kb TSS window), then caches AlphaGenome CHIP-TF
predictions for it. Dedups one row per TF (strongest biosample by max
signal in the enhancer ±2 kb), matching the live code path in the notebook.

Run:
    ALPHA_GENOME_API_KEY=<key> python preprocessing/refresh_alphagenome.py \
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
# TF panel matches the notebook's v2 live code path.
TARGET_TFS = ["EBF1","PAX5","SPIB","TCF3","POU2F2","SPI1","CEBPA","CEBPB","IRF8",
              "RUNX3","ETS1","TBX21","EOMES","TCF7","LEF1","GATA3","BCL11A","CTCF"]
# Primary PBMC immune-cell ontology terms. Preferred over cell lines (GM12878,
# K562) for PBMC interpretation — predictions come back conditioned on the
# right lineage. Add "EFO:0002784" (GM12878) as a fallback for high-coverage
# B-cell ChIP-TF if a TF lacks primary-cell training data.
BIOSAMPLES = [
    "CL:0000236",   # B cell
    "CL:0000623",   # natural killer cell
    "CL:0000624",   # CD4-positive, alpha-beta T cell
    "CL:0000625",   # CD8-positive, alpha-beta T cell
    "CL:0000576",   # monocyte
    "CL:0001054",   # CD14-positive monocyte
]


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
    # Some CL/EFO terms are not in AlphaGenome's vocabulary -- drop and retry on error.
    terms = list(BIOSAMPLES)
    while True:
        try:
            out = model.predict_interval(
                interval=context_iv,
                requested_outputs=[dna_client.OutputType.CHIP_TF],
                ontology_terms=terms,
            )
            break
        except Exception as exc:
            mt = re.search(r'Unsupported ontology: "([^"]+)"', str(exc))
            if not mt or not terms:
                raise
            bad = mt.group(1)
            print(f"  AlphaGenome does not support {bad}; dropping and retrying.")
            terms = [t for t in terms if t != bad]
    print(f"  biosamples used: {terms}")
    tf_pred = out.chip_tf
    tf_md = tf_pred.metadata

    # Dedup per TF: for each TF in the panel that has any returned tracks, keep
    # the (TF, biosample) row whose max signal inside the enhancer ±2 kb is highest.
    # Matches the notebook's load_ag_live dedup so cache and live agree.
    tf_names = tf_md["transcription_factor"].astype(str).values
    pos_ctx = np.linspace(int(context_iv.start), int(context_iv.end), tf_pred.values.shape[0])
    win_mask = (pos_ctx >= start - 2000) & (pos_ctx <= end + 2000)
    sig_per_row = tf_pred.values[win_mask].max(axis=0)
    best = {}                                       # TF -> (row_index, signal)
    for i, t in enumerate(tf_names):
        if t not in TARGET_TFS:
            continue
        if t not in best or sig_per_row[i] > best[t][1]:
            best[t] = (int(i), float(sig_per_row[i]))
    # Sort TFs by descending signal for stable, readable order.
    items = sorted(best.items(), key=lambda kv: -kv[1][1])
    keep_idx = [v[1][0] for v in items]
    kept     = [k for k, _ in items]
    print(f"Found {len(kept)} target TF tracks (dedup, by strongest biosample): {kept}")
    print(f"  (out of {len(tf_md)} total CHIP-TF tracks the model returned)")
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
