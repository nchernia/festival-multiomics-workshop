#!/usr/bin/env python3
"""Bake bias-corrected TF footprints into the workshop .h5mu (instructor post-step).

For each TF motif, scans the top peaks accessible in a target cell type for motif instances,
builds a Tn5 6-mer insertion-bias model from the fragments, and aggregates
bias-corrected Tn5 insertions around the motif centers per cell type. The
resulting per-cell-type footprint profiles are stored in mdata.uns so the
student notebook only has to plot them (no fragments, no compute, can't break).

Requires the ATAC fragments (tabix-indexed) + an hg38 FASTA (bgzipped+faidx).
CTCF is included as a positive control (its footprint is textbook-clean).

Run:
    python preprocessing/refresh_footprints.py \
        --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu \
        --fragments workshop_data/raw/atac_fragments.tsv.gz \
        --fasta /path/to/hg38.fa.gz
"""
import argparse, time
import numpy as np, pandas as pd, scipy.sparse as sp
import muon as mu, pysam
from pyjaspar import jaspardb

TFS = {"EBF1": "B-lineage", "CTCF": "control"}   # discovered TF + control
TARGET_CT = "B cell"            # cell type whose accessible peaks we scan
W, PVAL, K, TOPN = 150, 5e-4, 6, 6000
B = {"A": 0, "C": 1, "G": 2, "T": 3}


def kmer_id(seq):
    v = 0
    for ch in seq:
        b = B.get(ch, -1)
        if b < 0:
            return -1
        v = v * 4 + b
    return v


def pwm_logodds(name):
    mot = jaspardb(release="JASPAR2024").fetch_motifs_by_name(name)[0]
    counts = np.array([list(mot.counts[b]) for b in "ACGT"], float)
    lod = np.log2((counts + 0.25) / (counts.sum(0) + 1.0) / 0.25)
    bg = np.array([lod[np.random.default_rng(s).integers(0, 4, lod.shape[1]),
                       np.arange(lod.shape[1])].sum() for s in range(30000)])
    return lod, float(np.quantile(bg, 1 - PVAL)), mot.matrix_id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5mu", default="workshop_data/pbmc_10k_multiome_workshop.h5mu")
    ap.add_argument("--fragments", default="workshop_data/raw/atac_fragments.tsv.gz")
    ap.add_argument("--fasta", required=True, help="hg38 FASTA (bgzipped + .fai)")
    args = ap.parse_args()

    fa = pysam.FastaFile(args.fasta)
    frags = pysam.TabixFile(args.fragments)
    m = mu.read(args.h5mu)
    atac = m.mod["atac"]
    ct = dict(zip(atac.obs_names, atac.obs["cell_type"].astype(str)))
    order = list(atac.obs["cell_type"].cat.categories)
    n_cells = {c: sum(1 for v in ct.values() if v == c) for c in order}

    # top peaks accessible in the target cell type
    tgt_idx = np.array([i for i, b in enumerate(atac.obs_names) if ct[b] == TARGET_CT])
    X = atac.X.tocsr() if sp.issparse(atac.X) else sp.csr_matrix(atac.X)
    peakmean = np.asarray(X[tgt_idx].mean(0)).ravel()
    top = np.argsort(-peakmean)[:TOPN]
    pc = atac.var_names[top].str.extract(r"(chr\w+)[:\-](\d+)[:\-](\d+)")
    peaks = list(zip(pc[0], pc[1].astype(int), pc[2].astype(int)))
    print(f"scanning {len(peaks)} {TARGET_CT}-accessible peaks")

    # Tn5 6-mer bias model
    print("building Tn5 6-mer bias model ...")
    obs = np.zeros(4 ** K); n = 0; t0 = time.time()
    for c in ["chr1", "chr3", "chr11"]:
        for row in frags.fetch(c, 1, 60_000_000):
            f = row.split("\t"); fs, fe = int(f[1]), int(f[2])
            for ins in (fs + 4, fe - 5):
                kid = kmer_id(fa.fetch(c, ins - K // 2, ins - K // 2 + K).upper())
                if kid >= 0:
                    obs[kid] += 1
            n += 1
            if n >= 400_000:
                break
        if n >= 400_000:
            break
    exp = np.zeros(4 ** K)
    for c, s, e in peaks:
        seq = fa.fetch(c, s, e).upper()
        for i in range(len(seq) - K + 1):
            kid = kmer_id(seq[i:i + K])
            if kid >= 0:
                exp[kid] += 1
    bias = (obs / obs.sum()) / ((exp / exp.sum()) + 1e-9)
    bias[~np.isfinite(bias)] = 1.0
    print(f"  bias model from {n} fragments in {time.time()-t0:.0f}s")

    def scan(lod, thr):
        L = lod.shape[1]; lod_rc = lod[::-1, ::-1]; inst = []
        for c, s, e in peaks:
            seq = np.array([B.get(ch, -1) for ch in fa.fetch(c, s, e).upper()])
            for i in range(len(seq) - L + 1):
                w = seq[i:i + L]
                if (w < 0).any():
                    continue
                if lod[w, np.arange(L)].sum() >= thr:
                    inst.append((c, s + i + L // 2, +1))
                elif lod_rc[w, np.arange(L)].sum() >= thr:
                    inst.append((c, s + i + L // 2, -1))
        return inst

    foot = {
        "tf": np.array(list(TFS), dtype=object),
        "role": np.array([TFS[t] for t in TFS], dtype=object),
        "offsets": np.arange(-W, W + 1, dtype="int64"),
        "cell_types": np.array(order, dtype=object),
        "target_celltype": np.array([TARGET_CT], dtype=object),
        "n_instances": [],
        "motif_id": [],
    }
    for name in TFS:
        lod, thr, mid = pwm_logodds(name)
        inst = scan(lod, thr)
        print(f"{name} ({mid}): {len(inst)} motif instances")
        prof = {c: np.zeros(2 * W + 1) for c in order}
        expect = np.zeros(2 * W + 1)
        for c, center, strand in inst:
            wseq = fa.fetch(c, center - W - K, center + W + K).upper()
            ev = np.array([bias[kmer_id(wseq[j:j + K])] if kmer_id(wseq[j:j + K]) >= 0 else 1.0
                           for j in range(2 * W + 1)])
            if strand < 0:
                ev = ev[::-1]
            expect += ev
            for row in frags.fetch(c, max(0, center - W - 5), center + W + 5):
                f = row.split("\t"); cell = ct.get(f[3])
                if cell is None:
                    continue
                fs, fe = int(f[1]), int(f[2])
                for ins in (fs + 4, fe - 5):
                    d = ins - center
                    if strand < 0:
                        d = -d
                    if -W <= d <= W:
                        prof[cell][d + W] += 1
        expect = expect / len(inst)
        expect_norm = expect / expect.mean()
        # bias-corrected insertions per 1000 cells per motif, per cell type
        P = np.vstack([
            (prof[c] / max(n_cells[c], 1) / len(inst) * 1e3) / expect_norm
            for c in order
        ]).astype("float32")
        foot[f"profile_{name}"] = P
        foot["n_instances"].append(len(inst))
        foot["motif_id"].append(mid)
    foot["n_instances"] = np.array(foot["n_instances"], dtype="int64")
    foot["motif_id"] = np.array(foot["motif_id"], dtype=object)

    m.uns["footprint"] = foot
    m.write(args.h5mu)
    print(f"baked footprint for {list(TFS)} into {args.h5mu}")


if __name__ == "__main__":
    main()
