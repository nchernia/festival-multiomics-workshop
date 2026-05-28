#!/usr/bin/env python3
"""Make per-cell-type pseudobulk ATAC bigWigs for the igv-notebook browser view.
Tn5-insertion coverage, normalized per 1,000 cells. Two modes:

  --region chr11:60430000-60520000   one locus, base-pair resolution (tiny files)
  --genome-wide                       all main chromosomes, binned (default 25 bp)

Genome-wide is processed chrom-by-chrom (low memory) so IGV is fully navigable.

Run:
    python preprocessing/make_pseudobulk_bigwigs.py --genome-wide \
        --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu \
        --fragments workshop_data/raw/atac_fragments.tsv.gz --outdir workshop_data/tracks
"""
import argparse, os
import numpy as np, muon as mu, pysam, pyBigWig

# main hg38 chromosome sizes (chr1-22, X, Y)
HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def _meta(args):
    m = mu.read(args.h5mu)
    atac = m.mod["atac"]
    ct = dict(zip(atac.obs_names, atac.obs["cell_type"].astype(str)))
    order = list(atac.obs["cell_type"].cat.categories)
    n_cells = {c: sum(1 for v in ct.values() if v == c) for c in order}
    return ct, order, n_cells


def run_locus(args):
    os.makedirs(args.outdir, exist_ok=True)
    chrom, se = args.region.split(":"); start, end = map(int, se.split("-"))
    L = end - start
    frags = pysam.TabixFile(args.fragments)
    ct, order, n_cells = _meta(args)
    cov = {c: np.zeros(L, dtype=np.float32) for c in order}
    for row in frags.fetch(chrom, start, end):
        f = row.split("\t"); cell = ct.get(f[3])
        if cell is None:
            continue
        fs, fe = int(f[1]), int(f[2])
        for ins in (fs + 4, fe - 5):
            p = ins - start
            if 0 <= p < L:
                cov[cell][p] += 1
    if args.smooth > 1:
        k = np.ones(args.smooth) / args.smooth
        for c in order:
            cov[c] = np.convolve(cov[c], k, mode="same")
    for c in order:
        cov[c] *= 1000.0 / max(n_cells[c], 1)
    safe = {c: c.replace(" ", "_") for c in order}
    for c in order:
        path = os.path.join(args.outdir, f"{safe[c]}.bw")
        bw = pyBigWig.open(path, "w")
        bw.addHeader([(chrom, HG38[chrom])])
        bw.addEntries(chrom, start, values=cov[c].astype("float64"), span=1, step=1)
        bw.close()
        print(f"  wrote {path}  (max {cov[c].max():.1f})")
    print(f"Done (locus {args.region}).")


def run_genome_wide(args):
    os.makedirs(args.outdir, exist_ok=True)
    bin_sz = args.bin
    frags = pysam.TabixFile(args.fragments)
    ct, order, n_cells = _meta(args)
    safe = {c: c.replace(" ", "_") for c in order}
    chroms = [c for c in HG38 if c in set(frags.contigs)]
    writers = {}
    for c in order:
        bw = pyBigWig.open(os.path.join(args.outdir, f"{safe[c]}.bw"), "w")
        bw.addHeader([(ch, HG38[ch]) for ch in chroms])
        writers[c] = bw
    for ch in chroms:
        nb = HG38[ch] // bin_sz + 1
        cov = {c: np.zeros(nb, dtype=np.float32) for c in order}
        n = 0
        for row in frags.fetch(ch):
            f = row.split("\t"); cell = ct.get(f[3])
            if cell is None:
                continue
            fs, fe = int(f[1]), int(f[2])
            cov[cell][(fs + 4) // bin_sz] += 1
            cov[cell][(fe - 5) // bin_sz] += 1
            n += 1
        for c in order:
            v = cov[c] * (1000.0 / max(n_cells[c], 1))
            nz = np.nonzero(v)[0]
            if len(nz):
                writers[c].addEntries(
                    ch, (nz * bin_sz).astype(int).tolist(),
                    values=v[nz].astype("float64").tolist(),
                    span=bin_sz)
        print(f"  {ch}: {n:,} fragments")
    for c in order:
        writers[c].close()
    print(f"Done (genome-wide, {bin_sz} bp bins): {len(order)} bigWigs in {args.outdir}/")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5mu", default="workshop_data/pbmc_10k_multiome_workshop.h5mu")
    ap.add_argument("--fragments", default="workshop_data/raw/atac_fragments.tsv.gz")
    ap.add_argument("--region", help="chrom:start-end (locus mode)")
    ap.add_argument("--genome-wide", action="store_true")
    ap.add_argument("--bin", type=int, default=25, help="bin size (genome-wide)")
    ap.add_argument("--outdir", default="workshop_data/tracks")
    ap.add_argument("--smooth", type=int, default=50, help="smoothing (locus mode)")
    args = ap.parse_args()
    if args.genome_wide:
        run_genome_wide(args)
    elif args.region:
        run_locus(args)
    else:
        ap.error("specify --genome-wide or --region")


if __name__ == "__main__":
    main()
