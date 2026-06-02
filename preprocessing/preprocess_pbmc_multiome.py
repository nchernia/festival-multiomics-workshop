#!/usr/bin/env python3
"""Pre-process PBMC 10k Multiome data for the Festival of Genomics workshop.

Downloads the 10x Genomics PBMC 10k Multiome v2 dataset, processes RNA and
ATAC modalities (snapatac2 spectral + MACS3 per cell type), derives per-gene
TSS coordinates directly from a GENCODE GTF, and saves everything as a
single .h5mu file. A separate post-step (refresh_alphagenome.py) populates
the AlphaGenome CHIP-TF cache for the MS4A1 enhancer.

The workshop notebook's peak-to-gene
analysis is paired Spearman correlation over a ±250 kb candidate window,
which only needs gene TSS coordinates.

Run time: ~45-60 minutes (MACS3 per-cell-type is the bottleneck).
Requires: ~16 GB RAM, internet access.

Usage:
    python preprocess_pbmc_multiome.py --output-dir ./workshop_data
"""

import argparse
import gzip
import os
import re
import sys
import urllib.request

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc


# ---------------------------------------------------------------------------
# 1. Download 10x PBMC 10k Multiome data
# ---------------------------------------------------------------------------
DATA_URLS = {
    "filtered_feature_bc_matrix.h5": (
        "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/"
        "pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
    ),
    "atac_fragments.tsv.gz": (
        "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/"
        "pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
    ),
    "atac_fragments.tsv.gz.tbi": (
        "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/"
        "pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
    ),
    "atac_peak_annotation.tsv": (
        "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/"
        "pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv"
    ),
}

GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_44/gencode.v44.basic.annotation.gtf.gz"
)


def _http_download(url, path):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(path, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)


def download_data(output_dir):
    raw_dir = os.path.join(output_dir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    for filename, url in DATA_URLS.items():
        path = os.path.join(raw_dir, filename)
        if os.path.exists(path):
            print(f"  Already exists: {filename}")
            continue
        print(f"  Downloading {filename} ...")
        _http_download(url, path)
        print(f"    Done ({os.path.getsize(path) / 1e6:.1f} MB)")
    return raw_dir


# ---------------------------------------------------------------------------
# 2. RNA processing
# ---------------------------------------------------------------------------
def process_rna(rna):
    print("\n--- Processing RNA ---")
    print(f"  Input: {rna.n_obs} cells x {rna.n_vars} genes")

    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], inplace=True)
    rna = rna[rna.obs.pct_counts_mt < 20].copy()

    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)

    sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat")
    sc.tl.pca(rna, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(rna, n_pcs=30)
    sc.tl.leiden(rna, resolution=1.0)
    sc.tl.umap(rna)

    print(f"  Output: {rna.n_obs} cells, {rna.obs['leiden'].nunique()} clusters")
    return rna


# ---------------------------------------------------------------------------
# 3. ATAC processing (snapatac2: tile matrix -> spectral)
# ---------------------------------------------------------------------------
ATAC_MIN_FRAGMENTS = 1000
ATAC_MIN_TSSE      = 7.0


def process_atac(fragments_path, common_cells, output_dir):
    import snapatac2 as snap

    print("\n--- Processing ATAC (snapatac2 tile pipeline) ---")
    snap_dir = os.path.join(output_dir, "snapatac2")
    os.makedirs(snap_dir, exist_ok=True)
    backed_path = os.path.join(snap_dir, "atac_tiles.h5ad")
    if os.path.exists(backed_path):
        os.remove(backed_path)

    atac = snap.pp.import_fragments(
        fragments_path,
        chrom_sizes=snap.genome.hg38,
        file=backed_path,
        min_num_fragments=ATAC_MIN_FRAGMENTS,
        whitelist=list(common_cells),
        sorted_by_barcode=False,
    )
    print(f"  Imported: {atac.n_obs} cells (>={ATAC_MIN_FRAGMENTS} fragments)")

    snap.metrics.tsse(atac, snap.genome.hg38)
    n_before = atac.n_obs
    snap.pp.filter_cells(atac, min_counts=ATAC_MIN_FRAGMENTS, min_tsse=ATAC_MIN_TSSE)
    print(f"  After TSSe>={ATAC_MIN_TSSE} filter: {atac.n_obs} cells "
          f"({n_before - atac.n_obs} dropped)")

    snap.pp.add_tile_matrix(atac, bin_size=500)
    snap.pp.select_features(atac, n_features=50_000)

    n_before = atac.n_obs
    snap.pp.scrublet(atac)
    snap.pp.filter_doublets(atac)
    print(f"  After doublet removal: {atac.n_obs} cells "
          f"({n_before - atac.n_obs} doublets dropped)")

    snap.tl.spectral(atac, n_comps=30)
    snap.tl.umap(atac, use_rep="X_spectral")
    print(f"  Spectral embedding ready: {atac.n_obs} cells x {atac.n_vars} bins")
    return atac


# ---------------------------------------------------------------------------
# 4. MACS3 per cell type -> consensus peaks -> peak matrix
# ---------------------------------------------------------------------------
def call_peaks_per_celltype(atac_snap, cell_type_per_cell, output_dir):
    import snapatac2 as snap
    import scipy.sparse as sp

    print("\n--- Per-cell-type peak calling (MACS3) ---")

    barcodes = list(atac_snap.obs_names)
    labels = [cell_type_per_cell.get(b, None) for b in barcodes]
    keep_mask = np.array([l is not None for l in labels])
    print(f"  Cells with cell_type: {keep_mask.sum()} / {len(barcodes)}")
    if not keep_mask.all():
        atac_snap.subset(obs_indices=np.where(keep_mask)[0])
        barcodes = [b for b, k in zip(barcodes, keep_mask) if k]
        labels   = [l for l, k in zip(labels, keep_mask) if k]

    atac_snap.obs["cell_type"] = labels

    print(f"  Running MACS3 per cell type ({len(set(labels))} types) ...")
    snap.tl.macs3(atac_snap, groupby="cell_type", n_jobs=4)
    per_type_peaks = atac_snap.uns["macs3"]
    for ct, df in per_type_peaks.items():
        print(f"    {ct}: {len(df)} peaks")

    merged = snap.tl.merge_peaks(per_type_peaks, snap.genome.hg38, half_width=250)
    print(f"  Consensus peaks: {len(merged)}")
    atac_snap.uns["peaks_merged"] = merged

    peak_h5 = os.path.join(output_dir, "snapatac2", "atac_peaks.h5ad")
    if os.path.exists(peak_h5):
        os.remove(peak_h5)
    peak_mat = snap.pp.make_peak_matrix(
        atac_snap, use_rep="peaks_merged", file=peak_h5,
    )
    print(f"  Peak matrix: {peak_mat.n_obs} cells x {peak_mat.n_vars} peaks")

    X = peak_mat.X[:].astype(np.float32)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    atac = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=list(peak_mat.obs_names)),
        var=pd.DataFrame(index=list(peak_mat.var_names)),
    )
    if "X_spectral" in atac_snap.obsm:
        atac.obsm["X_spectral"] = np.asarray(atac_snap.obsm["X_spectral"][:])
    if "X_umap" in atac_snap.obsm:
        atac.obsm["X_umap"] = np.asarray(atac_snap.obsm["X_umap"][:])
    atac.obs["cell_type"] = pd.Categorical(labels)
    return atac


# ---------------------------------------------------------------------------
# 5. Hierarchical cell-type annotation
# ---------------------------------------------------------------------------
BROAD_MARKERS = {
    "T cell":         ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"],
    "NK cell":        ["NKG7", "GNLY", "KLRD1", "PRF1", "TYROBP"],
    "B cell":         ["MS4A1", "CD79A", "BANK1", "CD74", "HLA-DRA"],
    "CD14 Monocyte":  ["CD14", "VCAN", "FCN1", "S100A9", "LYZ"],
    "CD16 Monocyte":  ["FCGR3A", "MS4A7", "LST1", "IFITM3", "SAT1"],
    "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "CST3", "HLA-DRA"],
}
T_SUBTYPE_MARKERS = {
    "CD4 T cell": ["IL7R", "LTB", "LEF1", "TCF7", "MAL"],
    "CD8 T cell": ["CD8A", "CD8B", "CCL5", "GZMK"],
}
PLASMA_HI = ["MZB1", "XBP1", "JCHAIN"]
BCELL_LO  = ["MS4A1", "CD79A", "CD79B", "PAX5"]
PLASMA_THRESHOLD_HI = 0.7
PLASMA_THRESHOLD_LO = 0.7


def _score_marker_dict(rna, marker_dict):
    for ct, genes in marker_dict.items():
        avail = [g for g in genes if g in rna.var_names]
        col = f"score_{ct}"
        if len(avail) >= 3:
            sc.tl.score_genes(rna, gene_list=avail, score_name=col)
        else:
            sub = rna[:, avail].X
            sub = sub.toarray() if hasattr(sub, "toarray") else sub
            rna.obs[col] = sub.mean(axis=1)


def _mean_expr(rna, genes):
    rows = []
    for g in genes:
        if g not in rna.var_names: continue
        x = rna[:, g].X
        rows.append(x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel())
    return np.mean(rows, axis=0) if rows else np.zeros(rna.n_obs)


def annotate_cells(rna):
    """Hierarchical: leiden cluster -> broad lineage -> CD4/CD8 split for T cells.
    Then a per-cell pass splits the B-cell group into B cell + Plasma cell by
    marker expression (high MZB1/XBP1/JCHAIN, low MS4A1/CD79A/...).
    """
    import scipy.sparse as sp
    rna_t = rna.copy()
    _score_marker_dict(rna_t, BROAD_MARKERS)

    clusters = sorted(rna_t.obs["leiden"].unique(), key=int)
    broad_mat = pd.DataFrame(index=clusters, columns=list(BROAD_MARKERS), dtype=float)
    for cl in clusters:
        m = rna_t.obs["leiden"] == cl
        for ct in BROAD_MARKERS:
            broad_mat.loc[cl, ct] = rna_t.obs.loc[m, f"score_{ct}"].mean()
    cluster_to_broad = broad_mat.idxmax(axis=1).to_dict()
    rna_t.obs["broad_type"] = rna_t.obs["leiden"].map(cluster_to_broad)

    final = dict(cluster_to_broad)

    # T-cell subtypes
    t_mask = rna_t.obs["broad_type"] == "T cell"
    if t_mask.sum() > 0:
        t_rna = rna_t[t_mask].copy()
        sc.pp.neighbors(t_rna, n_pcs=30)
        sc.tl.leiden(t_rna, resolution=1.0, key_added="t_sub")
        cd8_genes = [g for g in ("CD8A", "CD8B") if g in t_rna.var_names]
        X = t_rna[:, cd8_genes].X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        cd8_expr = X.mean(axis=1)
        sub_to_label = {}
        for sub in sorted(t_rna.obs["t_sub"].unique(), key=int):
            m = (t_rna.obs["t_sub"] == sub).values
            sub_to_label[sub] = "CD8 T cell" if cd8_expr[m].mean() > 0.3 else "CD4 T cell"
        t_rna.obs["t_label"] = t_rna.obs["t_sub"].map(sub_to_label)
        for cl in [c for c in clusters if cluster_to_broad[c] == "T cell"]:
            cl_mask = (t_rna.obs["leiden"] == cl).values
            labels = t_rna.obs.loc[cl_mask, "t_label"]
            final[cl] = labels.mode().iloc[0] if len(labels) else "CD4 T cell"

    cell_type = rna_t.obs["leiden"].map(final).astype("category")
    cell_type_map = dict(zip(rna_t.obs_names, cell_type))

    # Per-cell plasma-cell split out of B cell
    plasma_score = _mean_expr(rna_t, PLASMA_HI)
    bcell_score  = _mean_expr(rna_t, BCELL_LO)
    is_plasma = ((cell_type.values == "B cell") &
                 (plasma_score > PLASMA_THRESHOLD_HI) &
                 (bcell_score  < PLASMA_THRESHOLD_LO))
    for bc, p in zip(rna_t.obs_names, is_plasma):
        if p:
            cell_type_map[bc] = "Plasma cell"
    print(f"  Plasma cells split out of B cell: {int(is_plasma.sum())}")

    return cell_type_map


# ---------------------------------------------------------------------------
# 6. Per-gene TSS from a GENCODE GTF (no GLUE needed)
# ---------------------------------------------------------------------------
GTF_FEATURE_GENE = "gene"
_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def derive_gene_tss(rna_var_names, output_dir):
    """Parse GENCODE GTF (gene-feature lines only), return per-gene TSS dict.

    Replaces the old GLUE-based gene-coordinate step. ~5 minutes for the
    GENCODE v44 basic annotation, including download.
    """
    print("\n--- Deriving per-gene TSS from GENCODE GTF ---")
    gtf_path = os.path.join(output_dir, "gencode.v44.basic.annotation.gtf.gz")
    if not os.path.exists(gtf_path):
        print(f"  Downloading {GTF_URL} ...")
        urllib.request.urlretrieve(GTF_URL, gtf_path)

    rna_set = set(rna_var_names)
    seen = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != GTF_FEATURE_GENE:
                continue
            chrom, start_s, end_s, strand = fields[0], fields[3], fields[4], fields[6]
            start, end = int(start_s) - 1, int(end_s)
            attrs = dict(_ATTR_RE.findall(fields[8]))
            gene = attrs.get("gene_name")
            if gene and gene in rna_set and gene not in seen:
                tss = end if strand == "-" else start
                seen[gene] = (chrom, int(tss), strand)
    print(f"  TSS resolved for {len(seen)} / {len(rna_var_names)} genes")
    return {
        "gene":   np.asarray(list(seen.keys()), dtype=object),
        "chrom":  np.asarray([v[0] for v in seen.values()], dtype=object),
        "tss":    np.asarray([v[1] for v in seen.values()], dtype="int64"),
        "strand": np.asarray([v[2] for v in seen.values()], dtype=object),
    }


# ---------------------------------------------------------------------------
# 7. Assemble + save MuData
# ---------------------------------------------------------------------------
def assemble_and_save(rna, atac, gene_tss, output_dir):
    print("\n--- Saving workshop data ---")
    for mod in (rna, atac):
        if "counts" in mod.layers:
            del mod.layers["counts"]

    mdata = mu.MuData({"rna": rna, "atac": atac})
    mdata.uns["gene_tss"] = gene_tss

    output_path = os.path.join(output_dir, "pbmc_10k_multiome_workshop.h5mu")
    mdata.write(output_path, compression="gzip")
    print(f"  Saved: {output_path}")
    print(f"  RNA:  {rna.n_obs} cells x {rna.n_vars} genes")
    print(f"  ATAC: {atac.n_obs} cells x {atac.n_vars} peaks")
    print(f"  File size: {os.path.getsize(output_path) / 1e6:.1f} MB")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Preprocess PBMC 10k Multiome")
    parser.add_argument("--output-dir", default="./workshop_data")
    parser.add_argument("--skip-alphagenome", action="store_true",
                        help="Skip the AlphaGenome CHIP-TF cache post-step")
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60); print("STEP 1: Download data"); print("=" * 60)
    raw_dir = download_data(output_dir)

    print("\n" + "=" * 60); print("STEP 2: Load 10x RNA matrix"); print("=" * 60)
    h5_path = os.path.join(raw_dir, "filtered_feature_bc_matrix.h5")
    mdata_raw = mu.read_10x_h5(h5_path)
    mdata_raw.var_names_make_unique()
    rna = mdata_raw.mod["rna"].copy()

    print("\n" + "=" * 60); print("STEP 3: Process RNA"); print("=" * 60)
    rna = process_rna(rna)

    print("\n" + "=" * 60); print("STEP 4: Annotate cells"); print("=" * 60)
    cell_type_map = annotate_cells(rna)
    rna.obs["cell_type"] = pd.Categorical([cell_type_map[b] for b in rna.obs_names])
    print("  Per-type totals:")
    print(rna.obs["cell_type"].value_counts().to_string())

    print("\n" + "=" * 60); print("STEP 5: Process ATAC (snapatac2)"); print("=" * 60)
    fragments_path = os.path.join(raw_dir, "atac_fragments.tsv.gz")
    atac_snap = process_atac(fragments_path, rna.obs_names, output_dir)

    print("\n" + "=" * 60); print("STEP 6: MACS3 per cell type"); print("=" * 60)
    atac = call_peaks_per_celltype(atac_snap, cell_type_map, output_dir)

    common = rna.obs_names.intersection(atac.obs_names)
    rna = rna[common].copy()
    atac = atac[common].copy()
    print(f"  Aligned RNA + ATAC: {rna.n_obs} cells")

    print("\n" + "=" * 60); print("STEP 7: Gene TSS from GTF"); print("=" * 60)
    gene_tss = derive_gene_tss(rna.var_names, output_dir)

    print("\n" + "=" * 60); print("STEP 8: Save workshop .h5mu"); print("=" * 60)
    assemble_and_save(rna, atac, gene_tss, output_dir)

    if not args.skip_alphagenome:
        print("\n" + "=" * 60); print("STEP 9: AlphaGenome CHIP-TF cache"); print("=" * 60)
        import subprocess
        h5mu_path = os.path.join(output_dir, "pbmc_10k_multiome_workshop.h5mu")
        refresh = os.path.join(os.path.dirname(__file__), "refresh_alphagenome.py")
        if not os.environ.get("ALPHA_GENOME_API_KEY"):
            print("  ALPHA_GENOME_API_KEY not set — skipping. Run later with:")
            print(f"    ALPHA_GENOME_API_KEY=<key> python {refresh} --h5mu {h5mu_path}")
        else:
            subprocess.run([sys.executable, refresh, "--h5mu", h5mu_path], check=False)
    else:
        print("\n  Skipping AlphaGenome (--skip-alphagenome)")

    print("\n" + "=" * 60); print("PREPROCESSING COMPLETE"); print("=" * 60)
    print(f"\nWorkshop data in {output_dir}/")
    for f in sorted(os.listdir(output_dir)):
        fp = os.path.join(output_dir, f)
        if os.path.isfile(fp):
            print(f"  {f} ({os.path.getsize(fp) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
