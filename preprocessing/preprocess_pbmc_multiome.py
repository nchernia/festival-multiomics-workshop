#!/usr/bin/env python3
"""
Pre-process PBMC 10k Multiome data for the Festival of Genomics workshop.

Downloads the 10x Genomics PBMC 10k multiome dataset, processes RNA and
ATAC modalities (snapatac2 spectral + MACS3 per cell type), trains a
GLUE integration model, pre-computes AlphaGenome predictions for the
canonical monocyte peak, and saves everything as a single .h5mu file.

Run time: ~3-5 hours depending on hardware (GLUE dominates).
Requires: ~16GB RAM, internet access.

Usage:
    python preprocess_pbmc_multiome.py --output-dir ./workshop_data
"""

import argparse
import os
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
# URLs for the 10x Genomics PBMC 10k Multiome v2 dataset
# These are the direct download links from the 10x datasets page:
# https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0
# Files this pipeline consumes:
#   - filtered_feature_bc_matrix.h5: paired RNA + ATAC count matrices (~200 MB)
#   - atac_fragments.tsv.gz (+ .tbi): fragments, needed by snapatac2 for
#       tile-matrix processing and MACS3 peak calling per cell type (~2.5 GB)
#   - atac_peak_annotation.tsv: CellRanger peak annotations (small; nice-to-have)
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


def _http_download(url, path):
    """urlretrieve that sets a User-Agent (10x's CDN blocks the urllib default)."""
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as r, open(path, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)


def download_data(output_dir):
    """Download raw data files if not already present."""
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
# 2. Process RNA modality
# ---------------------------------------------------------------------------
def process_rna(rna):
    """Standard scanpy RNA processing pipeline."""
    print("\n--- Processing RNA ---")
    print(f"  Input: {rna.n_obs} cells x {rna.n_vars} genes")

    # Basic QC filtering
    sc.pp.filter_cells(rna, min_genes=200)
    sc.pp.filter_genes(rna, min_cells=3)
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], inplace=True)
    rna = rna[rna.obs.pct_counts_mt < 20].copy()

    # Store raw counts for DEG analysis
    rna.layers["counts"] = rna.X.copy()

    # Normalize and log-transform
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)

    # HVGs (seurat flavor; operates on the log-normalized .X — no extra deps)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat")
    sc.tl.pca(rna, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(rna, n_pcs=30)
    sc.tl.leiden(rna, resolution=1.0)
    sc.tl.umap(rna)

    print(f"  Output: {rna.n_obs} cells, {rna.obs['leiden'].nunique()} clusters")
    return rna


# ---------------------------------------------------------------------------
# 3. Process ATAC modality (snapatac2: tile matrix -> spectral)
# ---------------------------------------------------------------------------
# ATAC QC thresholds. Stricter than the snapatac2 defaults: this teaching
# dataset should be clean enough that the downstream story isn't muddied by
# low-quality barcodes or doublets.
ATAC_MIN_FRAGMENTS = 1000   # minimum unique fragments per cell
ATAC_MIN_TSSE      = 7.0    # minimum TSS enrichment (real chromatin signal)


def process_atac(fragments_path, common_cells, output_dir):
    """Process ATAC fragments with snapatac2.

    Loads fragments, applies QC (fragment-count floor, TSS-enrichment floor,
    doublet removal), builds a 500 bp tile matrix, selects features, runs
    spectral dim-reduction + UMAP. Returns a snapatac2 (HDF5-backed) AnnData
    that's ready for per-cell-type peak calling.
    """
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
    print(f"  Imported: {atac.n_obs} cells "
          f"(whitelisted to RNA-passing cells, >={ATAC_MIN_FRAGMENTS} fragments)")

    # QC: TSS enrichment. Real chromatin signal piles up at TSSs; low-TSSe
    # barcodes are ambient/dead. Filter on both fragment count and TSSe.
    snap.metrics.tsse(atac, snap.genome.hg38)
    n_before = atac.n_obs
    snap.pp.filter_cells(
        atac, min_counts=ATAC_MIN_FRAGMENTS, min_tsse=ATAC_MIN_TSSE,
    )
    print(f"  After TSSe>={ATAC_MIN_TSSE} filter: {atac.n_obs} cells "
          f"({n_before - atac.n_obs} dropped)")

    snap.pp.add_tile_matrix(atac, bin_size=500)
    snap.pp.select_features(atac, n_features=50_000)

    # QC: doublet removal (snapatac2's scrublet on the selected tile features).
    n_before = atac.n_obs
    snap.pp.scrublet(atac)
    snap.pp.filter_doublets(atac)
    print(f"  After doublet removal: {atac.n_obs} cells "
          f"({n_before - atac.n_obs} doublets dropped)")

    snap.tl.spectral(atac, n_comps=30)
    snap.tl.umap(atac, use_rep="X_spectral")
    print(f"  After spectral: {atac.n_obs} cells x {atac.n_vars} bins, "
          f"X_spectral{atac.obsm['X_spectral'].shape}")
    return atac


# ---------------------------------------------------------------------------
# 4. Per-cell-type peak calling with MACS3 + consensus + peak matrix
# ---------------------------------------------------------------------------
def call_peaks_per_celltype(atac_snap, cell_type_per_cell, output_dir):
    """Run MACS3 per cell type, merge into consensus peaks, build cell x peak
    matrix. Returns a standard anndata.AnnData (cell x peak, integer counts)
    ready to feed into GLUE.

    cell_type_per_cell: dict{barcode -> str} for cells in atac_snap.obs_names.
    Any cell missing from this dict is dropped (e.g. ATAC-only barcodes).
    """
    import snapatac2 as snap
    import scipy.sparse as sp

    print("\n--- Per-cell-type peak calling (MACS3) ---")

    # Attach cell_type labels — drop unlabeled cells (paired multiome should
    # have a label for every ATAC barcode that's also in RNA).
    barcodes = list(atac_snap.obs_names)
    labels = [cell_type_per_cell.get(b, None) for b in barcodes]
    keep_mask = np.array([l is not None for l in labels])
    print(f"  Cells with cell_type: {keep_mask.sum()} / {len(barcodes)}")
    if not keep_mask.all():
        atac_snap.subset(obs_indices=np.where(keep_mask)[0])
        barcodes = [b for b, k in zip(barcodes, keep_mask) if k]
        labels   = [l for l, k in zip(labels, keep_mask) if k]

    # snapatac2 stores obs as polars; pass a plain list of strings (not pd.Categorical)
    atac_snap.obs["cell_type"] = labels

    print(f"  Running MACS3 per cell type ({len(set(labels))} types) ...")
    snap.tl.macs3(atac_snap, groupby="cell_type", n_jobs=4)
    per_type_peaks = atac_snap.uns["macs3"]
    for ct, df in per_type_peaks.items():
        print(f"    {ct}: {len(df)} peaks")

    merged = snap.tl.merge_peaks(
        per_type_peaks, snap.genome.hg38, half_width=250,
    )
    print(f"  Consensus peaks (merged): {len(merged)}")
    atac_snap.uns["peaks_merged"] = merged

    peak_h5 = os.path.join(output_dir, "snapatac2", "atac_peaks.h5ad")
    if os.path.exists(peak_h5):
        os.remove(peak_h5)
    peak_mat = snap.pp.make_peak_matrix(
        atac_snap, use_rep="peaks_merged", file=peak_h5,
    )
    print(f"  Peak matrix: {peak_mat.n_obs} cells x {peak_mat.n_vars} peaks")

    # Convert to standard anndata for downstream (GLUE, .h5mu serialization)
    X = peak_mat.X[:].astype(np.float32)
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    atac = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=list(peak_mat.obs_names)),
        var=pd.DataFrame(index=list(peak_mat.var_names)),
    )
    # Carry over spectral / UMAP from snap into the new anndata.
    if "X_spectral" in atac_snap.obsm:
        atac.obsm["X_spectral"] = np.asarray(atac_snap.obsm["X_spectral"][:])
    if "X_umap" in atac_snap.obsm:
        atac.obsm["X_umap"] = np.asarray(atac_snap.obsm["X_umap"][:])
    atac.obs["cell_type"] = pd.Categorical(labels)
    atac.layers["counts"] = atac.X.copy()
    return atac


# ---------------------------------------------------------------------------
# 5. Train GLUE integration model
# ---------------------------------------------------------------------------
def train_glue(rna, atac, output_dir, reuse=False):
    """Train GLUE model for RNA-ATAC integration. With reuse=True, load the
    pre-trained model from output_dir/glue_model_weights and skip training."""
    print("\n--- Training GLUE model ---" if not reuse else "\n--- Loading GLUE model ---")

    import scglue

    # Gene annotation: extract gene coordinates from the RNA var
    # For human data, we need genomic coordinates for genes
    print("  Downloading gene annotations ...")
    gtf_url = (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        "release_44/gencode.v44.basic.annotation.gtf.gz"
    )
    gtf_path = os.path.join(output_dir, "gencode.v44.basic.annotation.gtf.gz")
    if not os.path.exists(gtf_path):
        urllib.request.urlretrieve(gtf_url, gtf_path)

    # Annotate genes with genomic coordinates (rna.var_names are gene symbols)
    scglue.data.get_gene_annotation(rna, gtf=gtf_path, gtf_by="gene_name")
    rna = rna[:, rna.var["chromStart"].notna()].copy()
    print(f"  Genes with coordinates: {rna.n_vars}")

    # Stash per-gene TSS for the notebook's peak-to-gene window (no GLUE needed
    # downstream — this is just gene coordinates). TSS = chromEnd on the minus
    # strand, chromStart on the plus strand. Carried to mdata.uns at save time.
    _strand = rna.var["strand"].astype(str).values
    _tss = np.where(_strand == "-",
                    rna.var["chromEnd"].astype("int64").values,
                    rna.var["chromStart"].astype("int64").values)
    rna.uns["gene_tss"] = {
        "gene":   np.asarray(rna.var_names, dtype=object),
        "chrom":  np.asarray(rna.var["chrom"].astype(str).values, dtype=object),
        "tss":    _tss.astype("int64"),
        "strand": np.asarray(_strand, dtype=object),
    }

    # Parse peak coordinates into ATAC var (drop peaks on alt contigs, e.g. KI/GL)
    peak_info = atac.var_names.str.extract(
        r"(?P<chrom>chr\w+)[:\-](?P<chromStart>\d+)[:\-](?P<chromEnd>\d+)"
    )
    n_alt = int(peak_info["chromStart"].isna().sum())
    if n_alt:
        print(f"  Dropping {n_alt} peaks on alt contigs")
        keep = peak_info["chromStart"].notna().values
        atac = atac[:, keep].copy()
        peak_info = peak_info.loc[keep].reset_index(drop=True)
    atac.var["chrom"] = peak_info["chrom"].values
    atac.var["chromStart"] = peak_info["chromStart"].astype(int).values
    atac.var["chromEnd"] = peak_info["chromEnd"].astype(int).values

    # Build guidance graph linking genes to nearby peaks
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    print(f"  Guidance graph: {guidance.number_of_nodes()} nodes, "
          f"{guidance.number_of_edges()} edges")

    # Configure datasets for GLUE. NB likelihood needs raw integer counts —
    # rna.X / atac.X have been normalized, so point at layers["counts"].
    # ATAC representation: X_spectral (from snapatac2) — like LSI but built
    # on the tile matrix rather than CellRanger's peak matrix.
    # CRITICAL: use_obs_names=True is REQUIRED for PairedSCGLUEModel to detect
    # paired cells. Default is False, which generates random UUIDs per cell so
    # pmsk is all False, cos_loss/real_cross_loss stay at 0, and PairedSCGLUE
    # silently degrades to vanilla SCGLUE.
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True, use_rep="X_pca",
        use_layer="counts", use_obs_names=True,
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=False, use_rep="X_spectral",
        use_layer="counts", use_obs_names=True,
    )

    # Sanity-check that paired-cell signal will actually engage before we spend
    # hours training. With matched obs_names + use_obs_names=True, the trainer
    # builds a paired mask covering every cell — cos_loss should be nonzero.
    assert rna.obs_names.equals(atac.obs_names), (
        f"PairedSCGLUEModel needs matching obs_names (got rna={rna.n_obs}, "
        f"atac={atac.n_obs}, intersection={len(set(rna.obs_names) & set(atac.obs_names))})"
    )
    print(f"  Paired cells: {rna.n_obs} (rna.obs_names == atac.obs_names)")

    weights_path = os.path.join(output_dir, "glue_model_weights")
    if reuse:
        if not os.path.exists(weights_path):
            raise FileNotFoundError(
                f"--reuse-glue passed but {weights_path} not found"
            )
        print(f"  Loading saved GLUE model from {weights_path}")
        glue = scglue.models.load_model(weights_path)
    else:
        glue_dir = os.path.join(output_dir, "glue_model")
        # Use PairedSCGLUEModel for multiome (every cell has both modalities).
        # `lam_cos` directly drives paired RNA/ATAC embeddings of the same cell
        # together; bumped from 0.02 default to 0.20 for tighter alignment.
        # Also keep the strong adversarial term and early burn-in.
        glue = scglue.models.fit_SCGLUE(
            {"rna": rna, "atac": atac},
            guidance,
            model=scglue.models.PairedSCGLUEModel,
            init_kws={"random_seed": 0},
            compile_kws={
                "lam_align":       0.20,
                "lam_cos":         0.20,
                "lam_joint_cross": 0.10,
                "lam_real_cross":  0.10,
            },
            fit_kws={
                "directory": glue_dir,
                "align_burnin": 1,
                "safe_burnin": False,
                "max_epochs": 300,
                "patience": 20,
            },
        )

    # Encode both modalities into joint embedding
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)

    if not reuse:
        glue.save(weights_path)
        print(f"  GLUE model saved to {weights_path}")

    # Compute joint UMAP
    combined = ad.concat([rna, atac], label="modality", keys=["rna", "atac"])
    sc.pp.neighbors(combined, use_rep="X_glue")
    sc.tl.umap(combined)

    # Transfer UMAP coordinates back
    n_rna = rna.n_obs
    rna.obsm["X_umap_glue"] = combined.obsm["X_umap"][:n_rna]
    atac.obsm["X_umap_glue"] = combined.obsm["X_umap"][n_rna:]

    # Save guidance graph for peak-gene links (instructor-side artifact)
    import networkx as nx
    nx.write_graphml(guidance, os.path.join(output_dir, "guidance_graph.graphml"))

    # Extract bipartite peak<->gene edges for the workshop notebook
    # (peaks have names like "chr1:12345-67890"; genes are gene symbols).
    peak_set = set(atac.var_names)
    gene_set = set(rna.var_names)
    edges = []
    for u, v, d in guidance.edges(data=True):
        if u in peak_set and v in gene_set:
            peak, gene = u, v
        elif v in peak_set and u in gene_set:
            peak, gene = v, u
        else:
            continue
        edges.append((peak, gene, float(d.get("weight", 1.0))))
    guidance_edges = {
        "peak":   np.array([e[0] for e in edges], dtype=object),
        "gene":   np.array([e[1] for e in edges], dtype=object),
        "weight": np.array([e[2] for e in edges], dtype=np.float32),
    }
    print(f"  Guidance edges (peak<->gene): {len(edges)}")

    print("  GLUE training complete")
    return rna, atac, glue, guidance_edges


# ---------------------------------------------------------------------------
# 6. Pre-compute AlphaGenome results (cached fallback for participants)
# ---------------------------------------------------------------------------
# Subtype-discriminative markers — used for the cluster->celltype assignment.
# Pan-lineage markers like CD3D/CD3E (shared across T-cell subtypes) and NKG7
# (shared between NK and effector CD8) are dropped because they would inflate
# multiple scores equally and prevent CD4/CD8 separation. The notebook's
# DOTPLOT cell uses a richer marker dict; this one is just for assignment.
# Two-pass annotation: broad lineage first, then T-cell subtypes within
# T-cell clusters. Avoids the CD3D/CD3E shared-marker problem.
BROAD_MARKERS = {
    "T cell":         ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"],
    "NK cell":        ["NKG7", "GNLY", "KLRD1", "PRF1", "TYROBP"],
    "B cell":         ["MS4A1", "CD79A", "BANK1", "CD74", "HLA-DRA"],
    "CD14 Monocyte":  ["CD14", "VCAN", "FCN1", "S100A9", "LYZ"],
    "CD16 Monocyte":  ["FCGR3A", "MS4A7", "LST1", "IFITM3", "SAT1"],
    "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "CST3", "HLA-DRA"],
}
T_SUBTYPE_MARKERS = {
    "CD4 T cell":     ["IL7R", "LTB", "LEF1", "TCF7", "MAL"],
    "CD8 T cell":     ["CD8A", "CD8B", "CCL5", "GZMK"],
}


def _score_marker_dict(rna, marker_dict):
    """Compute and store sc.tl.score_genes for each cell type in marker_dict."""
    for ct, genes in marker_dict.items():
        avail = [g for g in genes if g in rna.var_names]
        col = f"score_{ct}"
        if len(avail) >= 3:
            sc.tl.score_genes(rna, gene_list=avail, score_name=col)
        else:
            sub = rna[:, avail].X
            sub = sub.toarray() if hasattr(sub, "toarray") else sub
            rna.obs[col] = sub.mean(axis=1)


def annotate_cells(rna):
    """Hierarchical two-pass annotation. Returns dict {barcode -> cell_type}.

    Pass 1: each leiden cluster -> broad lineage (T/NK/B/Mono/DC) by
            top mean module score across BROAD_MARKERS panels.
    Pass 2: T-cell clusters are re-clustered (leiden on PCA neighbors), each
            subcluster labelled CD4 / CD8 by mean CD8A+CD8B (>0.3 -> CD8).
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
    t_mask = rna_t.obs["broad_type"] == "T cell"
    if t_mask.sum() > 0:
        t_rna = rna_t[t_mask].copy()
        sc.pp.neighbors(t_rna, n_pcs=30)
        sc.tl.leiden(t_rna, resolution=1.0, key_added="t_sub")

        cd8_genes = [g for g in ("CD8A", "CD8B") if g in t_rna.var_names]
        X = t_rna[:, cd8_genes].X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        cd8_expr = X.mean(axis=1)
        t_rna.obs["_cd8_expr"] = cd8_expr

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
    return dict(zip(rna_t.obs_names, cell_type))




# ---------------------------------------------------------------------------
# 7. Assemble and save MuData
# ---------------------------------------------------------------------------
def _sanitize_var(adata):
    """Drop scglue-added columns that aren't needed by the workshop notebook
    and break h5mu serialization (mixed-type object columns like artif_dupl)."""
    drop = ["artif_dupl", "chrom", "chromStart", "chromEnd", "strand",
            "name", "score", "thickStart", "thickEnd", "itemRgb",
            "blockCount", "blockSizes", "blockStarts"]
    for c in drop:
        if c in adata.var.columns:
            del adata.var[c]


def assemble_and_save(rna, atac, output_dir,
                      guidance_edges=None, alphagenome_cache=None):
    """Combine into MuData and save. Auxiliary artifacts (peak<->gene edges,
    AlphaGenome cache) ride along in mdata.uns so participants only need
    this one file."""
    print("\n--- Saving workshop data ---")

    for mod in (rna, atac):
        if "counts" in mod.layers:
            del mod.layers["counts"]
        # scglue's __scglue__ uns dict contains internal training state with
        # 'null'-typed elements that older anndata versions can't read. Strip
        # it so the .h5mu is portable.
        if "__scglue__" in mod.uns:
            del mod.uns["__scglue__"]
    # Carry the per-gene TSS lookup (computed in train_glue) up to mdata.uns so
    # the notebook can build the peak-to-gene window without the GTF or GLUE.
    gene_tss = rna.uns.pop("gene_tss", None)

    _sanitize_var(rna)
    _sanitize_var(atac)

    # Keep cell_type on BOTH modalities. The main notebook consumes these
    # labels directly; the separate rna_annotation.ipynb re-derives them as a
    # teaching exercise (it strips and re-computes its own).
    mdata = mu.MuData({"rna": rna, "atac": atac})
    if gene_tss is not None:
        mdata.uns["gene_tss"] = gene_tss

    if guidance_edges is not None:
        mdata.uns["guidance_edges"] = guidance_edges
    if alphagenome_cache is not None:
        mdata.uns["alphagenome_cache"] = alphagenome_cache

    output_path = os.path.join(output_dir, "pbmc_10k_multiome_workshop.h5mu")
    mdata.write(output_path, compression="gzip")
    print(f"  Saved: {output_path}")
    print(f"  RNA:  {rna.n_obs} cells x {rna.n_vars} genes")
    print(f"  ATAC: {atac.n_obs} cells x {atac.n_vars} peaks")
    if guidance_edges is not None:
        print(f"  Guidance edges: {len(guidance_edges['peak'])}")
    if alphagenome_cache is not None:
        print(f"  AlphaGenome cache: yes ({alphagenome_cache['monocyte_peak']})")
    print(f"  File size: {os.path.getsize(output_path) / 1e6:.1f} MB")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Pre-process PBMC 10k Multiome for workshop"
    )
    parser.add_argument(
        "--output-dir", default="./workshop_data",
        help="Directory for all output files"
    )
    parser.add_argument(
        "--skip-glue", action="store_true",
        help="Skip GLUE training (useful for testing)"
    )
    parser.add_argument(
        "--reuse-glue", action="store_true",
        help="Skip the (slow) GLUE training; load pre-trained weights from "
             "<output-dir>/glue_model_weights and just re-encode."
    )
    parser.add_argument(
        "--skip-alphagenome", action="store_true",
        help="Skip AlphaGenome pre-computation"
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Download
    print("=" * 60)
    print("STEP 1: Download data")
    print("=" * 60)
    raw_dir = download_data(output_dir)

    # Step 2: Load 10x .h5 (just for RNA — ATAC goes through snapatac2)
    print("\n" + "=" * 60)
    print("STEP 2: Load multiome data (RNA from h5)")
    print("=" * 60)
    h5_path = os.path.join(raw_dir, "filtered_feature_bc_matrix.h5")
    mdata = mu.read_10x_h5(h5_path)
    mdata.var_names_make_unique()
    print(f"  Loaded: {mdata}")

    rna = mdata.mod["rna"].copy()

    # Step 3: Process RNA
    print("\n" + "=" * 60)
    print("STEP 3: Process RNA")
    print("=" * 60)
    rna = process_rna(rna)

    # Step 4: Annotate cells (broad lineage + T subcluster) — needed before
    # MACS3 per-cell-type peak calling.
    print("\n" + "=" * 60)
    print("STEP 4: Annotate cells")
    print("=" * 60)
    cell_type_map = annotate_cells(rna)
    rna.obs["cell_type"] = pd.Categorical([cell_type_map[b] for b in rna.obs_names])
    print(f"  Per-type totals:")
    print(rna.obs["cell_type"].value_counts().to_string())

    # Step 5: Process ATAC with snapatac2 (tile matrix -> spectral)
    print("\n" + "=" * 60)
    print("STEP 5: Process ATAC (snapatac2)")
    print("=" * 60)
    fragments_path = os.path.join(raw_dir, "atac_fragments.tsv.gz")
    atac_snap = process_atac(fragments_path, rna.obs_names, output_dir)

    # Step 6: Per-cell-type MACS3 peaks -> consensus -> peak matrix
    print("\n" + "=" * 60)
    print("STEP 6: Per-cell-type peak calling (MACS3)")
    print("=" * 60)
    atac = call_peaks_per_celltype(atac_snap, cell_type_map, output_dir)

    # Align RNA to ATAC's cell set
    common = rna.obs_names.intersection(atac.obs_names)
    rna = rna[common].copy()
    atac = atac[common].copy()
    print(f"  Aligned RNA + ATAC: {rna.n_obs} cells")

    # Step 7: GLUE integration
    guidance_edges = None
    if not args.skip_glue:
        print("\n" + "=" * 60)
        print("STEP 7: GLUE integration")
        print("=" * 60)
        rna, atac, glue, guidance_edges = train_glue(
            rna, atac, output_dir, reuse=args.reuse_glue
        )
    else:
        print("\n  Skipping GLUE (--skip-glue)")

    # Step 8: Save (gene_tss + guidance_edges ride along in mdata.uns)
    print("\n" + "=" * 60)
    print("STEP 8: Save workshop data")
    print("=" * 60)
    assemble_and_save(rna, atac, output_dir, guidance_edges=guidance_edges)

    # Step 9: AlphaGenome CHIP-TF cache for the CD8A enhancer. This is a
    # post-step that operates on the saved .h5mu (it needs gene_tss), so it
    # runs after the save above. Delegated to refresh_alphagenome_cd8a.py.
    if not args.skip_alphagenome:
        print("\n" + "=" * 60)
        print("STEP 9: AlphaGenome CHIP-TF cache (CD8A enhancer)")
        print("=" * 60)
        import subprocess
        h5mu_path = os.path.join(output_dir, "pbmc_10k_multiome_workshop.h5mu")
        refresh = os.path.join(os.path.dirname(__file__), "refresh_alphagenome_cd8a.py")
        if not os.environ.get("ALPHA_GENOME_API_KEY"):
            print("  ALPHA_GENOME_API_KEY not set — skipping. Run later with:")
            print(f"    ALPHA_GENOME_API_KEY=<key> python {refresh} --h5mu {h5mu_path}")
        else:
            subprocess.run([sys.executable, refresh, "--h5mu", h5mu_path], check=False)
    else:
        print("\n  Skipping AlphaGenome (--skip-alphagenome)")

    print("\n" + "=" * 60)
    print("PREPROCESSING COMPLETE")
    print("=" * 60)
    print(f"\nWorkshop data directory: {output_dir}/")
    print("Files:")
    for f in sorted(os.listdir(output_dir)):
        fpath = os.path.join(output_dir, f)
        if os.path.isfile(fpath):
            print(f"  {f} ({os.path.getsize(fpath) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
