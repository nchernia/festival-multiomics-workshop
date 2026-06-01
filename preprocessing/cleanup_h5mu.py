#!/usr/bin/env python3
"""One-shot cleanup of the workshop .h5mu before re-uploading to GCS.

Three things:
  1. Split "B cell" into "B cell" and "Plasma cell" via a per-cell marker rule
     (high mean of MZB1/XBP1/JCHAIN, low mean of MS4A1/CD79A/CD79B/PAX5).
  2. Strip unused GLUE residue: varm['X_glue'] on both modalities,
     obs['balancing_weight'] on both modalities.
  3. Save back to disk (gzip-compressed) so the file is ready for GCS re-upload.

Run AFTER refresh_alphagenome.py if you also want the AG cache regenerated
against the new labels; the AG cache only depends on the MS4A1 enhancer
peak coordinates and is independent of the cell labels, so the order
doesn't strictly matter — running cleanup first is fine.

Usage:
    python preprocessing/cleanup_h5mu.py \
        --h5mu pbmc_10k_multiome_workshop.h5mu
"""
import argparse
import warnings
import numpy as np

warnings.filterwarnings("ignore")
import muon as mu


PLASMA_HI = ["MZB1", "XBP1", "JCHAIN"]
BCELL_LO  = ["MS4A1", "CD79A", "CD79B", "PAX5"]
PLASMA_HI_THRESHOLD = 0.7
BCELL_LO_THRESHOLD  = 0.7


def _mean_expr(adata, genes):
    rows = []
    for g in genes:
        if g not in adata.var_names:
            continue
        x = adata[:, g].X
        rows.append(x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel())
    if not rows:
        return np.zeros(adata.n_obs)
    return np.mean(rows, axis=0)


def identify_plasma(rna):
    b_mask = (rna.obs["cell_type"] == "B cell").values
    plasma_score = _mean_expr(rna, PLASMA_HI)
    bcell_score  = _mean_expr(rna, BCELL_LO)
    return b_mask & (plasma_score > PLASMA_HI_THRESHOLD) & (bcell_score < BCELL_LO_THRESHOLD)


def relabel_plasma(mdata):
    rna  = mdata.mod["rna"]
    atac = mdata.mod["atac"]
    is_plasma = identify_plasma(rna)
    n = int(is_plasma.sum())
    print(f"Plasma cells (split from B cell): {n}")

    # Add the new category and assign cells in both modalities.
    for mod_name, mod in (("rna", rna), ("atac", atac)):
        cats = mod.obs["cell_type"]
        if "Plasma cell" not in cats.cat.categories:
            mod.obs["cell_type"] = cats.cat.add_categories(["Plasma cell"])
        mod.obs.loc[is_plasma, "cell_type"] = "Plasma cell"
        # Drop now-empty categories defensively (shouldn't be any, but safe).
        mod.obs["cell_type"] = mod.obs["cell_type"].cat.remove_unused_categories()
        print(f"  {mod_name}.obs.cell_type counts:")
        for ct, c in mod.obs["cell_type"].value_counts().sort_index().items():
            print(f"    {ct}: {c}")


def strip_glue_residue(mdata):
    removed = []
    for mod_name, mod in mdata.mod.items():
        for k in list(mod.varm.keys()):
            if "glue" in k.lower():
                del mod.varm[k]; removed.append(f"{mod_name}.varm['{k}']")
        for k in list(mod.obs.columns):
            if k.startswith("balancing_weight") or k.lower() == "balancing_weight":
                del mod.obs[k]; removed.append(f"{mod_name}.obs['{k}']")
        for k in ("X_umap_glue", "X_glue"):
            if k in mod.obsm:
                del mod.obsm[k]; removed.append(f"{mod_name}.obsm['{k}']")
    if removed:
        print("Stripped GLUE residue:")
        for r in removed:
            print(f"  {r}")
    else:
        print("No GLUE residue found.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5mu", default="pbmc_10k_multiome_workshop.h5mu")
    args = ap.parse_args()

    print(f"Reading {args.h5mu} ...")
    mdata = mu.read_h5mu(args.h5mu)

    print("\n--- 1. Plasma-cell separation ---")
    relabel_plasma(mdata)

    print("\n--- 2. GLUE residue ---")
    strip_glue_residue(mdata)

    print(f"\nWriting back to {args.h5mu} (gzip) ...")
    mdata.write(args.h5mu, compression="gzip")
    print("Done.")


if __name__ == "__main__":
    main()
