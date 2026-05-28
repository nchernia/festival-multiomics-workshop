# Hands-On Multi-Omics Workshop

Festival of Genomics workshop: integrating scRNA-seq + scATAC-seq with deep learning, ~45 minutes.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/USER/REPO/blob/main/workshop_multiomics_integration.ipynb)

> Replace `USER/REPO` in the badge above with your GitHub path before publishing.

## What students learn

The central skill is **peak-to-gene linking**: using the fact that multiome
measures RNA and ATAC in the *same cell* to ask whether a peak's accessibility
tracks a gene's expression across cells. This is what Signac (`LinkPeaks`),
ArchR (`addPeak2GeneLinks`), and Seurat compute — and it's the analysis that
*only* paired multiome enables.

The workshop is split into two notebooks:

**`workshop_multiomics_integration.ipynb`** — the main ~45-min session:

| Part | Tool | Output |
|---|---|---|
| Intro: what a multiome experiment measures | — | Conceptual framing (figure placeholder) |
| 1. scATAC + QC (~5 min) | snapatac2 `tsse` / `n_fragment` | ATAC UMAP by cell type, TSS enrichment, fragments |
| 2. Same cells (~2 min) | RNA + ATAC UMAPs, shared labels | identity is written in both layers (paired data — no integration model needed) |
| 3. Peak-to-gene & DORCs (~15 min) | **paired Spearman correlation** | genome-wide DORC discovery → zoom into MS4A1's distal B-cell enhancer; IGV browser view of per-cell-type tracks |
| 4. Which TF binds it (~6 min) | AlphaGenome CHIP-TF | TF binding predicted from sequence — EBF1/PAX5 (B-master TFs) at the enhancer |

**`rna_annotation.ipynb`** — self-paced sidebar: how the RNA cell-type labels
were derived (module scoring → hierarchical lineage → CD4/CD8 split → DEG
validation). The main notebook consumes the result.

## Repository layout

```
.
├── workshop_multiomics_integration.ipynb      ← main workshop notebook
├── rna_annotation.ipynb                        ← self-paced cell-annotation walkthrough
├── preprocessing/
│   ├── preprocess_pbmc_multiome.py            ← main pipeline (heavy: download, QC, MACS3, GLUE)
│   ├── make_pseudobulk_bigwigs.py             ← per-cell-type bigWigs for the IGV view
│   └── refresh_alphagenome.py                 ← AlphaGenome CHIP-TF cache (needs API key)
├── workshop_data/
│   ├── pbmc_10k_multiome_workshop.h5mu        ← single file students download
│   └── tracks/*.bw                            ← per-cell-type pseudobulk bigWigs (host for IGV)
├── README.md
└── STUDENT_SETUP.md                           ← pre-workshop instructions for students
```

Instructor post-steps (run once, after the main pipeline):
`make_pseudobulk_bigwigs.py` (IGV tracks) and `refresh_alphagenome.py` (TF-binding
cache). Everything they produce is baked into the `.h5mu` / `tracks/` so the
student notebook only loads and plots — nothing heavy runs live.

## For students

See [`STUDENT_SETUP.md`](STUDENT_SETUP.md). TL;DR:

1. Open `workshop_multiomics_integration.ipynb` in Colab via the badge above.
2. Run the first cell (`pip install ...`) — ~2 minutes.
3. The data file will auto-download (a few hundred MB) — no AlphaGenome API key required.
4. Optional: an AlphaGenome key lets the TF-footprinting step (Part 4) run live; otherwise pre-computed predictions are used.
5. `rna_annotation.ipynb` is an optional self-paced walkthrough of how the cell labels were derived.

## For instructors

The notebook depends on a single file `pbmc_10k_multiome_workshop.h5mu`. Generate it once with:

```bash
cd preprocessing
export ALPHA_GENOME_API_KEY=<your-key>     # for the CD8A CHIP-TF cache (Step 9)
python preprocess_pbmc_multiome.py --output-dir ../workshop_data
```

Pipeline steps: download → RNA QC/cluster → annotate → ATAC (snapatac2 QC +
MACS3 per cell type) → GLUE joint embedding → save (`gene_tss` + `guidance_edges`
ride along in `uns`) → **AlphaGenome CHIP-TF cache for the CD8A enhancer**
(delegated to `refresh_alphagenome_cd8a.py`, runs only if the API key is set).

Costs:
- ~200 MB 10x Genomics RNA + ATAC peak matrix download
- ~2.5 GB 10x ATAC fragments download (needed by snapatac2 / MACS3)
- ~3–5 hours GLUE training (`PairedSCGLUEModel`)
- A few AlphaGenome API calls (~1 min)

ATAC QC thresholds (set in `process_atac`): ≥1000 fragments/cell, TSS
enrichment ≥7, doublets removed via snapatac2's scrublet. Changing these
changes the cell set, so GLUE must be retrained (no `--reuse-glue`).

Output: `workshop_data/pbmc_10k_multiome_workshop.h5mu`. **Host this file at a public URL** (e.g. GCS bucket); students download it on demand.

### Refreshing just the AlphaGenome cache

The TF-binding cache is a standalone post-step on the saved `.h5mu` — no GLUE
retrain needed. It picks the featured MS4A1 distal enhancer (paired correlation
over a ±250 kb TSS window) and caches AlphaGenome CHIP-TF predictions for a
B-cell TF panel (EBF1, PAX5, SPIB, TCF3, POU2F2, …):

```bash
ALPHA_GENOME_API_KEY=<key> python preprocessing/refresh_alphagenome.py \
    --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu
```

### Iteration

If you tweak something and need to re-run preprocessing:

```bash
python preprocess_pbmc_multiome.py --output-dir ../workshop_data --reuse-glue   # skip the slow GLUE retrain
python preprocess_pbmc_multiome.py --output-dir ../workshop_data \
    --skip-glue --skip-alphagenome    # smoke-test ATAC + MACS3 only
```

### What's in the .h5mu

```
MuData
├── mod['rna']    cells × ~19k genes   (cell_type, X_pca, X_umap, leiden)
├── mod['atac']   cells × ~180k peaks  (cell_type, tsse, n_fragment, X_spectral, X_umap)
└── uns
    ├── gene_tss          — per-gene (chrom, TSS, strand) from the GTF; defines
    │                       the ±250 kb candidate window for peak-to-gene
    ├── dorc / dorc_links / dorc_params — genome-wide DORC table (per gene: # of
    │                       significantly correlated peaks) + the peak↔gene links
    └── alphagenome_cache — CHIP-TF predictions (B-cell TFs) for the MS4A1 enhancer
```

The notebook computes peak-to-gene links via paired Spearman correlation (peak
accessibility vs gene expression across cells), over peaks in a ±250 kb window
around each gene's TSS (`gene_tss`). No GLUE involved — the GLUE step in
preprocessing now only supplies gene coordinates; its embedding is unused.

(Cell counts depend on the ATAC QC thresholds above — stricter QC than before,
so the saved file is smaller than the original ~11.5k-cell version.)

Raw counts are dropped from the saved file (saves ~1 GB; not used by the notebook).

## Hosting `.h5mu` on GCS (one option)

```bash
gcloud storage cp workshop_data/pbmc_10k_multiome_workshop.h5mu \
    gs://your-bucket/festival-2026/pbmc_10k_multiome_workshop.h5mu \
    --project=YOUR-PROJECT
gcloud storage objects update \
    gs://your-bucket/festival-2026/pbmc_10k_multiome_workshop.h5mu \
    --add-acl-grant=entity=AllUsers,role=READER \
    --project=YOUR-PROJECT
```

Then in cell 3 of the notebook, students set:

```python
DATA_URL = "https://storage.googleapis.com/your-bucket/festival-2026/pbmc_10k_multiome_workshop.h5mu"
```

(or set the `DATA_URL` env-var via `os.environ` before running cell 3).

## Hosting the IGV tracks on GCS

Generate the genome-wide per-cell-type bigWigs (~460 MB total, ~10–20 min), so
IGV is navigable to any gene:

```bash
python preprocessing/make_pseudobulk_bigwigs.py --genome-wide \
    --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu \
    --fragments workshop_data/raw/atac_fragments.tsv.gz --outdir workshop_data/tracks
```

(Or `--region chr11:60430000-60520000` for a single locus, base-pair resolution.)

The embedded IGV view loads these bigWigs by URL. igv.js fetches them with
**HTTP range requests** and cross-origin, so the bucket needs public reads
**and a CORS policy** (without CORS, igv.js fails silently / "data view" errors):

```bash
# upload + make public
gcloud storage cp workshop_data/tracks/*.bw \
    gs://your-bucket/festival-2026/tracks/ --project=YOUR-PROJECT
gcloud storage objects update 'gs://your-bucket/festival-2026/tracks/*' \
    --add-acl-grant=entity=AllUsers,role=READER --project=YOUR-PROJECT

# CORS (required for igv.js)
echo '[{"origin":["*"],"method":["GET","HEAD"],"responseHeader":["*"],"maxAgeSeconds":3600}]' > /tmp/cors.json
gcloud storage buckets update gs://your-bucket --cors-file=/tmp/cors.json --project=YOUR-PROJECT
```

Then before the IGV cell, students set:

```python
os.environ["TRACKS_URL"] = "https://storage.googleapis.com/your-bucket/festival-2026/tracks"
```

(Local file paths don't work for igv.js bigWigs — range requests over the
Jupyter file-comm are unreliable, which is the "data view" error.)
