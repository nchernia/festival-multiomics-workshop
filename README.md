# Hands-On Multi-Omics Workshop

Festival of Genomics workshop: integrating scRNA-seq + scATAC-seq, ~45 minutes.
Runs entirely in **Google Colab** — no local install, no GPU.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nchernia/festival-multiomics-workshop/blob/main/workshop_multiomics_integration.ipynb)

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
| 1. The puzzle of gene regulation | — | Conceptual framing + the two matrices |
| 2. Two views of the same cells | independent PCA/spectral + Leiden, shared barcodes | Cluster UMAPs (RNA + ATAC), then cell-type reveal |
| 3. Peak-to-gene linking | paired Spearman correlation in ±250 kb window | Candidate enhancers for MS4A1, IGV view, scatter, UMAP overlay |
| 4. Which TF binds it | AlphaGenome CHIP-TF (primary PBMC ontologies + GM12878 fallback) | B-master TFs predicted to bind |

**`rna_annotation.ipynb`** — self-paced sidebar: how the RNA cell-type labels
were derived (module scoring → hierarchical lineage → CD4/CD8 split → plasma
cell separation → DEG validation). The main notebook consumes the result.

## For students

[Open the main notebook in Colab](https://colab.research.google.com/github/nchernia/festival-multiomics-workshop/blob/main/workshop_multiomics_integration.ipynb)
and **Runtime → Run all**. Total time: ~5–10 minutes (most of it is the
`.h5mu` download in the load cell).

See [`STUDENT_SETUP.md`](STUDENT_SETUP.md) for a step-by-step walkthrough and
an optional AlphaGenome API key setup (lets Part 4 query AlphaGenome live for
any peak you re-target).

## Repository layout

```
.
├── workshop_multiomics_integration.ipynb        ← main workshop notebook
├── rna_annotation.ipynb                         ← self-paced cell-annotation walkthrough
├── cell_type_scATAC.png                         ← Part 1 schematic (embedded in the notebook)
├── preprocessing/
│   ├── preprocess_pbmc_multiome.py              ← main pipeline (download, QC, MACS3, GTF, save)
│   ├── cleanup_h5mu.py                          ← one-shot fixer for an existing .h5mu
│   ├── make_pseudobulk_bigwigs.py               ← per-cell-type bigWigs for the IGV view
│   └── refresh_alphagenome.py                   ← AlphaGenome CHIP-TF cache (needs API key)
├── README.md
├── STUDENT_SETUP.md                             ← pre-workshop instructions for students
└── LICENSE
```

Heavy data (`.h5mu`, `tracks/*.bw`) is hosted separately (see below); the
repo holds only code + docs.

## For instructors

The notebook depends on a single hosted `.h5mu` file. Generate it once with:

```bash
cd preprocessing
export ALPHA_GENOME_API_KEY=<your-key>      # optional, only for the AG cache
python preprocess_pbmc_multiome.py --output-dir ../workshop_data
```

Pipeline: download → RNA QC/cluster/annotate (including B-cell vs Plasma-cell
split) → ATAC (snapatac2 QC + MACS3 per cell type) → derive per-gene TSS from
GENCODE GTF → save → **AlphaGenome CHIP-TF cache for the MS4A1 enhancer**
(delegated to `refresh_alphagenome.py`, runs only if `ALPHA_GENOME_API_KEY`
is set).

Costs:
- ~200 MB 10x Genomics RNA + ATAC peak matrix download
- ~2.5 GB 10x ATAC fragments download (needed by snapatac2 / MACS3)
- ~45–60 minutes total wall-clock (MACS3 per-cell-type is the bottleneck)
- A few AlphaGenome API calls (~1 min)

ATAC QC thresholds (in `process_atac`): ≥1,000 fragments/cell, TSS enrichment
≥7, doublets removed via snapatac2's scrublet.

**Note on GLUE.** Earlier versions of this pipeline trained a GLUE model for
joint embedding (3–5 hours). The workshop's peak-to-gene analysis is paired
Spearman correlation in a ±250 kb candidate window — it only needs per-gene
TSS coordinates, which the current pipeline reads directly from the GENCODE
GTF. GLUE is no longer used.

### Cleaning up an existing .h5mu (no re-preprocessing)

If you already have a `.h5mu` produced by an older pipeline and just want to
(a) split plasma cells out of the B-cell label and (b) strip GLUE residue, run:

```bash
python preprocessing/cleanup_h5mu.py --h5mu pbmc_10k_multiome_workshop.h5mu
```

~1 minute. Output overwrites the input.

### Post-steps (after the main pipeline)

```bash
# Per-cell-type bigWigs for the IGV view (~460 MB, ~10–20 min)
python preprocessing/make_pseudobulk_bigwigs.py --genome-wide \
    --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu \
    --fragments workshop_data/raw/atac_fragments.tsv.gz --outdir workshop_data/tracks
# (or `--region chr11:60430000-60520000` for a single-locus, bp-resolution version)

# AlphaGenome CHIP-TF cache for the featured enhancer (needs API key, ~1 min)
ALPHA_GENOME_API_KEY=<key> python preprocessing/refresh_alphagenome.py \
    --h5mu workshop_data/pbmc_10k_multiome_workshop.h5mu
```

### Hosting on GCS

The `.h5mu` and the IGV bigWigs both need a public, range-capable HTTPS
endpoint (CORS-enabled for the IGV tracks):

```bash
BUCKET=your-bucket
PROJ=YOUR-PROJECT

# .h5mu
gcloud storage cp workshop_data/pbmc_10k_multiome_workshop.h5mu \
    gs://$BUCKET/festival-2026/pbmc_10k_multiome_workshop.h5mu --project=$PROJ
gcloud storage objects update \
    gs://$BUCKET/festival-2026/pbmc_10k_multiome_workshop.h5mu \
    --add-acl-grant=entity=AllUsers,role=READER --project=$PROJ

# tracks/
gcloud storage cp workshop_data/tracks/*.bw \
    gs://$BUCKET/festival-2026/tracks/ --project=$PROJ
gcloud storage objects update 'gs://$BUCKET/festival-2026/tracks/*' \
    --add-acl-grant=entity=AllUsers,role=READER --project=$PROJ

# CORS (required for the IGV bigWigs)
echo '[{"origin":["*"],"method":["GET","HEAD"],"responseHeader":["*"],"maxAgeSeconds":3600}]' > /tmp/cors.json
gcloud storage buckets update gs://$BUCKET --cors-file=/tmp/cors.json --project=$PROJ
```

In the notebook:

- The load-data cell already reads `DATA_URL` from `os.environ` with a public
  default pointing at `broad-p16-calico/festival-2026/`.
- The IGV cell reads `TRACKS_URL` the same way.

Workshop default hosting:

- `.h5mu`: `https://storage.googleapis.com/broad-p16-calico/festival-2026/pbmc_10k_multiome_workshop.h5mu`
- `tracks/`: `https://storage.googleapis.com/broad-p16-calico/festival-2026/tracks`

### What's in the .h5mu

```
MuData
├── mod['rna']    cells × ~19k genes   (cell_type, X_pca, X_umap, leiden)
├── mod['atac']   cells × ~180k peaks  (cell_type, tsse, n_fragment, X_spectral, X_umap)
└── uns
    ├── gene_tss          — per-gene (chrom, TSS, strand) from GENCODE v44;
    │                       defines the ±250 kb candidate window for peak-to-gene
    └── alphagenome_cache — CHIP-TF predictions for the MS4A1 enhancer (dedup
                            to one row per TF: the biosample with the strongest
                            predicted binding signal inside the peak ±500 bp)
```

Cell types present in `obs['cell_type']` (both modalities, same labels):
`B cell`, `Plasma cell`, `CD4 T cell`, `CD8 T cell`, `NK cell`, `CD14 Monocyte`,
`CD16 Monocyte`, `Dendritic cell`.

The notebook computes peak-to-gene links via paired Spearman correlation
(peak accessibility vs gene expression across cells), over peaks in a
±250 kb window around each gene's TSS (`gene_tss`).

Raw counts are dropped from the saved file (saves ~1 GB; not used by the
notebook).

## License

MIT — see [`LICENSE`](LICENSE).
