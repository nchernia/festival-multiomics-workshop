# Workshop Setup Guide

Read this **before** the workshop. The setup itself is fast (~5 minutes), but the data download (a few hundred MB) takes a few minutes on Colab — don't leave it for the start.

## What you need

- A **Google account** (for Colab).
- A laptop with browser + reasonable internet.
- (Optional) An **AlphaGenome API key** — free, takes ~1 minute to get. See below.

You do **not** need a GPU, a local Python environment, or any genomics tools installed locally. Everything runs in Colab.

## Step 1 — open the notebook in Colab

Click here:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/USER/REPO/blob/main/workshop_multiomics_integration.ipynb)

> *(Instructor: replace `USER/REPO` with your GitHub repo before sending this link.)*

When prompted, sign in with Google and accept the "from GitHub" runtime warning.

## Step 2 — confirm the runtime

In Colab: **Runtime → Change runtime type → Hardware accelerator: CPU** (the default is fine; we don't need a GPU). Click **Connect** in the top-right.

## Step 3 — install dependencies

Run **cell 1** (`%%capture\n!pip install -q ...`). Takes ~2 minutes. You'll see no output — that's expected (`%%capture` hides it). When the cell finishes (green checkmark), move on.

If the install fails, restart the runtime (**Runtime → Restart runtime**) and re-run cell 1.

## Step 4 — set the data URL

Cell 3 downloads the workshop data file (a few hundred MB). Before running it, set the data URL.

In a Colab cell **above** cell 3 (or right at the top), run:

```python
import os
os.environ["DATA_URL"] = "https://storage.googleapis.com/<bucket>/<path>/pbmc_10k_multiome_workshop.h5mu"
```

> *(Instructor: provide the actual URL.)*

Then run cell 3. The download takes 1–3 minutes depending on your connection. You'll see "Downloading from ... Done!" when it's finished.

## Step 5 (optional) — get an AlphaGenome API key

Without a key, **the notebook still runs end-to-end** — Part 3 uses pre-computed predictions baked into the data file. With a key, Part 3 runs **live** against AlphaGenome's servers.

To get a key:

1. Go to <https://deepmind.google.com/science/alphagenome>
2. Click "Get API access" / "Try the API" and follow the form.
3. Free for academic / non-commercial use; takes ~1 minute.

Once you have a key, in Colab:

1. Click the **🔑 key icon** in the left sidebar.
2. Click **+ Add new secret**.
3. Name: `ALPHA_GENOME_API_KEY`. Value: paste your key.
4. Toggle **Notebook access**: ON.

The Part 3 setup cell will detect it automatically and switch to the live API path.

## Step 6 — run the notebook

Use **Runtime → Run all** for the full sweep, or step through cell-by-cell with **Shift+Enter** to follow along.

Total runtime: ~5–10 minutes for all cells (most of that is the data download in cell 3 and the dotplot rendering).

## Troubleshooting

**"No module named 'muon'" / similar import errors after restart**
→ Re-run cell 1. Restarting the runtime wipes installed packages; you have to reinstall.

**Cell 3 fails with `FileNotFoundError`**
→ You haven't set `DATA_URL` and the file isn't in the working directory. Set the env var as in Step 4.

**Part 3 raises `RuntimeError: No API key AND no cached AlphaGenome data`**
→ You're using a `.h5mu` that wasn't preprocessed with the cache step, *and* you don't have an API key. Either ask the instructor for the cached version, or set up an AlphaGenome key (Step 5).

**Out of RAM**
→ The default Colab tier has 12 GB RAM, which is enough. If you somehow hit a limit (e.g. you have many other notebooks open), close them and restart the runtime.

**Anything else**
→ Flag the instructor. Keep going on the parts that work; partial completion is fine.
