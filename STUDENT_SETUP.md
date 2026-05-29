# Workshop Setup Guide

Read this **before** the workshop. The notebook runs in Google Colab — no
local install needed. Two clicks and a short download and you're ready.

## What you need

- A **Google account** (for Colab).
- A laptop with a browser and reasonable internet.
- *(Optional)* An **AlphaGenome API key** — free, ~1 min to get. With a key
  Part 4 queries AlphaGenome live for any gene you re-target. Without one
  you still get the workshop's featured example (MS4A1) from a baked cache.

You do **not** need a GPU, a local Python install, or any genomics tools.

## Step 1 — open the notebook in Colab

Click here:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nchernia/festival-multiomics-workshop/blob/main/workshop_multiomics_integration.ipynb)

> *(Instructor: replace `nchernia/festival-multiomics-workshop` with your GitHub repo before sending this link.)*

Sign in with Google and accept the "from GitHub" runtime warning.

## Step 2 — confirm the runtime

In Colab: **Runtime → Change runtime type → Hardware accelerator: CPU** (the
default is fine; we don't need a GPU). Click **Connect** in the top-right.

## Step 3 — run all cells

**Runtime → Run all.**

- The first cell installs dependencies (`%%capture` hides the output — that's
  expected; a green checkmark when it finishes is what to look for).
- The load-data cell downloads `pbmc_10k_multiome_workshop.h5mu` (~900 MB)
  from the workshop's public URL. Takes 1–3 minutes on a typical connection.
- Everything after that runs in seconds — the IGV browser tracks stream by
  range request, the AlphaGenome panel reads precomputed predictions.

Total runtime: ~5–10 minutes.

## (Optional) AlphaGenome API key for live re-targeting

Without a key, Part 4 (which TF binds the enhancer) shows precomputed
predictions for the featured MS4A1 enhancer. With a key, you can change the
gene at the top of Part 3 and re-run — the AlphaGenome panel queries live
for the new peak.

To get a key:

1. Go to <https://deepmind.google.com/science/alphagenome>.
2. Click "Get API access" / "Try the API" and fill in the (~1-min) form.
3. Free for academic / non-commercial use.

Once you have a key, in Colab:

1. Click the 🔑 **key icon** in the left sidebar.
2. Click **+ Add new secret**.
3. Name: `ALPHA_GENOME_API_KEY`. Value: paste your key.
4. Toggle **Notebook access**: ON.

The AlphaGenome cell detects it automatically.

## Troubleshooting

**`No module named 'muon'`** (or similar) after a restart
→ Re-run the install cell. Restarting the runtime wipes installed packages.

**The data download stalls or fails**
→ Restart the runtime and try the load cell again; the workshop file is on
public GCS so transient network issues are usually all it is.

**The IGV browser shows tracks but nothing displays**
→ Trust the notebook (Colab will prompt) and try again. The IGV widget
loads bigWigs by URL with range requests; if it still fails, flag the
instructor.

**Part 4 (AlphaGenome) raises a missing-cache error**
→ You're running a notebook that wasn't preprocessed with the cache step
*and* you don't have an API key. Either ask the instructor for the cached
version of the `.h5mu`, or set up an AlphaGenome key (above).

**Out of RAM**
→ The free Colab tier has 12 GB which is enough. If you somehow hit a limit,
close other notebooks and restart.

**Anything else**
→ Flag the instructor. Keep going on the parts that work; partial
completion is fine.
