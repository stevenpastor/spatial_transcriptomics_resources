# Visium HD Spatial Transcriptomics Tutorial

[Open In Colab](https://colab.research.google.com/github/stevenpastor/spatial_transcriptomics_resources/blob/main/notebook/visium_hd_crc_p1_tutorial.ipynb)

Analyze a published 10x Genomics Visium HD spatial transcriptomics dataset at one resolution with QC, filtering, cell type annotation, neighborhood analysis, spatially variable genes, and ligand-receptor communication.

## Quick Start

Click the **Colab** badge above. The notebook auto-installs dependencies and downloads pre-processed data (~150 MB) from Figshare. No local setup needed.

## Dataset

Human Colorectal Cancer (CRC), Patient 1 from the [Nature Genetics 2025 publication](https://www.nature.com/articles/s41588-025-02193-3). Pre-processed to 50K bins at 8 µm resolution with ground-truth cell type annotations.

## What you will learn

1. Load and explore Visium HD spatial transcriptomics data
2. Quality control with filtering (per-bin metrics + spatial outlier detection)
3. Compare marker-based annotations against published ground-truth labels, for your convenience
4. Spatial neighborhood analysis, spatially variable genes, cell-cell communication

## 🕹️ Bonus: Game Tutorial

[Open In Colab: Game](https://colab.research.google.com/github/stevenpastor/spatial_transcriptomics_resources/blob/main/notebook/visium_hd_game_tutorial.ipynb)

A playful companion notebook. Hunt seven retro arcade sprites (Pac-Man, Ghost, Space Invader, ...) hidden inside a synthetic tissue; each one unlocks when you correctly perform a core analysis step (QC, clustering, annotation, neighborhoods, SVGs). Runs on Colab in under a minute, no external data download.

## Repository structure

```
notebook/
  visium_hd_crc_p1_tutorial.ipynb   # Main tutorial (run on Colab)
  visium_hd_game_tutorial.ipynb     # Game: find hidden sprites by doing the steps right
scripts/
  utils.py                          # Helper functions used by the main notebook
  game_utils.py                     # Synthetic-tissue generator + validators for the game
  generate_precomputed.py           # Script used to create the precomputed data
```

