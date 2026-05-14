# Spatial Transcriptomics Resources

Tutorials, tools, and references for spatial transcriptomics analysis. Still under construction. Over time it will grow into a one-stop set of links and resources for spatial work, with the goal of helping you decide how to analyze your own data along with some examples of analysis.

---

## Start here

New to spatial transcriptomics? Read [`intro_to_spatial_background.md`](intro_to_spatial_background/intro_to_spatial_background.md) first. It covers why spatial data matters, how some of the technologies work, and the analysis vocabulary you will encounter in this repo and its tutorials. The current focus is Visium HD but other technologies will soon follow. Once you underdstand the basics, go to the quickstart below.

---

## Tutorials

### Visium HD Spatial Transcriptomics

#### Quick Start

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/stevenpastor/spatial_transcriptomics_resources/blob/main/notebook/visium_hd_quickstart.ipynb)

Walk through a Visium HD analysis on a 10x Genomics dataset: QC, filtering, cell type annotation, neighborhood analysis, spatially variable genes, and ligand-receptor communication.

Click the **Colab** badge above to launch the notebook. Packages install themselves and the pre-processed data (~150 MB) downloads from Figshare on first run. Nothing to set up locally.

##### Dataset

Human Colorectal Cancer (CRC), Patient 1, from the [Nature Genetics 2025 publication](https://www.nature.com/articles/s41588-025-02193-3). The notebook works with ~40K 8 um bins (a random subsample of the full tissue) plus the published ground-truth cell type labels from the original authors.

#### What you will learn

1. How to load and poke around a Visium HD dataset
2. QC and filtering, both per-bin and via spatial outlier detection
3. Annotating bins with marker genes
4. Neighborhood analysis, spatially variable genes, and cell-cell communication

## Repository structure

```
notebook/
  visium_hd_quickstart.ipynb          # Quickstart (run on Colab)
scripts/
  utils.py                            # Helper functions used by the notebook
  generate_precomputed.py             # Script used to create the precomputed data: no need to touch this
intro_to_spatial_background/
  intro_to_spatial_background.md      # Background reading on spatial transcriptomics
```
