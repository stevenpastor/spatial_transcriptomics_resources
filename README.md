# Spatial Transcriptomics Resources

A collection of tutorials, tools, and references for spatial transcriptomics analysis. Under construction. In the future, this repository will include links and resources for all things spatial to help users determine how best to analyze their data.

---

## Start here

New to spatial transcriptomics? Read [`intro_to_spatial_background.md`](intro_to_spatial_background/intro_to_spatial_background.md) first for a primer on why spatial data matters, how Visium HD works, and the analysis vocabulary used downstream. Currently the focus is Visium HD and I will expand to other techs in the near future. Once you understand the background, move on to the quickstart tutorial below.

---

## Tutorials

### Visium HD Spatial Transcriptomics

#### Quick Start

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/stevenpastor/spatial_transcriptomics_resources/blob/main/notebook/visium_hd_quickstart.ipynb)

Analyze a 10x Genomics Visium HD spatial transcriptomics dataset with QC, filtering, cell type annotation, neighborhood analysis, spatially variable genes, and ligand-receptor communication.

Click the **Colab** badge above to run the quickstart notebook. It auto-installs dependencies and downloads pre-processed data (~150 MB) from Figshare. No local setup needed.

##### Dataset

Human Colorectal Cancer (CRC), Patient 1 from the [Nature Genetics 2025 publication](https://www.nature.com/articles/s41588-025-02193-3). Pre-processed to ~50K bins at 8 um resolution with ground-truth cell type annotations from the original study group.

#### What you will learn

1. Load and explore Visium HD spatial transcriptomics data
2. Quality control with filtering (per-bin metrics + spatial outlier detection)
3. Compare marker-based annotations against published ground-truth labels
4. Spatial neighborhood analysis, spatially variable genes, cell-cell communication

## Repository structure

```
notebook/
  visium_hd_quickstart.ipynb          # Quickstart (run on Colab)
scripts/
  utils.py                            # Helper functions used by the notebook
  generate_precomputed.py             # Script used to create the precomputed data
intro_to_spatial_background/
  intro_to_spatial_background.md      # Background reading on spatial transcriptomics
```
