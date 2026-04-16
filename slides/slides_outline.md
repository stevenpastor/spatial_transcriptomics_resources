# Visium HD Spatial Transcriptomics: Slide Outline (4 slides)

## Slide 1: Why Spatial Transcriptomics

scRNA-seq dissociates tissue, so you lose where cells sit. Spatial transcriptomics keeps the coordinates. In a tumor, immune cell position relative to tumor cells predicts therapy response, and that signal is gone after dissociation.

![Visium HD slides](https://cdn.10xgenomics.com/image/upload/v1713545308/blog/Visium%20HD%20NPI/Visium%20HD%20slides.png)

### Two technology families

| Feature        | Sequencing-based (Visium HD, Slide-seq) | Imaging-based (MERFISH, Xenium, CosMx) |
|----------------|----------------------------------------|----------------------------------------|
| Gene coverage  | Whole transcriptome (~18k+ genes)      | Panels (100 to 6,000 genes)            |
| Resolution     | 2 to 55 um                             | Subcellular (~100 nm)                  |
| Tissue         | FFPE and fresh frozen                  | Mostly fresh frozen                    |

---

## Slide 2: Visium HD Technology

Continuous 2 um barcoded oligos over a 6.5 mm square, transferred onto tissue by the CytAssist instrument. Probe-based chemistry, works on FFPE and fresh frozen.

### Binning

- **2 um**: max resolution, very sparse, rarely used directly
- **8 um**: ~1 to 3 cells per bin, primary analysis resolution
- **16 um**: ~5 to 15 cells, good for QC and exploration
- **Segmented (v4.0+)**: StarDist nuclei detection gives one row per cell

### Visium v1 vs. Visium HD

| Feature        | Visium v1              | Visium HD             |
|----------------|------------------------|-----------------------|
| Resolution     | 55 um spots            | 2 um continuous       |
| Cells / unit   | ~5 to 30 per spot      | ~1 to 3 per 8 um bin  |
| Gaps           | 100 um between spots   | None                  |
| FFPE           | Yes (v2)               | Yes (primary)         |
| Segmentation   | No                     | Yes (StarDist)        |

Space Ranger processes raw data into binned matrices + spatial coordinates in standard 10x format.

---

## Slide 3: When to Use Visium HD

### Good fit

- Whole-transcriptome discovery with no prior gene list
- FFPE archives and clinical biobanks
- Near-single-cell resolution (8 um is ~1 to 3 cells)
- Heterogeneous tissue with multiple compartments

### Not a good fit

- You already have a gene panel (Xenium/CosMx are faster and cheaper)
- You need subcellular localization (use MERFISH or seqFISH+)
- Live-cell dynamics (fixed tissue only)
- Large cohorts above 50 samples (cost and compute add up)

---

## Slide 4: Analysis Pipeline

1. QC and filtering (spatial-aware, not just global thresholds)
2. Normalization and highly variable gene selection
3. PCA and UMAP
4. Cell-type annotation (marker scoring, deconvolution, or label transfer)
5. Spatial analysis: neighborhood enrichment, spatially variable genes
6. Cell-cell communication (ligand-receptor)
7. Validation against ground truth or orthogonal data

Spatial QC catches artifacts that global QC misses: tissue folds, air bubbles, edge effects. Always check QC metrics on the tissue coordinates, not just in histograms.
