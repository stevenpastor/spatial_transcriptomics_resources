# Visium HD Spatial Transcriptomics: Slide Outline (4 slides)

## Slide 1: Why Spatial Transcriptomics

scRNA-seq dissociates tissue, so you lose where cells sit. Spatial transcriptomics keeps the coordinates. In a tumor, immune cell position relative to tumor cells predicts therapy response, and that signal is gone after dissociation.

### Two technology families

| Feature        | Sequencing-based (Visium HD, Slide-seq) | Imaging-based (MERFISH, Xenium, CosMx) |
|----------------|----------------------------------------|----------------------------------------|
| Gene coverage  | Whole transcriptome                    | Panels (100 to 6,000 genes)            |
| Resolution     | 2 to 55 um                             | Subcellular (~100 nm)                  |
| Tissue         | FFPE and fresh frozen                  | Mostly fresh frozen                    |

---

## Slide 2: How Visium HD Works

### The slide: a continuous oligo lawn

The Visium HD slide contains a 6.5 x 6.5 mm capture area covered in millions of 2 x 2 um barcoded oligonucleotide squares with no gaps. Each square has a unique spatial barcode, so any molecule captured at that position is tagged with its location.

![Visium HD slides](https://cdn.10xgenomics.com/image/upload/v1713545308/blog/Visium%20HD%20NPI/Visium%20HD%20slides.png)

### CytAssist: linking tissue to oligos

Tissue is sectioned and stained (H&E) on a standard glass slide, then placed into the CytAssist instrument alongside the Visium HD slide. CytAssist transfers hybridized probes from the tissue section onto the oligo lawn, where they bind to the spatially barcoded capture probes. This links each transcript to its x,y coordinate on the tissue.

![Slides to CytAssist](https://cdn.10xgenomics.com/image/upload/v1713545305/blog/Visium%20HD%20NPI/Slides%20to%20CytAssist.png)

### From tissue to sequencing data

![CytAssist Workflow](https://cdn.10xgenomics.com/image/upload/v1678748956/products/CytAssist/CytAssist_Workflow.png)

1. **Sample prep and imaging**: section tissue onto glass slide, stain with H&E, image at high resolution
2. **CytAssist transfer**: probes hybridize to mRNA in the tissue, then CytAssist transfers them onto the Visium HD capture slide
3. **Probe extension**: captured probes are extended, incorporating the spatial barcode and UMI
4. **Library construction**: standard Illumina-compatible library is built from the barcoded cDNA
5. **Sequencing**: each read carries the spatial barcode (location), UMI (molecule identity), and the gene sequence

### Binning and segmentation

Raw data is at 2 um resolution but very sparse. Space Ranger aggregates into larger bins:

- **2 um**: max resolution, very sparse, rarely used directly
- **8 um**: ~1 to 3 cells per bin, primary analysis resolution
- **16 um**: ~5 to 15 cells, good for QC and exploration

Even at 8 um, each bin may contain parts of multiple cells, so gene expression is a mixture. **Segmentation** (Space Ranger v4.0+) solves this by running StarDist nuclei detection on the H&E image to draw a polygon around each nucleus, then assigning barcodes to individual cells. The result is one row per cell instead of one row per bin, giving true single-cell resolution. This matters for cell type annotation and any analysis where mixing signal from neighboring cells would blur the result.

---

## Slide 3: Space Ranger - From Raw Data to Spatial Matrices

### Inputs

Space Ranger needs three inputs to produce spatial gene expression data:

| Input | What it is | Why it is needed |
|-------|-----------|-----------------|
| **FASTQ files** | Raw sequencing reads from Illumina | Contain the spatial barcode, UMI, and gene sequence for every captured molecule |
| **H&E image** | High-resolution microscope image of the stained tissue section | Provides the tissue morphology context; bins are mapped back to pixel coordinates on this image for all downstream spatial plots |
| **CytAssist image** | Image taken by the CytAssist instrument during probe transfer | Used to align the tissue position on the capture slide; connects where the tissue sat on the oligo lawn to the barcode grid |

Additional required inputs: reference transcriptome, probe set CSV, and slide/area IDs.

### Image registration: Loupe Browser

Before running Space Ranger, the H&E image and CytAssist image need to be aligned so that barcodes map to the correct tissue coordinates. Space Ranger can auto-detect this alignment, but for difficult cases (rotated sections, low contrast) you can manually align in Loupe Browser:

1. Upload the CytAssist image and identify the fiducial markers (the circular dots around the capture area frame)
2. Upload the H&E microscope image and pin matching landmarks on both images
3. Assess alignment quality and optionally auto-refine

![Loupe Browser: Identify fiducials](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/06_identify_fiducials.png)

![Loupe Browser: Assess alignment](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/18_assess_alignment.png)

The exported alignment file (JSON) is passed to Space Ranger via `--loupe-alignment` so the pipeline uses your manual registration instead of auto-detection.

### Outputs

```
CRC_P1/outs/
├── binned_outputs/
│   ├── square_002um/                     # 2 um bins (highest res, sparsest)
│   ├── square_008um/                     # 8 um bins (primary analysis)
│   │   ├── filtered_feature_bc_matrix/   # gene x bin count matrix (MEX format)
│   │   ├── filtered_feature_bc_matrix.h5 # same matrix in HDF5
│   │   └── spatial/
│   │       ├── tissue_positions.parquet  # bin coordinates on the H&E image
│   │       ├── scalefactors_json.json    # pixel-to-um conversion factors
│   │       ├── tissue_hires_image.png    # downscaled H&E
│   │       └── tissue_lowres_image.png   # further downscaled H&E
│   └── square_016um/                     # 16 um bins (exploration/QC)
├── segmented_outputs/                    # v4.0+: one row per segmented cell
│   └── filtered_feature_cell_matrix.h5
├── cloupe_008um.cloupe                   # interactive viewer file
├── web_summary.html                      # QC report
└── metrics_summary.csv
```

The key spatial file is `tissue_positions.parquet`, which maps every barcode to its pixel coordinates on the H&E image. This is what puts the "spatial" in spatial transcriptomics: without it, the count matrix is just a standard scRNA-seq-like matrix with no location information.

### What we load in the tutorial

In the tutorial notebooks, we load pre-processed data from these Space Ranger outputs:

- **`adata_8um.h5ad`**: the 8 um binned count matrix (`filtered_feature_bc_matrix.h5`) already converted to an AnnData object with spatial coordinates from `tissue_positions.parquet` stored in `adata.obsm["spatial"]`. This is the primary object for all analysis.
- **`adata_8um_annotated.h5ad`**: same object after normalization, PCA, UMAP, Leiden clustering, and cell type annotation have been run. Pre-computed so the Colab notebook runs in minutes instead of hours.
- **`ground_truth_labels.csv`**: published cell type labels from the original study (de Oliveira et al., Nature Genetics 2025) for validation.
- **H&E image**: the tissue image for spatial overlay plots.

We use the 8 um resolution because it balances signal density (~1 to 3 cells per bin) with manageable data size.

---

## Slide 4: Analysis Pipeline

### 1. Quality control and filtering

**Why**: raw spatial data contains empty bins, damaged tissue regions, and technical artifacts. Bins under tissue folds or air bubbles have inflated or depleted counts that do not reflect real biology. Filtering these out prevents them from distorting normalization, clustering, and every step downstream.

- Per-bin metrics: total counts, number of genes detected, mitochondrial %
- Spatial QC: plot these metrics on tissue coordinates to catch artifacts that look normal in histograms but are spatially clustered (folds, edge effects, bubbles)
- Spatial outlier detection: flag bins whose metrics deviate from their local neighbors

### 2. Normalization and feature selection

**Why**: raw counts are confounded by sequencing depth. A bin with 2x more total counts will appear to express every gene at 2x the level. Normalization (e.g., total-count scaling + log transform) removes this technical variation so that differences between bins reflect biology, not capture efficiency. Selecting highly variable genes focuses downstream analysis on genes that actually vary across the tissue and reduces noise from housekeeping genes.

### 3. Dimensionality reduction (PCA and UMAP)

**Why**: with thousands of genes per bin, you cannot cluster or visualize directly. PCA compresses the data into 30 to 50 principal components that capture the major axes of variation. UMAP projects those components into 2D for visualization, placing bins with similar expression profiles near each other. This lets you see whether the data separates into distinct cell populations before assigning labels.

### 4. Cell type annotation

**Why**: clusters from step 3 are just numbered groups. Annotation gives them biological meaning. Multiple approaches:

- **Marker gene scoring**: score each bin for known marker gene sets (e.g., EPCAM/KRT for tumor, CD3D/CD3E for T cells). Fast, interpretable, no reference dataset needed.
- **Label transfer**: project bins onto an annotated scRNA-seq reference and transfer labels. More automated, but depends on reference quality.
- **Deconvolution**: estimate the fraction of each cell type within a bin. Useful at coarser resolutions where bins contain multiple cells.

### 5. Spatial neighborhood analysis

**Why**: knowing cell types is not enough. The same T cell has different functional implications depending on whether it sits next to a tumor cell or in a lymphoid aggregate. Spatial analysis asks: which cell types co-localize more (or less) than expected by chance?

- **Neighborhood enrichment**: for each cell type pair, count how often they are neighbors vs. a random permutation baseline. Produces a z-score matrix showing which pairs are enriched or depleted.
- **Co-occurrence analysis**: measures how the probability of two cell types co-occurring changes as a function of distance. Captures gradients that neighborhood enrichment (which uses a fixed neighbor definition) can miss.

### 6. Spatially variable genes

**Why**: standard differential expression finds genes that differ between groups (e.g., tumor vs. stroma) but ignores spatial pattern. A gene can be spatially structured, expressed in a gradient or confined to a boundary, without being differentially expressed between predefined groups. Moran's I measures spatial autocorrelation: whether nearby bins have more similar expression than expected. High Moran's I genes reveal spatial organization that group-based DE would miss entirely.

### 7. Cell-cell communication (ligand-receptor)

**Why**: cells communicate through ligand-receptor interactions, and spatial transcriptomics lets you test whether a ligand and its receptor are actually expressed in neighboring cell types, not just co-expressed somewhere in the dataset. Squidpy's ligrec analysis tests each ligand-receptor pair across each cell type pair, using permutation tests to assess significance. This identifies specific signaling axes between spatially adjacent populations (e.g., tumor-secreted ligands binding receptors on nearby immune cells).

### 8. Validation

**Why**: computational annotations and spatial patterns need independent confirmation. Comparing marker-based annotations against published ground truth labels checks whether your pipeline recovers known biology. Examining marker gene expression across annotated clusters (e.g., dotplots) confirms that annotations are biologically coherent and not artifacts of the clustering or scoring approach.
