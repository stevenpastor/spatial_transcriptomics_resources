# Visium HD Spatial Transcriptomics

## 1. Why Spatial Transcriptomics

### The problem with dissociation

Single-cell and single-nucleus RNA-seq (scRNA-seq, snRNA-seq) measure what each cell expresses, but the tissue has to be dissociated first. Once the cells come apart from the section, the spatial picture is gone. You no longer know where a given cell sat, what its neighbors were, or which tissue compartment it belonged to. For a lot of biological questions, that gap is the whole story.

Cell function depends on context. A fibroblast at a wound margin behaves nothing like one buried in deep stroma. Immune cells touching tumor cells run different programs than the ones sitting in a distant lymphoid aggregate. Signaling gradients, physical cell-cell contacts, and local oxygen all shape gene expression, and every one of them is a spatial property. Once you have dissociated the tissue, you cannot recover any of that after the fact.

### What spatial transcriptomics adds

Spatial transcriptomics measures gene expression on intact tissue sections, so you keep the expression profile and the x,y coordinate of every measurement. That opens up questions you cannot ask with scRNA-seq:

- Which cell types live next to each other, and is that co-localization more frequent than chance would predict?
- Do certain genes form spatial gradients, sit along boundaries, or stay confined to specific tissue regions?
- For a candidate ligand-receptor pair, are the two cells actually close enough to talk?
- How does the local microenvironment shift across a tumor-stroma boundary or along the length of an organ?

### Two technology families

Spatial methods fall into two camps. The sequencing-based methods capture mRNA on a patterned slide of barcoded oligonucleotides and then sequence the recovered cDNA, giving you the whole transcriptome at whatever spatial resolution the slide pattern supports. Imaging-based methods take a different route: they detect a predefined panel of genes directly in the tissue with rounds of fluorescent probe hybridization and imaging. The trade is precision for breadth. You get subcellular resolution on the imaging side, but only on the genes you chose ahead of time.

| Method            | Family          | Resolution        | Gene coverage          | Tissue type         | Typical use case                                              |
|-------------------|-----------------|-------------------|------------------------|---------------------|---------------------------------------------------------------|
| Visium (standard) | Sequencing      | 55 um spots       | Whole transcriptome    | FFPE, fresh frozen  | Pilot studies, regional expression, lower per-sample cost     |
| Visium HD         | Sequencing      | 2 um bins         | Whole transcriptome    | FFPE, fresh frozen  | Near single-cell resolution on archival samples               |
| Slide-seqV2       | Sequencing      | 10 um beads       | Whole transcriptome    | Fresh frozen        | Higher resolution sequencing on fresh tissue                  |
| Stereo-seq        | Sequencing      | 220 nm to ~1 um   | Whole transcriptome    | Fresh frozen        | Large tissue areas (whole organs, embryos) at high resolution |
| Xenium            | Imaging         | Subcellular       | 100 to ~5,000 genes    | FFPE, fresh frozen  | Validated targeted panels with cell boundaries                |
| MERFISH/MERSCOPE  | Imaging         | Subcellular       | 100 to ~1,000 genes    | Mostly fresh frozen | Custom panels for circuit-level or developmental questions    |
| CosMx SMI         | Imaging         | Subcellular       | 1,000 to 6,000 genes   | FFPE, fresh frozen  | Larger gene panels including proteins on the same section     |

### How to choose

A few rules of thumb. For whole-transcriptome discovery on archival FFPE blocks, reach for Visium HD first. If you have a focused gene list and want true subcellular resolution with cell boundaries, the choice between Xenium, MERFISH/MERSCOPE, and CosMx comes down to which panel fits your biology. Standard Visium still pulls its weight for pilot work, or when you care more about sample numbers than within-section resolution. Stereo-seq and Slide-seqV2 sit in their own niche: high-resolution whole-transcriptome work on fresh frozen material.

This tutorial leans on Visium HD because it handles FFPE, captures the full transcriptome, and ships with a polished pipeline (Space Ranger) and viewer (Loupe Browser) that most labs already have on hand. The analysis ideas covered here (QC, clustering, neighborhood analysis, spatially variable genes, ligand-receptor) port over to the other platforms with only small adjustments.

---

## 2. How Visium HD Works

### The slide: a continuous oligo lawn

The Visium HD slide has a 6.5 x 6.5 mm capture area tiled edge-to-edge with millions of 2 x 2 um barcoded oligonucleotide squares. There are no gaps between them. Every square carries its own spatial barcode, so anything captured at that position picks up the barcode of that location.

![Visium HD slides](https://cdn.10xgenomics.com/image/upload/v1713545308/blog/Visium%20HD%20NPI/Visium%20HD%20slides.png)

### CytAssist: linking tissue to oligos

Tissue gets sectioned and H&E-stained on a standard glass slide, then loaded into the CytAssist instrument alongside the Visium HD slide. CytAssist transfers hybridized probes from the tissue section onto the oligo lawn, where they bind to the spatially barcoded capture probes. That handoff is what ties each transcript to an x,y coordinate on the tissue.

![Slides to CytAssist](https://cdn.10xgenomics.com/image/upload/v1713545305/blog/Visium%20HD%20NPI/Slides%20to%20CytAssist.png)

### From tissue to sequencing data

![CytAssist Workflow](https://cdn.10xgenomics.com/image/upload/v1678748956/products/CytAssist/CytAssist_Workflow.png)

1. **Sample prep and imaging**: section the tissue onto a glass slide, stain with H&E, and image at high resolution.
2. **CytAssist transfer**: probes hybridize to mRNA in the tissue, then CytAssist moves them onto the Visium HD capture slide.
3. **Probe extension**: captured probes get extended, picking up the spatial barcode and UMI in the process.
4. **Library construction**: a standard Illumina-compatible library is built from the barcoded cDNA.
5. **Sequencing**: every read carries three things: the spatial barcode (where), the UMI (which molecule), and the gene sequence (what).

### Binning and segmentation

Raw data sits at 2 um resolution, but at that scale it is also very sparse. Space Ranger aggregates the squares into larger bins:

- **2 um**: max resolution, very sparse, almost never used directly
- **8 um**: ~1 to 3 cells per bin; the resolution most analyses run at
- **16 um**: ~5 to 15 cells; handy for QC and quick exploration

Even at 8 um, a single bin can clip parts of more than one cell, so the gene expression you see is a mixture. Segmentation (Space Ranger v4.0+) handles this by running StarDist nuclei detection on the H&E image, drawing a polygon around each nucleus, and assigning the underlying barcodes to individual cells. You end up with one row per cell rather than one row per bin, which is closer to true single-cell resolution. That matters for cell type annotation, and for anything where signal from neighboring cells would smear the result. The segmented output is its own deliverable from Space Ranger and comes with a few caveats worth knowing before you commit to it, covered in the "Segmentation output" subsection below.

---

## 3. Space Ranger: From Raw Data to Spatial Matrices

### Inputs

Space Ranger needs three inputs to produce spatial gene expression data:

| Input | What it is | Why it is needed |
|-------|-----------|-----------------|
| **FASTQ files** | Raw sequencing reads from Illumina | Contain the spatial barcode, UMI, and gene sequence for every captured molecule |
| **H&E image** | High-resolution microscope image of the stained tissue section | Provides the tissue morphology context; bins are mapped back to pixel coordinates on this image for all downstream spatial plots |
| **CytAssist image** | Image taken by the CytAssist instrument during probe transfer | Used to align the tissue position on the capture slide; connects where the tissue sat on the oligo lawn to the barcode grid |

Three more inputs are required, and they are easy to get wrong:

- **Reference transcriptome**: a prebuilt 10x reference for your species (human GRCh38 or mouse GRCm39). Grab it from the 10x Genomics support site rather than building your own, unless you have a real reason to. The version should line up with the probe set release notes.
- **Probe set CSV**: the list of probe sequences for your assay. Visium HD uses its own probe set, distinct from standard Visium, so the CSV has to match the chemistry. 10x ships it with the kit and also posts it for download. If the wrong CSV gets passed in, Space Ranger will happily run and hand back nonsense counts.
- **Slide serial number and capture area ID**: both printed on the slide itself. Space Ranger uses them to look up the correct fiducial layout. Mistype the serial and registration falls over.

### Image registration: Loupe Browser

Before Space Ranger can do anything useful, the H&E image and CytAssist image need to be aligned so barcodes line up with the right tissue coordinates. Space Ranger can auto-detect that alignment in most cases, and the right move is to let it try first. Manual alignment in Loupe Browser is the fallback for when auto-detection misfires.

**When do you need manual alignment?** A handful of situations push you off the auto path. The tissue section may have been rotated, or placed off-center on the slide. Faint or uneven H&E staining can throw the detector off. Debris and air bubbles sometimes obscure the fiducial frame. And Space Ranger's web summary may flat-out report low alignment confidence after a first auto-attempt. Spatial plots that look smeared, shifted, or rotated relative to the H&E are another tell that the alignment needs another look.

The Loupe Browser workflow runs in three steps:

1. **Upload the CytAssist image and identify the fiducial markers.** Fiducials are the small circular dots arranged in a frame around the capture area, and they are the ground truth for where the barcode grid sits. Loupe will try to find them automatically. Any that land on tissue or debris instead of on the dots need a manual fix.
2. **Upload the H&E microscope image and pin matching landmarks on both images.** Good landmark pairs are features visible in both images, like a vessel cross-section, a duct, a fold edge, or the corner of a tissue region. Aim for at least four well-spread pairs across the capture area. Landmarks all clustered in one spot underdetermine the transform.
3. **Assess alignment quality and optionally auto-refine.** Loupe shows the overlay with a quality readout. You want H&E features sitting cleanly on top of their CytAssist counterparts. If one corner of the section is off, drop a landmark there and re-refine.

![Loupe Browser: Identify fiducials](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/06_identify_fiducials.png)

![Loupe Browser: Assess alignment](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/18_assess_alignment.png)

The exported alignment file is a small JSON. Pass it to Space Ranger with `--loupe-alignment=path/to/alignment.json` and the pipeline uses your manual registration instead of auto-detection. The JSON is portable and worth committing to a project repo so the registration step stays reproducible. Space Ranger refuses to run if the JSON's slide serial does not match the FASTQs or the CytAssist image, so a mismatched file fails fast instead of silently producing wrong coordinates.

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

The key spatial file is `tissue_positions.parquet`, which maps every barcode to its pixel coordinates on the H&E image. That file is what puts the "spatial" in spatial transcriptomics. Without it, the count matrix is just an scRNA-seq-shaped matrix with no location information.

### Which outputs to actually use

Space Ranger drops a lot of files. Most analyses only touch a handful. Here are the ones to grab, and why:

| File | What it is | When to use it |
|------|-----------|----------------|
| `binned_outputs/square_008um/filtered_feature_bc_matrix.h5` | Gene-by-bin count matrix in HDF5 | Default load target. This is what `scanpy.read_10x_h5` opens into AnnData. |
| `binned_outputs/square_008um/spatial/tissue_positions.parquet` | Per-bin pixel coordinates | Join to AnnData on barcode to populate `adata.obsm["spatial"]`. |
| `binned_outputs/square_008um/spatial/scalefactors_json.json` | Scale factors mapping bin coords to full-res image pixels | Needed when overlaying spots on the H&E image, especially when using the lowres/hires downsamples. |
| `binned_outputs/square_008um/spatial/tissue_hires_image.png` | Downscaled H&E ready to plot | The background image for `sc.pl.spatial` and similar overlays. |
| `segmented_outputs/filtered_feature_cell_matrix.h5` | Gene-by-cell matrix from nuclei segmentation | Use when you want single-cell resolution (see "Segmentation output" below). |
| `cloupe_008um.cloupe` | Interactive viewer file | Open in Loupe Browser for exploration; not for Python analysis. |
| `web_summary.html` | QC report (alignment, saturation, reads in tissue, fiducial detection) | Read this before doing anything else. If alignment or saturation are red, your downstream conclusions are at risk. |
| `metrics_summary.csv` | Same metrics in machine-readable form | Use when you need to track QC across many samples programmatically. |

A few outputs are usually safe to skip. The `square_002um/` matrix is the raw resolution, but it is so sparse that most analyses just aggregate it back up to 8 um (or rely on segmentation). The `raw_feature_bc_matrix` (no `filtered_` prefix) includes bins outside the tissue boundary, so reach for it only when you need to debug the in-tissue filter. The `square_016um/` matrix is fine for quick spatial QC plots but throws away too much resolution for serious cell typing.

### Segmentation output

Segmentation hands you a cell-level matrix in place of a bin-level matrix. Space Ranger v4.0+ runs StarDist on the H&E image to find nuclei, draws a polygon around each one, and assigns the underlying 2 um barcodes to whichever nucleus polygon contains them (with some local rules for cytoplasmic capture). One row per nucleus is the result.

Inside `segmented_outputs/` you get:

- `filtered_feature_cell_matrix.h5`: the cell-level count matrix.
- `nucleus_segmentations.geojson`: the nucleus polygons. Useful for plotting cell outlines and for computing centroids when you want spatial coordinates in your AnnData. `utils.load_visium_hd_segmented` reads this file and treats the polygon centroids as `obsm["spatial"]`.
- A per-cell metrics CSV with counts, gene counts, and nucleus area.

A couple of trade-offs worth knowing up front. Segmentation produces cleaner per-cell expression and matches scRNA-seq conventions, so cell type annotation tools tend to behave better on the segmented output. The cost is that any nucleus the segmenter misses loses its transcripts (or hands them off to a neighbor), and the totals you measure are not absolute, since they depend on the segmenter as much as on the tissue. Lean on segmentation when cell identity is what you care about. Stay on bins when you want quantitative tissue-level expression patterns.

When the bundled StarDist segmentation does not work well on your tissue, a few alternatives are available:

- **bin2cell** (Python): takes the Space Ranger 2 um output plus an external nucleus mask and aggregates bins into cells. This is the standard escape hatch for swapping in a different segmenter.
- **StarDist (general)**: a nucleus segmentation model. Run it on the H&E image to produce a nucleus mask, then feed that mask to bin2cell.
- **Cellpose**: a general-purpose cell and nucleus segmentation model. People often use it for whole-cell masks from H&E or membrane stains, and the resulting mask plugs into bin2cell the same way.

The choice comes down to your tissue. StarDist is the default for H&E nuclei. Cellpose tends to handle harder cases and other stains better.

Once you have a cell-level AnnData (either from `segmented_outputs/` or from bin2cell), downstream analysis looks the same as scRNA-seq: normalize, find highly variable genes, PCA, neighbors, UMAP, Leiden, annotate. The spatial layer comes back through `squidpy.gr.spatial_neighbors`, only now it operates on cell centroids instead of bin centroids. Every analysis covered in Section 4 still applies.

### What you will learn in the quickstart

The quickstart notebook (`notebook/visium_hd_quickstart.ipynb`) walks through:

1. Loading and poking around Visium HD data
2. Quality control with filtering, both per-bin metrics and spatial outlier detection
3. Marker-based annotation, compared against the published ground-truth labels
4. Spatial neighborhood analysis, spatially variable genes, and cell-cell communication

The notebook loads pre-processed data instead of running Space Ranger end-to-end, so Colab runtimes stay short. The files it pulls from Figshare are:

- **`adata_8um.h5ad`**: the 8 um binned count matrix (`filtered_feature_bc_matrix.h5`) already converted to an AnnData object, with spatial coordinates from `tissue_positions.parquet` stored in `adata.obsm["spatial"]`. This is the main object the analysis works against.
- **`adata_8um_annotated.h5ad`**: the same object after normalization, PCA, UMAP, Leiden clustering, and cell type annotation. Pre-computed so the Colab notebook finishes in minutes rather than hours.
- **`ground_truth_labels.csv`**: published cell type labels from the original study (de Oliveira et al., Nature Genetics 2025) for validation.
- **H&E image**: the tissue image used for spatial overlay plots.

We work at 8 um resolution because it balances signal density (~1 to 3 cells per bin) against manageable data size.

---

## 4. Python Data Structures and Tooling

Before opening the notebook, it helps to know the three libraries that show up in every cell: `anndata`, `scanpy`, and `squidpy`. The R equivalents (Seurat, SpatialExperiment) will get their own primer in a follow-up doc. Everything here is Python.

### AnnData: one object, many slots

scRNA-seq and spatial datasets need the same things bundled together: a count matrix, per-cell metadata, per-gene metadata, reduced-dimension embeddings, and a place to dump arbitrary side data like images, neighbor graphs, and parameters. AnnData is the standard container that holds all of that in one object with aligned indices.

The slots you will touch:

```python
adata.X              # cells x genes count matrix (usually sparse)
adata.obs            # per-cell (or per-bin) metadata DataFrame
adata.var            # per-gene metadata DataFrame
adata.obsm["spatial"]   # cells x 2 array of x,y pixel coordinates
adata.obsm["X_pca"]     # PCA embedding
adata.obsm["X_umap"]    # UMAP embedding
adata.uns            # unstructured dict (e.g., spatial image and scale factors)
adata.layers         # alternative matrices, e.g. raw counts kept after normalization
adata.obsp           # per-cell pairwise matrices (e.g., neighbor graphs)
```

Indexing works like a DataFrame. `adata[adata.obs["cell_type"] == "T cell"]` returns a view of only those cells, and `adata[:, "EPCAM"]` selects a single gene.

A spatial AnnData and an scRNA-seq AnnData share the same structure. The differences are pretty small: the spatial one fills in `obsm["spatial"]` with x,y coordinates, and usually stashes the tissue image and scale factors under `uns["spatial"]`. Every other operation (normalization, clustering, DE, embeddings) behaves identically on either kind of object. That is the whole point of the format. Tools written for one work on the other.

Reading and writing is `sc.read_h5ad("file.h5ad")` and `adata.write_h5ad("file.h5ad")`. The tutorial data ships as `.h5ad` files, which is HDF5 with an AnnData schema on top.

### scanpy: the analysis workhorse

`scanpy` runs the standard pipeline (QC, normalization, highly variable gene selection, PCA, neighbors, UMAP, Leiden, differential expression, plotting) directly on AnnData. The canonical flow looks like this:

```python
import scanpy as sc

sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
```

For plotting, `sc.pl.umap(adata, color="leiden")` draws the UMAP scatter, `sc.pl.violin` handles per-cluster expression plots, and `sc.pl.spatial(adata, color="cell_type")` builds the spatial overlay from `obsm["spatial"]` and the H&E in `uns["spatial"]`.

### squidpy: the spatial layer on top

`squidpy` bolts spatial-aware analysis onto the scanpy/AnnData stack: neighbor graphs built from coordinates, neighborhood enrichment, spatial autocorrelation (Moran's I), and image-based features. The foundation is the spatial neighbor graph:

```python
import squidpy as sq

sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=8)
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
sq.gr.spatial_autocorr(adata, mode="moran")
```

After `spatial_neighbors` runs, the connectivity graph lives in `adata.obsp["spatial_connectivities"]` and every downstream squidpy call uses it. This is the spatial analog of `sc.pp.neighbors`, which builds an expression-space neighbor graph in `adata.obsp["connectivities"]`. Both can sit on the same AnnData without issue.

### Adjacent libraries worth knowing

- `anndata`: the library that defines the AnnData class. You rarely import it directly. Scanpy re-exports what you need.
- `bin2cell`: bin-to-cell aggregation for Visium HD, useful when you want to swap in a custom nucleus segmenter.
- `cell2location`, `tangram`, `RCTD` (Python ports): deconvolution methods that estimate cell type fractions per bin when resolution drops below single-cell.
- `spatialdata` and `spatialdata-io`: an emerging unified format for multi-modal spatial datasets (images, points, shapes, tables). Useful once you start mixing modalities. AnnData is still the practical default for plain Visium HD analysis.

---

## 5. Analysis Pipeline

### Before analysis: upstream processing

By the time we open a notebook, a handful of steps are already in the rearview:

1. **Tissue preparation and imaging**: the tissue was sectioned onto a glass slide, stained with H&E, and imaged at high resolution. That H&E image plays double duty as the morphology reference and the coordinate system for all spatial plots.

2. **CytAssist and sequencing**: the section ran through CytAssist to transfer probes onto the Visium HD capture slide, libraries were built, and the sample was sequenced on an Illumina instrument. What lands on disk is FASTQs.

3. **Image registration**: the H&E microscope image and the CytAssist image need to be lined up so barcodes map to the correct tissue positions. Space Ranger can auto-detect that alignment in most cases. For tricky samples (rotated sections, low contrast, or poor fiducial detection), Loupe Browser offers a manual workflow where you identify fiducials and pin landmarks between the two images. If you did the alignment manually, the exported JSON file gets passed to Space Ranger via `--loupe-alignment`.

4. **Space Ranger**: takes the FASTQs, H&E image, CytAssist image (plus the alignment file if applicable), reference transcriptome, and probe set as inputs. It demultiplexes reads, maps them to the genome, registers images, and writes out binned count matrices at multiple resolutions (2, 8, 16 um), spatial coordinate files, and a QC report. This is the last "black box" step before interactive analysis begins.

Everything past this point is what we cover in the tutorial, starting from the Space Ranger outputs loaded as AnnData objects (see Section 4 for the data structures and libraries involved).

---

### 1. Quality control and filtering

Raw spatial data is full of empty bins, damaged tissue regions, and technical junk. Bins sitting under tissue folds or air bubbles end up with inflated or depleted counts that have nothing to do with biology. Filter them out early or they will warp normalization, clustering, and everything that follows.

- Per-bin metrics: total counts, number of genes detected, mitochondrial %
- Spatial QC: plot those same metrics on tissue coordinates to catch artifacts that look fine in a histogram but cluster spatially (folds, edge effects, bubbles)
- Spatial outlier detection: flag bins whose metrics drift away from their local neighbors

### 2. Normalization and feature selection

Sequencing depth is a confound, plain and simple. A bin with 2x more total counts will look like it expresses every gene at 2x the level. Normalization (total-count scaling plus a log transform, for example) strips that technical variation out so the differences you see between bins come from biology, not capture efficiency. Picking highly variable genes then narrows downstream analysis to the genes that actually vary across the tissue, and cuts the noise from housekeeping genes.

### 3. Dimensionality reduction (PCA and UMAP)

Thousands of genes per bin is too much to cluster or visualize directly. PCA shrinks the data into 30 to 50 principal components that capture the major axes of variation. UMAP then projects those components into 2D, placing bins with similar expression profiles near each other. You get a chance to see whether the data already separates into distinct populations before you start labeling anything.

### 4. Cell type annotation

The clusters from step 3 are just numbered groups. Annotation is where they pick up biological meaning. There are a few common approaches:

- **Marker gene scoring**: score each bin against known marker gene sets (EPCAM/KRT for tumor, CD3D/CD3E for T cells, and so on). Fast, interpretable, no reference dataset needed.
- **Label transfer**: project bins onto an annotated scRNA-seq reference and transfer the labels across. More automated, but only as good as the reference.
- **Deconvolution**: estimate the fraction of each cell type within a bin. Most useful at coarser resolutions where bins contain multiple cells.

### 5. Spatial neighborhood analysis

Knowing cell types is not enough on its own. The same T cell means different things biologically depending on whether it sits next to a tumor cell or buried inside a lymphoid aggregate. Spatial analysis asks the next question: which cell types co-localize more (or less) than chance?

- **Neighborhood enrichment**: for each cell type pair, count how often they show up as neighbors and compare that to a random permutation baseline. The output is a z-score matrix showing which pairs are enriched or depleted.
- **Co-occurrence analysis**: measures how the probability of two cell types co-occurring changes with distance. Catches gradients that neighborhood enrichment, with its fixed neighbor definition, would miss.

### 6. Spatially variable genes

Standard differential expression finds genes that differ between groups (tumor vs. stroma, say), but it has no sense of spatial pattern. A gene can be spatially structured (a gradient, a boundary, a regional pocket) without being differentially expressed between any of the groups you defined. Moran's I measures spatial autocorrelation, which is whether nearby bins have more similar expression than you would expect. High Moran's I picks out spatial organization that group-based DE would never even look for.

### 7. Cell-cell communication (ligand-receptor)

Cells communicate through ligand-receptor interactions, and spatial transcriptomics gives you a way to ask whether a ligand and its receptor actually show up in neighboring cell types, not just somewhere in the dataset. We use **LIANA**, a meta-method that wraps five established L-R scoring methods (CellPhoneDB, CellChat, NATMI, Connectome, SingleCellSignalR) and aggregates their per-pair ranks into a single consensus score. Each method defines a "strong" interaction a little differently (mean expression, expression product, specificity, logFC, scaled weight), so aggregating across them cancels out the per-method biases. The Visium HD colorectal cancer paper this tutorial follows ([de Oliveira et al., Nature Genetics 2025](https://www.nature.com/articles/s41588-025-02193-3)) ran LIANA in its Figure 5 cell-cell communication analysis, so this matches the published workflow.

### 8. Validation

Computational annotations and spatial patterns are only as good as the independent confirmation you can put behind them. Comparing marker-based annotations against published ground-truth labels checks whether your pipeline recovers known biology. Plotting marker gene expression across annotated clusters (with a dotplot, say) is the other half of the test, confirming the annotations are biologically coherent and not just an artifact of the clustering or scoring approach.
