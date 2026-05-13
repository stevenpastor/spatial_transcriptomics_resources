# Visium HD Spatial Transcriptomics

## 1. Why Spatial Transcriptomics

### The problem with dissociation

Single-cell and single-nucleus RNA-seq (scRNA-seq, snRNA-seq) measure what each cell expresses, but they require dissociating the tissue first. Once cells are separated from the section, all spatial information is lost: you no longer know where a cell sat, what it was next to, or which tissue compartment it belonged to. For many biological questions, this is a critical gap.

Cell function depends on context. The same fibroblast behaves differently in a wound margin vs. deep stroma. Immune cells adjacent to tumor cells engage different programs than those in a distant lymphoid aggregate. Gradients of signaling molecules, physical cell-cell contacts, and local oxygen levels all shape gene expression, and all of these are spatial. Dissociation-based methods cannot recover this information after the fact.

### What spatial transcriptomics adds

Spatial transcriptomics measures gene expression in intact tissue sections, preserving both the expression profile and the x,y coordinates of each measurement. This enables analyses that are impossible with scRNA-seq:

- Which cell types are neighbors, and does that co-localization happen more than expected by chance?
- Are specific genes expressed in spatial gradients, boundaries, or confined to particular tissue regions?
- Which ligand-receptor signaling pairs are active between spatially adjacent cell types?
- How does the cellular microenvironment change across a tumor-stroma boundary or along an organ axis?

### Two technology families

Spatial methods split into two camps. **Sequencing-based** methods capture mRNA on a patterned slide of barcoded oligonucleotides, then sequence the recovered cDNA. You get the whole transcriptome at whatever spatial resolution the slide pattern allows. **Imaging-based** methods detect a predefined panel of genes directly in the tissue using rounds of fluorescent probe hybridization and imaging. You get subcellular spatial precision but only on the genes you chose up front.

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

A few rules of thumb. If you need whole-transcriptome discovery on archival FFPE blocks, Visium HD is the default. If you have a focused gene list and need true subcellular resolution with cell boundaries, pick Xenium, MERFISH/MERSCOPE, or CosMx based on which panel fits your biology. Standard Visium is still useful for pilot work or when sample numbers matter more than within-section resolution. Stereo-seq and Slide-seqV2 fill the niche of high-resolution whole-transcriptome work on fresh frozen material.

This tutorial uses Visium HD because it works on FFPE, captures the full transcriptome, and ships with a polished pipeline (Space Ranger) and viewer (Loupe Browser) that most labs already have access to. The same analysis ideas (QC, clustering, neighborhood analysis, spatially variable genes, ligand-receptor) carry over to the other technologies with small adjustments.

---

## 2. How Visium HD Works

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

Even at 8 um, each bin may contain parts of multiple cells, so gene expression is a mixture. **Segmentation** (Space Ranger v4.0+) solves this by running StarDist nuclei detection on the H&E image to draw a polygon around each nucleus, then assigning barcodes to individual cells. The result is one row per cell instead of one row per bin, giving true single-cell resolution. This matters for cell type annotation and any analysis where mixing signal from neighboring cells would blur the result. The segmented output is its own deliverable from Space Ranger and has a few caveats worth knowing before you commit to it, covered in the "Segmentation output" subsection below.

---

## 3. Space Ranger: From Raw Data to Spatial Matrices

### Inputs

Space Ranger needs three inputs to produce spatial gene expression data:

| Input | What it is | Why it is needed |
|-------|-----------|-----------------|
| **FASTQ files** | Raw sequencing reads from Illumina | Contain the spatial barcode, UMI, and gene sequence for every captured molecule |
| **H&E image** | High-resolution microscope image of the stained tissue section | Provides the tissue morphology context; bins are mapped back to pixel coordinates on this image for all downstream spatial plots |
| **CytAssist image** | Image taken by the CytAssist instrument during probe transfer | Used to align the tissue position on the capture slide; connects where the tissue sat on the oligo lawn to the barcode grid |

Three more inputs are required and easy to get wrong:

- **Reference transcriptome**: a prebuilt 10x reference matching your species (human GRCh38 or mouse GRCm39). Download from the 10x Genomics support site; do not build your own unless you have a specific reason. Pick the version that matches the probe set release notes.
- **Probe set CSV**: the list of probe sequences used in your assay. The probe set must match the chemistry (Visium HD has its own probe set, distinct from standard Visium). 10x ships the CSV with each kit and posts it for download. If the CSV does not match the kit, Space Ranger will silently produce nonsense counts.
- **Slide serial number and capture area ID**: printed on the slide. Space Ranger uses these to look up the correct fiducial layout. If you transcribe the serial wrong, registration will fail.

### Image registration: Loupe Browser

Before running Space Ranger, the H&E image and CytAssist image need to be aligned so that barcodes map to the correct tissue coordinates. Space Ranger can auto-detect this alignment in most cases, and you should let it try first. Manual alignment in Loupe Browser is the fallback when auto-alignment goes wrong.

**When do you need manual alignment?** A few situations push you out of the auto path: the tissue section was rotated or placed off-center on the slide, the H&E staining is faint or has uneven contrast, debris or air bubbles obscure the fiducial frame, or Space Ranger's web summary reports low alignment confidence after a first auto-attempt. If your spatial plots look smeared, shifted, or rotated relative to the H&E, that is also a flag to go back and check alignment.

The Loupe Browser workflow has three steps:

1. **Upload the CytAssist image and identify the fiducial markers.** Fiducials are the small circular dots arranged in a frame around the capture area; they are the ground truth for where the barcode grid sits. Loupe will try to detect them automatically. Fix any that landed on tissue or debris instead of on the dots.
2. **Upload the H&E microscope image and pin matching landmarks on both images.** A good landmark pair is a feature visible in both images: a vessel cross-section, a duct, a fold edge, or the corner of a tissue region. Aim for at least four well-spread pairs covering different parts of the capture area; clustered landmarks underdetermine the transform.
3. **Assess alignment quality and optionally auto-refine.** Loupe shows the overlay with a quality readout. Look for the H&E features sitting cleanly on top of the CytAssist counterparts. If a corner of the section is off, add a landmark there and re-refine.

![Loupe Browser: Identify fiducials](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/06_identify_fiducials.png)

![Loupe Browser: Assess alignment](https://cdn.10xgenomics.com/image/upload/v1708554526/software-support/Spatial-GEX/LB-v8.0/hd-manual-alignment/18_assess_alignment.png)

The exported alignment file is a small JSON. You pass it to Space Ranger with `--loupe-alignment=path/to/alignment.json` and the pipeline uses your manual registration instead of auto-detection. The JSON is portable and worth committing to a project repo so the registration step is reproducible. Space Ranger refuses to run if the JSON references a slide serial that does not match the FASTQs or the CytAssist image, so a mismatched file fails fast rather than producing silently wrong coordinates.

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

### Which outputs to actually use

Space Ranger drops a lot of files. Most analyses only touch a handful. Here is what to grab and why:

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

A few outputs to skip by default. The `square_002um/` matrix is the raw resolution but very sparse, so most analyses aggregate it back up to 8 um (or use segmentation) rather than working with it directly. The `raw_feature_bc_matrix` (without `filtered_` in the name) includes bins outside the tissue boundary; only reach for it if you need to debug the in-tissue filter. The `square_016um/` matrix is good for quick spatial QC plots but loses too much resolution for serious cell typing.

### Segmentation output

Segmentation gives you a cell-level matrix instead of a bin-level matrix. Space Ranger v4.0+ runs StarDist on the H&E image to detect nuclei, draws a polygon around each one, and assigns each underlying 2 um barcode to the cell whose nucleus polygon contains it (with some local rules for cytoplasmic capture). The result is one row per nucleus.

Inside `segmented_outputs/` you get:

- `filtered_feature_cell_matrix.h5`: the cell-level count matrix.
- `nucleus_segmentations.geojson`: the nucleus polygons. Useful for plotting cell outlines and for computing centroids when you want spatial coordinates in your AnnData. `utils.load_visium_hd_segmented` reads this file and uses the polygon centroids as `obsm["spatial"]`.
- A per-cell metrics CSV with counts, gene counts, and nucleus area.

The trade-offs to know up front. Segmentation produces cleaner per-cell expression and matches scRNA-seq conventions, so cell type annotation tools tend to work better on it. The cost is that any cell whose nucleus the segmenter missed gets its transcripts dropped or assigned to a neighbor, and the totals you measure are not absolute (they depend on the segmenter as much as the tissue). Use segmentation when you need cell identity. Stick to bins when you care about quantitative tissue-level expression patterns.

If the bundled StarDist segmentation does not work well on your tissue, you have alternatives:

- **bin2cell** (Python): takes the Space Ranger 2 um output plus an external nucleus mask and aggregates bins into cells. This is the standard escape hatch when you want a different segmenter.
- **StarDist (general)**: a nucleus segmentation model. Run it on the H&E image to get a nucleus mask, then feed that mask to bin2cell.
- **Cellpose**: a general-purpose cell and nucleus segmentation model. Often used for whole-cell masks from H&E or membrane stains; the resulting mask plugs into bin2cell the same way.

The choice depends on your tissue: StarDist is the default for H&E nuclei, Cellpose handles harder cases and other stains better.

Once you have a cell-level AnnData (whether from `segmented_outputs/` or from bin2cell), downstream analysis looks the same as scRNA-seq: normalize, find highly variable genes, PCA, neighbors, UMAP, Leiden, annotate. The spatial layer comes back in via `squidpy.gr.spatial_neighbors`, which now operates on cell centroids instead of bin centroids. Every analysis covered in Section 4 still applies.

### What you will learn in the quickstart

The quickstart notebook (`notebook/visium_hd_quickstart.ipynb`) walks through:

1. Loading and exploring Visium HD spatial transcriptomics data
2. Quality control with filtering (per-bin metrics plus spatial outlier detection)
3. Comparing marker-based annotations against published ground-truth labels
4. Spatial neighborhood analysis, spatially variable genes, and cell-cell communication

To keep Colab runtimes short, the notebook loads pre-processed data rather than running Space Ranger end-to-end. The files it pulls from Figshare are:

- **`adata_8um.h5ad`**: the 8 um binned count matrix (`filtered_feature_bc_matrix.h5`) already converted to an AnnData object with spatial coordinates from `tissue_positions.parquet` stored in `adata.obsm["spatial"]`. This is the primary object for all analysis.
- **`adata_8um_annotated.h5ad`**: same object after normalization, PCA, UMAP, Leiden clustering, and cell type annotation have been run. Pre-computed so the Colab notebook runs in minutes instead of hours.
- **`ground_truth_labels.csv`**: published cell type labels from the original study (de Oliveira et al., Nature Genetics 2025) for validation.
- **H&E image**: the tissue image for spatial overlay plots.

We use the 8 um resolution because it balances signal density (~1 to 3 cells per bin) with manageable data size.

---

## 4. Python Data Structures and Tooling

Before opening the notebook, it helps to know the three libraries you will see on every cell: `anndata`, `scanpy`, and `squidpy`. The R equivalents (Seurat, SpatialExperiment) will get their own primer in a follow-up doc; everything here is Python.

### AnnData: one object, many slots

scRNA-seq and spatial datasets both need the same things stored together: a count matrix, per-cell metadata, per-gene metadata, reduced-dimension embeddings, and arbitrary side data (images, neighbor graphs, parameters). AnnData is the standard container that holds all of that in one object with aligned indices.

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

Indexing works like a DataFrame: `adata[adata.obs["cell_type"] == "T cell"]` returns a view of only those cells, and `adata[:, "EPCAM"]` selects a single gene.

A spatial AnnData and an scRNA-seq AnnData have the same structure. The only differences are that the spatial one populates `obsm["spatial"]` with x,y coordinates and usually stores the tissue image plus scale factors under `uns["spatial"]`. Every other operation (normalization, clustering, DE, embeddings) works identically on either kind of object. That is the point of the format: tools written for one work on the other.

Reading and writing is `sc.read_h5ad("file.h5ad")` and `adata.write_h5ad("file.h5ad")`. The tutorial data ships as `.h5ad` files, which is just HDF5 with an AnnData schema on top.

### scanpy: the analysis workhorse

`scanpy` provides the standard pipeline (QC, normalization, highly variable gene selection, PCA, neighbors, UMAP, Leiden, differential expression, plotting) operating directly on AnnData. The canonical flow looks like:

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

For plotting, `sc.pl.umap(adata, color="leiden")` makes the UMAP scatter, `sc.pl.violin` makes per-cluster expression plots, and `sc.pl.spatial(adata, color="cell_type")` makes the spatial overlay using `obsm["spatial"]` and the H&E in `uns["spatial"]`.

### squidpy: the spatial layer on top

`squidpy` adds spatial-aware analysis on top of the scanpy/AnnData stack: neighbor graphs from coordinates, neighborhood enrichment, spatial autocorrelation (Moran's I), and image-based features. The foundation is the spatial neighbor graph:

```python
import squidpy as sq

sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=8)
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
sq.gr.spatial_autocorr(adata, mode="moran")
```

Once `spatial_neighbors` runs, the connectivity graph lives in `adata.obsp["spatial_connectivities"]` and every downstream squidpy call uses it. This is the spatial analog of `sc.pp.neighbors`, which builds an expression-space neighbor graph in `adata.obsp["connectivities"]`. Both can coexist on the same AnnData.

### Adjacent libraries worth knowing

- `anndata`: the library that defines the AnnData class. You rarely import it directly; scanpy re-exports what you need.
- `bin2cell`: bin-to-cell aggregation for Visium HD when you want to swap in a custom nucleus segmenter.
- `cell2location`, `tangram`, `RCTD` (Python ports): deconvolution methods that estimate cell type fractions per bin when resolution is below single-cell.
- `spatialdata` and `spatialdata-io`: an emerging unified format for multi-modal spatial datasets (images, points, shapes, tables). Useful when you start mixing modalities; AnnData is still the practical default for plain Visium HD analysis.

---

## 5. Analysis Pipeline

### Before analysis: upstream processing

By the time we open a notebook, several steps have already happened:

1. **Tissue preparation and imaging**: the tissue was sectioned onto a glass slide, stained with H&E, and imaged at high resolution. This H&E image is both the morphology reference and the coordinate system for all spatial plots.

2. **CytAssist and sequencing**: the tissue section was run through CytAssist to transfer probes onto the Visium HD capture slide, libraries were built, and the sample was sequenced on an Illumina instrument. The output is FASTQ files.

3. **Image registration**: the H&E microscope image and the CytAssist image need to be aligned so barcodes map to the correct tissue positions. Space Ranger can auto-detect this alignment in most cases. For difficult samples (rotated sections, low contrast, or poor fiducial detection), Loupe Browser provides a manual alignment workflow where you identify fiducials and pin landmarks between the two images. If manual alignment was used, the exported JSON file is passed to Space Ranger via `--loupe-alignment`.

4. **Space Ranger**: takes the FASTQs, H&E image, CytAssist image (and alignment file if applicable), reference transcriptome, and probe set as inputs. It demultiplexes reads, maps them to the genome, registers images, and outputs binned count matrices at multiple resolutions (2, 8, 16 um) along with spatial coordinate files and a QC report. This is the last "black box" step before interactive analysis begins.

Everything from this point forward is what we do in the tutorial, starting from the Space Ranger outputs loaded as AnnData objects (see Section 4 for the data structures and libraries involved).

---

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

**Why**: cells communicate through ligand-receptor interactions, and spatial transcriptomics lets you test whether a ligand and its receptor are actually expressed in neighboring cell types, not just co-expressed somewhere in the dataset. We use **LIANA**, a meta-method that wraps five established L-R scoring methods (CellPhoneDB, CellChat, NATMI, Connectome, SingleCellSignalR) and aggregates their per-pair ranks into a single consensus score. Each method captures a slightly different definition of a "strong" interaction (mean expression, expression product, specificity, logFC, scaled weight), so aggregating reduces method-specific biases. The Visium HD colorectal cancer paper this tutorial uses ([de Oliveira et al., Nature Genetics 2025](https://www.nature.com/articles/s41588-025-02193-3)) ran LIANA in its Figure 5 cell-cell communication analysis, so this matches the published workflow.

### 8. Validation

**Why**: computational annotations and spatial patterns need independent confirmation. Comparing marker-based annotations against published ground truth labels checks whether your pipeline recovers known biology. Examining marker gene expression across annotated clusters (e.g., dotplots) confirms that annotations are biologically coherent and not artifacts of the clustering or scoring approach.
