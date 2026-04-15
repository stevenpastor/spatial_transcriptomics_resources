# Visium HD Spatial Transcriptomics — Slide Outline (5 slides)

## Slide 1: Why Spatial Transcriptomics + Technology Landscape

**Key Message**: Spatial transcriptomics preserves *where* cells are while measuring *what* they express — closing the gap that scRNA-seq opens by dissociating tissue. Two technology families dominate, with different trade-offs.

### Why it matters

- **Problem**: scRNA-seq dissociates tissue → you know WHAT cells express, but not WHERE. Cell function depends on location, neighbors, and microenvironment.
- **Solution**: Measure gene expression *in situ*, on an intact section.
- **Key insight**: WHERE + WHAT = mechanism. Example: in a tumor, immune-cell position relative to tumor cells predicts therapy response — unrecoverable after dissociation.
- **Recognition**: Nature *Method of the Year 2020* (lineage: ISH → FISH → ST).

### Two technology families

| Feature           | Sequencing-based                                  | Imaging-based                              |
|-------------------|---------------------------------------------------|--------------------------------------------|
| **Examples**      | Visium, Visium HD, Slide-seq, Stereo-seq          | MERFISH, seqFISH+, Xenium, CosMx, CODEX    |
| **Gene coverage** | Whole transcriptome (~18,000+ genes)              | Targeted panels (100–6,000 genes)          |
| **Resolution**    | 2–55 µm (Visium HD: 2 µm)                         | Subcellular (~100 nm)                      |
| **Tissue**        | FFPE + Fresh Frozen                               | Primarily Fresh Frozen (FFPE expanding)    |
| **Quantification**| UMI-based, digital                                | Spot counting, analog                      |
| **Best for**      | Discovery, FFPE archives, unbiased transcriptome  | Validation, subcellular signal, panels     |

### Advantages ↔ Limitations (at a glance)

| ✅ Advantages                                        | ⚠️ Limitations                                          |
|------------------------------------------------------|---------------------------------------------------------|
| Preserves tissue architecture                        | Resolution ↔ gene-coverage trade-off                    |
| Enables neighborhood / L-R analysis in context       | 2D sections miss 3D organization                        |
| Integrates with histology (H&E, IF)                  | Large datasets (millions of bins), compute-hungry       |
| FFPE-compatible (Visium HD) → unlocks clinical archives | ~$1,500–2,500 per sample; no temporal info           |

---

## Slide 2: Visium HD Deep Dive — Arrays, Binning, Segmentation

**Key Message**: Visium HD reaches near-single-cell resolution with whole-transcriptome coverage by combining a **continuous 2 µm barcoded array** with **computational binning and nuclei segmentation** in Space Ranger v4.0+.

### The instrument and the array

- **CytAssist**: automated instrument that stamps tissue onto the barcoded array, compatible with standard H&E / IF workflows.
- **Array**: continuous 2 µm barcoded oligos over a 6.5 mm × 6.5 mm capture area — *no gaps* between capture units (unlike v1 Visium spots).
- **Chemistry**: probe-based, works on **FFPE and fresh frozen**; standard Illumina short-read sequencing downstream.

### How Space Ranger turns raw data into analyzable output

- **Binning** — aggregates 2 µm barcodes into square bins:
  - **2 µm**: maximum resolution, very sparse (<10 UMIs most bins, rarely used directly)
  - **8 µm**: ~1–3 cells per bin, **primary analysis resolution**
  - **16 µm**: ~5–15 cells, robust signal, good for exploration / QC
- **Segmentation** (new in v4.0) — runs StarDist on the H&E image to detect nuclei, then assigns barcodes to individual cells. Output is one row per segmented cell, enabling true single-cell spatial analysis instead of per-bin aggregates.

### Visium (v1) vs Visium HD

| Feature             | Visium v1                      | Visium HD                         |
|---------------------|--------------------------------|-----------------------------------|
| Resolution          | 55 µm spots                    | 2 µm continuous array             |
| Cells per unit      | ~5–30 per spot                 | ~1–3 per 8 µm bin                 |
| Gaps                | 100 µm center-to-center        | Continuous (none)                 |
| FFPE support        | Yes (v2)                       | Yes (primary)                     |
| Segmentation        | No                             | Yes (StarDist, v4.0+)             |

→ This whole slide sets up slide 3, where we actually look at what Space Ranger writes to disk.

---

## Slide 3: Space Ranger in Practice — Command → Outputs → Files

**Key Message**: Everything spatial you ever do downstream is just the standard 10x MEX matrix plus two extra files that tell you **where each barcode lives on the tissue**. This slide shows them literally.

> *Terminal blocks below are representative — file layout, column names, and formats match a real Space Ranger v4.0 run, but numbers are illustrative. Real local data lives on Figshare for this tutorial.*

### 3a. Running Space Ranger

```bash
spaceranger count \
    --id=CRC_P1 \
    --transcriptome=/refs/refdata-gex-GRCh38-2024-A \
    --probe-set=/refs/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv \
    --fastqs=/data/fastqs/CRC_P1 \
    --image=/data/histology/CRC_P1_HE.tif \
    --cytaimage=/data/cytassist/CRC_P1_cyta.tif \
    --slide=H1-XXXXXXX --area=A1 \
    --create-bam=false \
    --localcores=16 --localmem=128
```

One command. It reads the FASTQs, aligns against the probe set, registers the H&E image, and writes a `CRC_P1/outs/` directory that scanpy / squidpy can load directly. Typical runtime: **2–8 hours**.

### 3b. What `tree CRC_P1/outs/` looks like

```text
CRC_P1/outs/
├── binned_outputs/                              ←  spatial-specific
│   ├── square_002um/                            ←  2 µm bins  (highest res, sparsest)
│   ├── square_008um/                            ←  8 µm bins  (primary analysis)
│   │   ├── filtered_feature_bc_matrix/          ←  standard 10x MEX triplet  (§3c)
│   │   ├── filtered_feature_bc_matrix.h5        ←  same data, HDF5
│   │   └── spatial/                             ←  spatial-specific
│   │       ├── tissue_positions.parquet         ←  THE spatial join key  (§3d)
│   │       ├── scalefactors_json.json
│   │       ├── tissue_hires_image.png
│   │       └── tissue_lowres_image.png
│   └── square_016um/                            ←  16 µm bins (exploration)
├── segmented_outputs/                           ←  spatial-specific (v4.0+)
│   ├── filtered_feature_cell_matrix.h5
│   └── nucleus_segmentations.geojson            ←  StarDist polygons  (§3d)
├── cloupe_008um.cloupe                          ←  Loupe browser
├── web_summary.html                             ←  QC report (open in browser)
└── metrics_summary.csv
```

**Compared to plain 10x scRNA-seq**, the only *new* things are: `binned_outputs/` (multiple resolutions), `spatial/tissue_positions.*` (bin → pixel mapping), and `segmented_outputs/*.geojson` (nucleus polygons). Everything else is the standard scRNA-seq file set.

### 3c. MEX format — what `head` actually shows

```text
$ head square_008um/filtered_feature_bc_matrix/features.tsv
ENSG00000187634   SAMD11      Gene Expression
ENSG00000188976   NOC2L       Gene Expression
ENSG00000187961   KLHL17      Gene Expression
ENSG00000187583   HES4        Gene Expression
ENSG00000187642   ISG15       Gene Expression
ENSG00000188157   AGRN        Gene Expression
```
→ Ensembl ID | gene symbol | feature type. `scanpy.read_10x_h5` uses column 2 (symbol) as `var_names` — that's why `adata_8um.var_names[:5]` in the notebook returns `['SAMD11', 'NOC2L', 'KLHL17', 'HES4', 'ISG15']`.

```text
$ head square_008um/filtered_feature_bc_matrix/barcodes.tsv
s_008um_00000_00001-1
s_008um_00000_00002-1
s_008um_00000_00003-1
s_008um_00000_00004-1
s_008um_00000_00005-1
```
→ Barcode format: `s_{binsize}_{row}_{col}-{slide}`. The `(row, col)` integers index the 2 µm array grid — *this is the hook that lets Space Ranger reassemble spatial coordinates later*.

```text
$ head square_008um/filtered_feature_bc_matrix/matrix.mtx
%%MatrixMarket matrix coordinate integer general
%
18085 40000 2431207
1 1 3
1 14 1
2 1 1
3 1 2
```
→ Sparse triplet format. The third header line reads `n_genes n_barcodes n_nonzero`. Each data row is `gene_idx barcode_idx umi_count` (1-indexed). `scanpy.read_10x_mtx` parses this into `adata.X` (CSR sparse matrix, shape `n_barcodes × n_genes`).

**None of this is spatial yet** — it's the exact same format a 10x Chromium scRNA-seq run writes. Spatial information lives in the next two files.

### 3d. What makes it spatial

**`tissue_positions.parquet` — the bin → pixel map**

```text
$ python -c "import pandas as pd; print(pd.read_parquet('spatial/tissue_positions.parquet').head())"
                 barcode  in_tissue  array_row  array_col  pxl_row_in_fullres  pxl_col_in_fullres
0  s_008um_00000_00001-1          1          0          0              4213.0              2187.0
1  s_008um_00000_00002-1          1          0          1              4213.0              2195.0
2  s_008um_00000_00003-1          1          0          2              4213.0              2203.0
3  s_008um_00000_00004-1          0          0          3              4213.0              2211.0
4  s_008um_00000_00005-1          1          0          4              4213.0              2219.0
```

**This is the spatial join key.** For every barcode in the MEX matrix, it gives:
- `in_tissue` — 1 if the bin sits under the tissue section, 0 if off-tissue (the loader drops `0` rows; see `scripts/utils.py:101–103`)
- `array_row`, `array_col` — grid coordinates on the 2 µm array
- `pxl_row_in_fullres`, `pxl_col_in_fullres` — **pixel coordinates on the full-resolution H&E image**, which is what `adata.obsm["spatial"]` stores and what every downstream plot uses (`scripts/utils.py:91–93`).

**`nucleus_segmentations.geojson` — the StarDist output**

```text
$ head -30 segmented_outputs/nucleus_segmentations.geojson
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": { "cell_id": 1, "area_px": 412 },
      "geometry": {
        "type": "Polygon",
        "coordinates": [[
          [4231.5, 2198.0], [4232.1, 2199.7], [4233.4, 2200.8],
          [4234.0, 2201.9], [4233.7, 2203.1], [4232.6, 2204.0],
          ...
        ]]
      }
    },
    {
      "type": "Feature",
      "properties": { "cell_id": 2, "area_px": 387 },
      ...
    }
  ]
}
```

One GeoJSON feature per segmented nucleus — **a polygon, not a point**. The loader in `scripts/utils.py:164–170` walks `features[]`, averages each polygon's vertices to get a centroid, and uses `(centroid_y, centroid_x)` as the cell's spatial coordinate. Barcodes in the segmented matrix look like `cellid_000000001-1`, and the integer after `cellid_` matches `properties.cell_id` — that's the join.

### 3e. Bottom line

Everything downstream — scanpy, squidpy, this tutorial — is just:

1. Read the MEX matrix → `AnnData`
2. Read `tissue_positions.parquet` (or the geojson centroids) → `adata.obsm["spatial"]`
3. Filter to `in_tissue == 1`
4. Proceed with standard analysis, now with coordinates attached

That's the whole "spatial" part. Two extra files.

---

## Slide 4: When to Use Visium HD

**Key Message**: Visium HD is the right call for unbiased whole-transcriptome spatial discovery on FFPE at near-single-cell resolution. It's the wrong call if you already have a panel or need live/subcellular measurements.

### Ideal use cases

- **Whole-transcriptome discovery** — no prior gene list, want to find novel spatial markers
- **FFPE archives** — unlock spatial information from clinical biobanks that can't go on fresh-frozen platforms
- **Near-single-cell resolution needed** — 8 µm bins ≈ 1–3 cells; segmentation gets you to one
- **Heterogeneous tissue** — where cell-type composition varies across space (tumor ↔ stroma ↔ immune infiltrate)
- **Integration with histology** — direct H&E overlay, pathologist-friendly

### Not ideal for

- **Targeted validation** — if you already have a panel, Xenium / CosMx are faster and cheaper
- **Subcellular localization** — use imaging-based methods (MERFISH, seqFISH+)
- **Live-cell dynamics** — fixed tissue only, no time course
- **3D reconstruction** — 2D sections; serial-section 3D is possible but complex
- **Large cohorts (>50 samples)** — cost and compute scale become dominant

### Practical numbers

| Axis             | Typical                                                   |
|------------------|-----------------------------------------------------------|
| Data size        | ~5–50 GB per sample (depth- and bin-size dependent)       |
| RAM for analysis | 16+ GB for 8 µm; 64+ GB for 2 µm                          |
| Space Ranger runtime | 2–8 hours                                             |
| Downstream runtime   | Hours to days                                         |
| Cost per sample  | ~$1,500–2,500 (reagents + CytAssist + sequencing)         |

---

## Slide 5: Analysis Pipeline & QC Considerations

**Key Message**: Visium HD pipelines are standard scanpy plus spatial-aware QC and spatial-aware analysis steps. A few gotchas are specific to the probe panel and the binning strategy.

### Pipeline at a glance

1. **QC & filtering** — spatial-aware, not just per-bin thresholds
2. **Normalization & HVG selection** — standard scanpy
3. **Dimensionality reduction** — PCA, UMAP for visualization
4. **Cell-type annotation** — deconvolution (binned) or label transfer / marker scoring (segmented)
5. **Spatial analysis** — neighborhood enrichment, spatially variable genes (Moran's I)
6. **Cell-cell communication** — ligand-receptor analysis in spatial context
7. **Validation** — compare with ground-truth annotations, orthogonal data, known biology

### Choosing a resolution

- Start with **16 µm** for exploration and QC (best signal-to-noise).
- **8 µm** for primary analysis (best balance; what this tutorial uses).
- **2 µm** only for specific high-resolution questions — it's very sparse.
- **Segmented** outputs when you need true single-cell analysis.

### QC considerations unique to Visium HD

- **Spatial QC is essential.** Tissue folds, debris, edge effects create spatially-correlated artifacts that per-bin thresholds miss — use Moran's I on QC metrics and spatial neighbor-based outlier detection (see `scripts/utils.py` `spatial_outlier_detection`).
- **High-MT regions ≠ automatic "bad".** If high `pct_counts_mt` clusters at tissue edges or folds → damage artifact, filter. If dispersed in the tumor → likely real hypoxia / necrosis, keep.
- **⚠️ The 10x Human WTA probe panel *excludes ribosomal protein genes* (`RPL*`/`RPS*`) to prevent probe saturation.** `pct_counts_ribo` will be **flat zero** — this is expected, not a bug. Treat `pct_counts_mt` as your primary metabolic QC axis and don't spend time debugging the empty ribo panel. (`compute_qc_metrics` in `scripts/utils.py:204` now logs this explicitly.)
- **Resolution-dependent thresholds.** QC cutoffs differ across 2 / 8 / 16 µm — don't reuse numbers across resolutions.
- **Segmentation artifacts.** Over/under-segmentation, nuclear vs. cytoplasmic capture — sanity-check cell counts against the H&E.

### Deconvolution vs segmentation

- **Binned data (8/16 µm)**: multiple cells per bin → **deconvolution** (cell2location, RCTD, STdeconvolve)
- **Segmented data**: one cell per observation → **direct annotation** (label transfer, marker scoring)
- **Best practice**: run both where possible and compare for concordance.

### What this tutorial demonstrates

Real QC on a clinically relevant human CRC dataset → ground-truth cell-type comparison → neighborhood enrichment → SVGs → ligand-receptor interactions, all with the file formats you just saw on slide 3.
