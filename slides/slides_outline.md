# Visium HD Spatial Transcriptomics: Slide Outline (5 slides)

## Slide 1: Why Spatial Transcriptomics, and the Landscape

scRNA-seq dissociates tissue, so you know what cells express but not where they sit. Spatial transcriptomics measures expression in situ and keeps the coordinates. In a tumor, immune-cell position relative to tumor cells predicts therapy response. That signal is gone after dissociation. (Nature Method of the Year 2020.)

### Two technology families

| Feature        | Sequencing-based                                 | Imaging-based                          |
|----------------|--------------------------------------------------|----------------------------------------|
| Examples       | Visium, Visium HD, Slide-seq, Stereo-seq         | MERFISH, seqFISH+, Xenium, CosMx       |
| Gene coverage  | Whole transcriptome (~18k+ genes)                | Panels (100 to 6,000 genes)            |
| Resolution     | 2 to 55 um (Visium HD: 2 um)                     | Subcellular (~100 nm)                  |
| Tissue         | FFPE and fresh frozen                            | Mostly fresh frozen                    |
| Best for       | Discovery, FFPE archives, unbiased transcriptome | Validation, panels, subcellular signal |

### Trade-offs

| Advantages                                   | Limitations                                 |
|----------------------------------------------|---------------------------------------------|
| Preserves tissue architecture                | Resolution vs. gene-coverage trade-off      |
| Neighborhood and L-R analysis in context     | 2D sections miss 3D organization            |
| Integrates with H&E / IF                     | Millions of bins, compute-hungry            |
| FFPE unlocks clinical archives (Visium HD)   | ~$1,500 to $2,500 per sample, no time axis  |

---

## Slide 2: Visium HD (Arrays, Binning, Segmentation)

Continuous 2 um barcoded oligos over a 6.5 mm square, transferred onto tissue by the **CytAssist** instrument. Probe-based chemistry, works on FFPE and fresh frozen, Illumina sequencing downstream. No gaps between capture units.

### Space Ranger converts the raw array into analyzable output

- **Binning** aggregates 2 um barcodes into square bins:
  - 2 um: max resolution, very sparse. Rarely used directly.
  - **8 um: ~1 to 3 cells per bin. Primary analysis resolution.**
  - 16 um: ~5 to 15 cells. Robust, good for QC and exploration.
- **Segmentation (v4.0+)**: StarDist runs on the H&E, detects nuclei, and assigns barcodes to individual cells. Output is one row per cell instead of per bin.

### Visium v1 vs. Visium HD

| Feature        | Visium v1              | Visium HD             |
|----------------|------------------------|-----------------------|
| Resolution     | 55 um spots            | 2 um continuous       |
| Cells / unit   | ~5 to 30 per spot      | ~1 to 3 per 8 um bin  |
| Gaps           | 100 um between spots   | None                  |
| FFPE           | Yes (v2)               | Yes (primary)         |
| Segmentation   | No                     | Yes (StarDist)        |

---

## Slide 3: Space Ranger in Practice (Command, Outputs, Formats)

Everything spatial downstream is the standard 10x MEX matrix plus two extra files that tell you where each barcode sits on the tissue. (Terminal blocks are representative; numbers illustrative.)

### 3a. Run

```bash
spaceranger count \
    --id=CRC_P1 \
    --transcriptome=/refs/refdata-gex-GRCh38-2024-A \
    --probe-set=/refs/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2024-A.csv \
    --fastqs=/data/fastqs/CRC_P1 \
    --image=/data/histology/CRC_P1_HE.tif \
    --cytaimage=/data/cytassist/CRC_P1_cyta.tif \
    --slide=H1-XXXXXXX --area=A1 \
    --create-bam=false --localcores=16 --localmem=128
```

Writes `CRC_P1/outs/`. Runtime 2 to 8 hours.

### 3b. `tree CRC_P1/outs/`

```text
CRC_P1/outs/
├── binned_outputs/                         [spatial]
│   ├── square_002um/                       2 um bins (sparse)
│   ├── square_008um/                       8 um bins (primary)
│   │   ├── filtered_feature_bc_matrix/     standard 10x MEX (see 3c)
│   │   ├── filtered_feature_bc_matrix.h5   same, HDF5
│   │   └── spatial/                        [spatial]
│   │       ├── tissue_positions.parquet    join key (see 3d)
│   │       ├── scalefactors_json.json
│   │       └── tissue_{hires,lowres}_image.png
│   └── square_016um/                       16 um bins (exploration)
├── segmented_outputs/                      [spatial, v4.0+]
│   ├── filtered_feature_cell_matrix.h5
│   └── nucleus_segmentations.geojson       StarDist polygons (see 3d)
├── cloupe_008um.cloupe
├── web_summary.html
└── metrics_summary.csv
```

The only new pieces vs. plain 10x scRNA-seq: `binned_outputs/`, `spatial/tissue_positions.*`, and `segmented_outputs/*.geojson`.

### 3c. MEX triplet (standard 10x, not spatial)

```text
$ head features.tsv
ENSG00000187634   SAMD11    Gene Expression
ENSG00000188976   NOC2L     Gene Expression
ENSG00000187961   KLHL17    Gene Expression
```
Ensembl ID, symbol, feature type. `scanpy.read_10x_h5` uses column 2 as `var_names`.

```text
$ head barcodes.tsv
s_008um_00000_00001-1
s_008um_00000_00002-1
s_008um_00000_00003-1
```
Format: `s_{binsize}_{row}_{col}-{slide}`. The `(row, col)` integers index the 2 um grid.

```text
$ head matrix.mtx
%%MatrixMarket matrix coordinate integer general
%
18085 40000 2431207
1 1 3
1 14 1
2 1 1
```
Sparse triplet. Header line: `n_genes n_barcodes n_nonzero`. Data rows: `gene_idx barcode_idx umi_count` (1-indexed). `read_10x_mtx` loads this into `adata.X` as CSR, shape `n_barcodes x n_genes`.

None of 3c is spatial. Same format as a 10x Chromium scRNA-seq run.

### 3d. What makes it spatial

**`tissue_positions.parquet`** (bin to pixel map):

```text
                 barcode  in_tissue  array_row  array_col  pxl_row_in_fullres  pxl_col_in_fullres
0  s_008um_00000_00001-1          1          0          0              4213.0              2187.0
1  s_008um_00000_00002-1          1          0          1              4213.0              2195.0
2  s_008um_00000_00003-1          1          0          2              4213.0              2203.0
3  s_008um_00000_00004-1          0          0          3              4213.0              2211.0
```
`in_tissue=0` bins are dropped by the loader. `pxl_row/col_in_fullres` become `adata.obsm["spatial"]`. See `scripts/utils.py:91-103`.

**`nucleus_segmentations.geojson`** (one polygon per nucleus):

```text
{ "type": "FeatureCollection",
  "features": [
    { "type": "Feature",
      "properties": { "cell_id": 1, "area_px": 412 },
      "geometry": { "type": "Polygon",
        "coordinates": [[ [4231.5,2198.0],[4232.1,2199.7],[4233.4,2200.8], ... ]] }},
    ...
  ]}
```
Loader averages each polygon's vertices to get a centroid, which becomes the cell's coordinate. Barcodes are `cellid_XXXXXXXXX-1`, the integer matches `properties.cell_id`. See `scripts/utils.py:164-170`.

### 3e. Downstream

1. MEX matrix into `AnnData`.
2. `tissue_positions` or geojson centroids into `adata.obsm["spatial"]`.
3. Filter `in_tissue == 1`.
4. Standard scanpy / squidpy analysis with coordinates attached.

---

## Slide 4: When to Use Visium HD

Right call for unbiased whole-transcriptome spatial discovery on FFPE at near-single-cell resolution. Wrong call if you already have a gene panel or need live / subcellular data.

### Ideal

- **Whole-transcriptome discovery**: no prior gene list, finding novel spatial markers.
- **FFPE archives**: clinical biobanks that cannot go on fresh-frozen platforms.
- **Near-single-cell**: 8 um is ~1 to 3 cells; segmentation is one.
- **Heterogeneous tissue**: tumor, stroma, immune infiltrate in the same section.
- **Histology integration**: direct H&E overlay.

### Not ideal

- **You already have a panel**: Xenium or CosMx are faster and cheaper.
- **Subcellular localization**: use MERFISH or seqFISH+.
- **Live-cell dynamics**: fixed tissue only.
- **3D reconstruction**: 2D sections only.
- **Large cohorts (>50 samples)**: cost and compute dominate.

### Numbers

| Axis                 | Typical                                     |
|----------------------|---------------------------------------------|
| Data size            | ~5 to 50 GB per sample                      |
| RAM                  | 16+ GB for 8 um, 64+ GB for 2 um            |
| Space Ranger runtime | 2 to 8 hours                                |
| Downstream runtime   | Hours to days                               |
| Cost per sample      | ~$1,500 to $2,500                           |

---

## Slide 5: Analysis Pipeline and QC Gotchas

Standard scanpy, plus spatial-aware QC and analysis. A few things trip people up on this platform.

### Pipeline

1. QC and filtering (spatial-aware).
2. Normalization and HVG selection.
3. PCA and UMAP.
4. Cell-type annotation: deconvolution (binned) or label transfer / marker scoring (segmented).
5. Spatial analysis: neighborhood enrichment, SVGs (Moran's I).
6. Cell-cell communication (ligand-receptor).
7. Validation against ground truth / orthogonal data.

### Resolution

- **16 um**: start here for exploration and QC.
- **8 um**: primary analysis (what this tutorial uses).
- **2 um**: only for specific high-resolution questions. Very sparse.
- **Segmented**: when you need true single-cell resolution.

### QC gotchas specific to Visium HD

- **Use spatial QC, not just per-bin thresholds.** Tissue folds, debris, and edges create spatially-correlated artifacts. Use Moran's I on QC metrics and neighbor-based outlier detection (`spatial_outlier_detection` in `scripts/utils.py`).
- **High-mt is not automatically bad.** Clustered at edges or folds: damage, filter. Dispersed through the tumor: likely real hypoxia or necrosis, keep.
- **`pct_counts_ribo` is flat zero on this platform.** The 10x Human WTA probe panel excludes ribosomal protein genes (RPL\*, RPS\*) to avoid probe saturation. Expected, not a bug. Use `pct_counts_mt` as the metabolic QC axis. `compute_qc_metrics` in `scripts/utils.py:204` logs this.
- **Thresholds are resolution-dependent.** Do not reuse cutoffs across 2, 8, and 16 um.
- **Watch segmentation artifacts.** Over / under-segmentation, nuclear vs. cytoplasmic capture. Sanity-check cell counts against the H&E.

### Deconvolution vs. segmentation

- Binned (8 or 16 um): multiple cells per bin, use deconvolution (cell2location, RCTD, STdeconvolve).
- Segmented: one cell per row, use direct annotation (label transfer, marker scoring).
- If possible, run both and check concordance.

### What this tutorial demonstrates

Real QC on a human CRC dataset, ground-truth comparison, neighborhood enrichment, SVGs, and ligand-receptor analysis, using the file formats from slide 3.
