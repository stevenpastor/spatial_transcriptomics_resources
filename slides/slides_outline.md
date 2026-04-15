# Visium HD Spatial Transcriptomics: Slide Outline (5 slides)

## Slide 1: Why Spatial Transcriptomics and Technology Landscape

**Key message**: Spatial transcriptomics measures gene expression in situ, preserving the tissue coordinates that scRNA-seq loses during dissociation. Two technology families dominate, with different trade-offs.

### Why it matters

- Problem: scRNA-seq dissociates tissue. You know what cells express, not where they sit. Cell function depends on location, neighbors, and microenvironment.
- Solution: measure gene expression in situ on an intact section.
- Why it changes the question you can ask: in a tumor, immune-cell position relative to tumor cells predicts therapy response. That signal is gone after dissociation.
- Recognition: Nature Method of the Year 2020. Lineage: ISH, then FISH, then spatial transcriptomics.

### Two technology families

| Feature           | Sequencing-based                                  | Imaging-based                              |
|-------------------|---------------------------------------------------|--------------------------------------------|
| Examples          | Visium, Visium HD, Slide-seq, Stereo-seq          | MERFISH, seqFISH+, Xenium, CosMx, CODEX    |
| Gene coverage     | Whole transcriptome (~18,000+ genes)              | Targeted panels (100 to 6,000 genes)       |
| Resolution        | 2 to 55 um (Visium HD: 2 um)                      | Subcellular (~100 nm)                      |
| Tissue            | FFPE and fresh frozen                             | Primarily fresh frozen (FFPE expanding)    |
| Quantification    | UMI-based, digital                                | Spot counting, analog                      |
| Best for          | Discovery, FFPE archives, unbiased transcriptome  | Validation, subcellular signal, panels     |

### Advantages and limitations

| Advantages                                          | Limitations                                             |
|------------------------------------------------------|---------------------------------------------------------|
| Preserves tissue architecture                        | Resolution vs. gene-coverage trade-off                  |
| Enables neighborhood and L-R analysis in context     | 2D sections miss 3D organization                        |
| Integrates with histology (H&E, IF)                  | Large datasets (millions of bins), compute-hungry       |
| FFPE-compatible (Visium HD), unlocks clinical archives | ~$1,500 to $2,500 per sample, no temporal information |

---

## Slide 2: Visium HD Deep Dive (Arrays, Binning, Segmentation)

**Key message**: Visium HD reaches near-single-cell resolution with whole-transcriptome coverage by combining a continuous 2 um barcoded array with computational binning and nuclei segmentation in Space Ranger v4.0 and later.

### The instrument and the array

- **CytAssist**: automated instrument that transfers tissue onto the barcoded array. Compatible with standard H&E and IF workflows.
- **Array**: continuous 2 um barcoded oligos over a 6.5 mm by 6.5 mm capture area. No gaps between capture units (unlike v1 Visium spots).
- **Chemistry**: probe-based. Works on FFPE and fresh frozen. Standard Illumina short-read sequencing downstream.

### How Space Ranger turns raw data into analyzable output

- **Binning**: aggregates 2 um barcodes into square bins.
  - 2 um: maximum resolution, very sparse (<10 UMIs in most bins, rarely used directly).
  - 8 um: roughly 1 to 3 cells per bin. Primary analysis resolution.
  - 16 um: roughly 5 to 15 cells. Robust signal, good for exploration and QC.
- **Segmentation** (new in v4.0): runs StarDist on the H&E image to detect nuclei, then assigns barcodes to individual cells. Output is one row per segmented cell, enabling true single-cell spatial analysis instead of per-bin aggregates.

### Visium v1 vs. Visium HD

| Feature             | Visium v1                      | Visium HD                         |
|---------------------|--------------------------------|-----------------------------------|
| Resolution          | 55 um spots                    | 2 um continuous array             |
| Cells per unit      | ~5 to 30 per spot              | ~1 to 3 per 8 um bin              |
| Gaps                | 100 um center-to-center        | Continuous (none)                 |
| FFPE support        | Yes (v2)                       | Yes (primary)                     |
| Segmentation        | No                             | Yes (StarDist, v4.0+)             |

This slide sets up slide 3, which shows what Space Ranger actually writes to disk.

---

## Slide 3: Space Ranger in Practice (Command, Outputs, File Formats)

**Key message**: Everything spatial you do downstream is the standard 10x MEX matrix plus two extra files that tell you where each barcode lives on the tissue. This slide shows those files literally.

Note: terminal blocks below are representative. File layout, column names, and formats match a real Space Ranger v4.0 run, but the numbers are illustrative. The tutorial itself loads precomputed `.h5ad` files from Figshare, not raw Space Ranger output.

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

One command. It reads the FASTQs, aligns against the probe set, registers the H&E image, and writes a `CRC_P1/outs/` directory that scanpy and squidpy can load directly. Typical runtime: 2 to 8 hours.

### 3b. What `tree CRC_P1/outs/` looks like

```text
CRC_P1/outs/
├── binned_outputs/                              [spatial-specific]
│   ├── square_002um/                            2 um bins  (highest res, sparsest)
│   ├── square_008um/                            8 um bins  (primary analysis)
│   │   ├── filtered_feature_bc_matrix/          standard 10x MEX triplet  (see 3c)
│   │   ├── filtered_feature_bc_matrix.h5        same data, HDF5
│   │   └── spatial/                             [spatial-specific]
│   │       ├── tissue_positions.parquet         THE spatial join key  (see 3d)
│   │       ├── scalefactors_json.json
│   │       ├── tissue_hires_image.png
│   │       └── tissue_lowres_image.png
│   └── square_016um/                            16 um bins (exploration)
├── segmented_outputs/                           [spatial-specific, v4.0+]
│   ├── filtered_feature_cell_matrix.h5
│   └── nucleus_segmentations.geojson            StarDist polygons  (see 3d)
├── cloupe_008um.cloupe                          Loupe browser
├── web_summary.html                             QC report (open in browser)
└── metrics_summary.csv
```

Compared to a plain 10x scRNA-seq run, the only new things are `binned_outputs/` (multiple resolutions), `spatial/tissue_positions.*` (bin to pixel mapping), and `segmented_outputs/*.geojson` (nucleus polygons). Everything else is the standard scRNA-seq file set.

### 3c. MEX format: what `head` actually shows

```text
$ head square_008um/filtered_feature_bc_matrix/features.tsv
ENSG00000187634   SAMD11      Gene Expression
ENSG00000188976   NOC2L       Gene Expression
ENSG00000187961   KLHL17      Gene Expression
ENSG00000187583   HES4        Gene Expression
ENSG00000187642   ISG15       Gene Expression
ENSG00000188157   AGRN        Gene Expression
```

Columns: Ensembl ID, gene symbol, feature type. `scanpy.read_10x_h5` uses column 2 (symbol) as `var_names`. That is why `adata_8um.var_names[:5]` in the notebook returns `['SAMD11', 'NOC2L', 'KLHL17', 'HES4', 'ISG15']`.

```text
$ head square_008um/filtered_feature_bc_matrix/barcodes.tsv
s_008um_00000_00001-1
s_008um_00000_00002-1
s_008um_00000_00003-1
s_008um_00000_00004-1
s_008um_00000_00005-1
```

Barcode format: `s_{binsize}_{row}_{col}-{slide}`. The `(row, col)` integers index the 2 um array grid. This is the hook that lets Space Ranger reassemble spatial coordinates later.

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

Sparse triplet format. The third header line reads `n_genes n_barcodes n_nonzero`. Each data row is `gene_idx barcode_idx umi_count` (1-indexed). `scanpy.read_10x_mtx` parses this into `adata.X` as a CSR sparse matrix with shape `n_barcodes by n_genes`.

None of this is spatial yet. It is the exact same format a 10x Chromium scRNA-seq run writes. Spatial information lives in the next two files.

### 3d. What makes it spatial

**`tissue_positions.parquet`: the bin to pixel map**

```text
$ python -c "import pandas as pd; print(pd.read_parquet('spatial/tissue_positions.parquet').head())"
                 barcode  in_tissue  array_row  array_col  pxl_row_in_fullres  pxl_col_in_fullres
0  s_008um_00000_00001-1          1          0          0              4213.0              2187.0
1  s_008um_00000_00002-1          1          0          1              4213.0              2195.0
2  s_008um_00000_00003-1          1          0          2              4213.0              2203.0
3  s_008um_00000_00004-1          0          0          3              4213.0              2211.0
4  s_008um_00000_00005-1          1          0          4              4213.0              2219.0
```

This is the spatial join key. For every barcode in the MEX matrix it gives:

- `in_tissue`: 1 if the bin sits under the tissue section, 0 if off-tissue. The loader drops rows where this is 0 (see `scripts/utils.py:101-103`).
- `array_row`, `array_col`: grid coordinates on the 2 um array.
- `pxl_row_in_fullres`, `pxl_col_in_fullres`: pixel coordinates on the full-resolution H&E image. This is what `adata.obsm["spatial"]` stores and what every downstream plot uses (`scripts/utils.py:91-93`).

**`nucleus_segmentations.geojson`: the StarDist output**

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

One GeoJSON feature per segmented nucleus, and it is a polygon, not a point. The loader in `scripts/utils.py:164-170` walks `features[]`, averages each polygon's vertices to get a centroid, and uses `(centroid_y, centroid_x)` as the cell's spatial coordinate. Barcodes in the segmented matrix look like `cellid_000000001-1`, and the integer after `cellid_` matches `properties.cell_id`. That is the join.

### 3e. Summary

Everything downstream (scanpy, squidpy, this tutorial) is:

1. Read the MEX matrix into an `AnnData`.
2. Read `tissue_positions.parquet` or the GeoJSON centroids into `adata.obsm["spatial"]`.
3. Filter to `in_tissue == 1`.
4. Run standard analysis, now with coordinates attached.

That is the entire "spatial" part. Two extra files on top of the ordinary 10x output.

---

## Slide 4: When to Use Visium HD

**Key message**: Visium HD is the right call for unbiased whole-transcriptome spatial discovery on FFPE at near-single-cell resolution. It is the wrong call if you already have a gene panel or need live-cell or subcellular measurements.

### Ideal use cases

- **Whole-transcriptome discovery**: no prior gene list, looking for novel spatial markers.
- **FFPE archives**: spatial information from clinical biobanks that cannot go on fresh-frozen platforms.
- **Near-single-cell resolution needed**: 8 um bins contain roughly 1 to 3 cells. Segmentation gets you to one.
- **Heterogeneous tissue**: cell-type composition varies across space (tumor, stroma, immune infiltrate).
- **Integration with histology**: direct H&E overlay, pathologist-friendly.

### Not ideal for

- **Targeted validation**: if you already have a gene panel, Xenium or CosMx are faster and cheaper.
- **Subcellular localization**: use imaging-based methods (MERFISH, seqFISH+).
- **Live-cell dynamics**: fixed tissue only, no time course.
- **3D reconstruction**: 2D sections. Serial-section 3D is possible but complex.
- **Large cohorts (>50 samples)**: cost and compute scale become the dominant constraint.

### Practical numbers

| Axis                 | Typical                                                 |
|----------------------|---------------------------------------------------------|
| Data size            | ~5 to 50 GB per sample (depth and bin-size dependent)   |
| RAM for analysis     | 16+ GB for 8 um, 64+ GB for 2 um                        |
| Space Ranger runtime | 2 to 8 hours                                            |
| Downstream runtime   | Hours to days                                           |
| Cost per sample      | ~$1,500 to $2,500 (reagents, CytAssist, sequencing)     |

---

## Slide 5: Analysis Pipeline and QC Considerations

**Key message**: Visium HD pipelines are standard scanpy plus spatial-aware QC and spatial-aware analysis steps. A few gotchas are specific to the probe panel and the binning strategy.

### Pipeline overview

1. **QC and filtering**: spatial-aware, not just per-bin thresholds.
2. **Normalization and HVG selection**: standard scanpy.
3. **Dimensionality reduction**: PCA and UMAP for visualization.
4. **Cell-type annotation**: deconvolution (binned) or label transfer and marker scoring (segmented).
5. **Spatial analysis**: neighborhood enrichment, spatially variable genes (Moran's I).
6. **Cell-cell communication**: ligand-receptor analysis in spatial context.
7. **Validation**: compare with ground-truth annotations, orthogonal data, known biology.

### Choosing a resolution

- Start with 16 um for exploration and QC (best signal-to-noise).
- Use 8 um for primary analysis (best balance, and what this tutorial uses).
- Use 2 um only for specific high-resolution questions. It is very sparse.
- Use segmented outputs when you need true single-cell analysis.

### QC considerations unique to Visium HD

- **Spatial QC is essential.** Tissue folds, debris, and edge effects create spatially-correlated artifacts that per-bin thresholds miss. Use Moran's I on QC metrics and spatial neighbor-based outlier detection (see `spatial_outlier_detection` in `scripts/utils.py`).
- **High-mt regions are not automatically bad.** If high `pct_counts_mt` clusters at tissue edges or folds, it is damage, and you should filter. If it is dispersed across the tumor, it is likely real hypoxia or necrosis, and you should keep it.
- **The 10x Human WTA probe panel excludes ribosomal protein genes (RPL\*, RPS\*) to prevent probe saturation.** This means `pct_counts_ribo` will be flat zero on this data. It is expected, not a bug. Use `pct_counts_mt` as the primary metabolic QC axis and do not spend time debugging the empty ribo panel. `compute_qc_metrics` in `scripts/utils.py:204` now logs this explicitly.
- **Resolution-dependent thresholds.** QC cutoffs differ across 2, 8, and 16 um bins. Do not reuse numbers across resolutions.
- **Segmentation artifacts.** Over-segmentation, under-segmentation, nuclear vs. cytoplasmic capture. Sanity-check cell counts against the H&E.

### Deconvolution vs. segmentation

- **Binned data (8 or 16 um)**: multiple cells per bin, so use deconvolution (cell2location, RCTD, STdeconvolve).
- **Segmented data**: one cell per observation, so use direct annotation (label transfer, marker scoring).
- **Best practice**: run both where possible and check for concordance.

### What this tutorial demonstrates

Real QC on a clinically relevant human CRC dataset, ground-truth cell-type comparison, neighborhood enrichment, spatially variable genes, and ligand-receptor interactions, using the file formats shown in slide 3.
