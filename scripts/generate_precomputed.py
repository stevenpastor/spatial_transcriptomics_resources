"""
Generate pre-processed .h5ad files for Colab / Binder distribution.

Run this ONCE on a machine where the full CRC P1 dataset has already been
processed (i.e. the checkpoint files exist in data/checkpoints_crc/).

Usage:
    python scripts/generate_precomputed.py

Outputs:
    data/precomputed/crc_p1_raw_50k.h5ad         (< 100 MB, 50K subsample)
    data/precomputed/crc_p1_annotated_50k.h5ad   (< 100 MB, 50K subsample)
    data/precomputed/crc_p1_liana.parquet        (~0.1 - few MB, FULL tissue)

The h5ads are 50K-bin subsamples so the notebook can download them
quickly on Colab. The LIANA parquet, however, is computed on the full
~500K-bin annotated checkpoint so the L-R statistics reflect the
whole tissue and small populations (SELENOP+ macrophages, sparse
chemokines) clear LIANA's expression-proportion thresholds.
"""

import sys
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
TARGET_BINS = 50_000          # subsample target
MAX_FILE_MB = 95              # abort if larger (GitHub limit is 100 MB)
FALLBACK_BINS = 40_000        # if 50K is too big, retry with this
RANDOM_STATE = 42

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
CHECKPOINT_DIR = DATA_DIR / "checkpoints_crc"
OUTPUT_DIR = DATA_DIR / "precomputed"
SPATIAL_OUT = OUTPUT_DIR / "spatial"

METADATA_PATH = DATA_DIR / "P1CRC_Metadata.parquet"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def file_mb(path):
    return path.stat().st_size / 1e6


def subsample_adata(adata, n_obs, seed=RANDOM_STATE):
    """Subsample to n_obs, preserving sparse format."""
    if adata.shape[0] <= n_obs:
        print(f"  Already <= {n_obs:,} bins, no subsampling needed.")
        return adata.copy()
    np.random.seed(seed)
    idx = np.random.choice(adata.shape[0], size=n_obs, replace=False)
    idx.sort()
    return adata[idx].copy()


def filter_zero_genes(adata):
    """Remove genes with zero counts across all bins."""
    if sparse.issparse(adata.X):
        gene_counts = np.asarray(adata.X.sum(axis=0)).flatten()
    else:
        gene_counts = adata.X.sum(axis=0)
    keep = gene_counts > 0
    n_removed = (~keep).sum()
    adata = adata[:, keep].copy()
    print(f"  Removed {n_removed:,} zero-count genes -> {adata.shape[1]:,} genes remain")
    return adata


def embed_annotations(adata, meta_path):
    """Attach ground-truth annotations from the metadata parquet."""
    if not meta_path.exists():
        print(f"  WARNING: Metadata not found at {meta_path}, skipping annotations.")
        return adata

    meta = pd.read_parquet(meta_path)
    meta["tissue"] = meta["tissue"].astype(str)
    meta = meta[meta["tissue"] == "1"].copy()
    meta = meta.set_index("barcode")

    common = adata.obs_names.intersection(meta.index)
    print(f"  Matched {len(common):,} / {adata.shape[0]:,} bins to metadata")

    for col in ["DeconvolutionClass", "DeconvolutionLabel1", "DeconvolutionLabel2",
                "UnsupervisedL1", "UnsupervisedL2", "Periphery"]:
        if col in meta.columns:
            series = meta.reindex(adata.obs_names)[col]
            if hasattr(series, "cat"):
                series = series.astype(object)
            adata.obs[col] = series.values

    n_labeled = adata.obs["UnsupervisedL1"].notna().sum()
    print(f"  UnsupervisedL1 coverage: {n_labeled:,} / {adata.shape[0]:,}")
    return adata


def embed_spatial_image(adata, outs_dir):
    """Find and embed the H&E image + scalefactors in adata.uns."""
    for search_dir in [outs_dir / "spatial",
                       outs_dir / "binned_outputs" / "square_008um" / "spatial"]:
        img_path = search_dir / "tissue_lowres_image.png"
        sf_path = search_dir / "scalefactors_json.json"
        if img_path.exists() and sf_path.exists():
            from skimage import io as skio
            img = skio.imread(str(img_path))
            with open(sf_path) as f:
                sf = json.load(f)
            adata.uns["spatial_image"] = img
            adata.uns["scalefactors"] = sf
            print(f"  Embedded H&E image ({img.shape}) and scalefactors")

            # Also copy to output spatial dir
            SPATIAL_OUT.mkdir(parents=True, exist_ok=True)
            shutil.copy2(img_path, SPATIAL_OUT / "tissue_lowres_image.png")
            shutil.copy2(sf_path, SPATIAL_OUT / "scalefactors_json.json")
            # Also grab hires if available
            hires = search_dir / "tissue_hires_image.png"
            if hires.exists():
                shutil.copy2(hires, SPATIAL_OUT / "tissue_hires_image.png")
            print(f"  Copied spatial files to {SPATIAL_OUT}")
            return adata

    print("  WARNING: H&E image not found, skipping embedding.")
    return adata


def save_with_size_check(adata, path, max_mb=MAX_FILE_MB):
    """Save h5ad and return True if under size limit."""
    adata.write_h5ad(path)
    size = file_mb(path)
    print(f"  Saved: {path.name} ({size:.1f} MB)")
    if size > max_mb:
        print(f"  WARNING: {size:.1f} MB exceeds {max_mb} MB limit!")
        return False
    return True


def subcluster_macrophages(adata, mac_labels=("Macrophage", "Proliferating Macrophages"),
                           n_top_genes=2000, n_pcs=30, n_neighbors=15,
                           resolution=0.8, min_subtype_size=50):
    """Re-cluster the macrophage subset and label SPP1+ vs SELENOP+ subclusters.

    Returns a Series (indexed by adata.obs_names) labeling each macrophage bin
    as 'Macrophage_SPP1+', 'Macrophage_SELENOP+', or 'Macrophage_other'.
    Non-macrophage bins are mapped to None.

    This recovers the paper's headline cell-cell-communication finding (Fig 5):
    SPP1+ macrophages at the tumor edge driving ECM/integrin signaling vs.
    SELENOP+ macrophages elsewhere in the TME.
    """
    base = adata.obs["DeconvolutionLabel1"].astype(str)
    mac_mask = base.isin(mac_labels).values
    print(f"  Macrophage bins to subcluster: {int(mac_mask.sum()):,}")

    adata_mac = adata[mac_mask].copy()
    if "counts" in adata_mac.layers:
        adata_mac.X = adata_mac.layers["counts"].copy()
    sc.pp.normalize_total(adata_mac, target_sum=1e4)
    sc.pp.log1p(adata_mac)
    sc.pp.highly_variable_genes(adata_mac, n_top_genes=n_top_genes, flavor="seurat")
    adata_mac_h = adata_mac[:, adata_mac.var["highly_variable"]].copy()
    sc.pp.scale(adata_mac_h, max_value=10)
    sc.pp.pca(adata_mac_h, n_comps=n_pcs)
    sc.pp.neighbors(adata_mac_h, n_neighbors=n_neighbors)
    sc.tl.leiden(adata_mac_h, resolution=resolution, key_added="leiden_mac")

    def expr(a, gene):
        X = a[:, gene].X
        return (X.toarray() if sparse.issparse(X) else np.asarray(X)).flatten()

    spp1 = expr(adata_mac, "SPP1")
    selp = expr(adata_mac, "SELENOP")
    clusters = adata_mac_h.obs["leiden_mac"].astype(int).values

    candidates = []
    for c in sorted(np.unique(clusters)):
        m = (clusters == c)
        if m.sum() < min_subtype_size:
            continue
        candidates.append((c, spp1[m].mean(), selp[m].mean(), int(m.sum())))

    if len(candidates) < 2:
        print(f"  WARNING: fewer than 2 candidate macrophage subclusters >= "
              f"{min_subtype_size} bins; subtype labels NOT created.")
        return pd.Series(index=adata.obs_names, dtype=object)

    spp1_top = max(candidates, key=lambda s: s[1])
    selp_top = max(candidates, key=lambda s: s[2])
    if spp1_top[0] == selp_top[0]:
        print("  WARNING: same cluster wins both SPP1+ and SELENOP+; "
              "subtype labels NOT created.")
        return pd.Series(index=adata.obs_names, dtype=object)
    print(f"  SPP1+ subcluster {spp1_top[0]}: n={spp1_top[3]}, "
          f"mean SPP1={spp1_top[1]:.2f}, mean SELENOP={spp1_top[2]:.2f}")
    print(f"  SELENOP+ subcluster {selp_top[0]}: n={selp_top[3]}, "
          f"mean SPP1={selp_top[1]:.2f}, mean SELENOP={selp_top[2]:.2f}")

    mac_label = np.array(["Macrophage_other"] * len(clusters), dtype=object)
    mac_label[clusters == spp1_top[0]] = "Macrophage_SPP1+"
    mac_label[clusters == selp_top[0]] = "Macrophage_SELENOP+"

    result = pd.Series(index=adata.obs_names, dtype=object)
    result.iloc[np.where(mac_mask)[0]] = mac_label
    return result


def run_liana(adata, groupby, min_group_size=50,
              expr_prop=0.1, min_cells=20, n_perms=1000, seed=RANDOM_STATE):
    """Run liana.mt.rank_aggregate and return the result DataFrame."""
    import liana as li

    obs = adata.obs[groupby].astype(str)
    keep_mask = (obs != "nan") & obs.notna()
    sub = adata[keep_mask.values].copy()

    counts = sub.obs[groupby].value_counts()
    keep_types = counts[counts >= min_group_size].index.tolist()
    sub = sub[sub.obs[groupby].isin(keep_types)].copy()
    sub.obs[groupby] = sub.obs[groupby].astype("category")

    print(f"  LIANA on {sub.n_obs:,} bins x {sub.n_vars:,} genes, "
          f"{len(keep_types)} cell types (>= {min_group_size} bins each)")

    li.mt.rank_aggregate(
        sub,
        groupby=groupby,
        resource_name="consensus",
        expr_prop=expr_prop,
        min_cells=min_cells,
        use_raw=False,
        verbose=False,
        n_perms=n_perms,
        seed=seed,
    )
    return sub.uns["liana_res"].copy(), sub.obs[groupby].copy()


def generate_liana(annot_path, out_path):
    """Compute LIANA + macrophage subtypes on the FULL annotated checkpoint.

    `annot_path` is kept in the signature for caller compatibility but is
    not read. LIANA always runs on the full, unsubsampled checkpoint at
    `CHECKPOINT_DIR / "crc_p1_8um_annotated.h5ad"` so the shipped parquet
    reflects whole-tissue statistics rather than the 50K subsample used
    for in-notebook visualization.

    Adds SPP1+ vs SELENOP+ macrophage labels via subclustering on the
    macrophage subset, then runs LIANA at that granularity. Writes the
    L-R results plus the macrophage subtype labels (as an embedded side
    table) into the parquet schema metadata, so the notebook can color
    the spatial overlay without recomputing.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    full_ckpt = CHECKPOINT_DIR / "crc_p1_8um_annotated.h5ad"
    if not full_ckpt.exists():
        raise FileNotFoundError(
            f"Full annotated checkpoint missing: {full_ckpt}. "
            "LIANA must run on the full ~500K-bin tissue, not the 50K "
            "subsample. Build this checkpoint by running the notebook's "
            "QC + normalization + annotation cells on the full Space Ranger "
            "output with subsampling disabled."
        )

    print(f"\n--- LIANA precompute (full tissue) ---")
    adata = sc.read_h5ad(full_ckpt)
    print(f"  Loaded annotated (full): {adata.shape}")

    if "DeconvolutionLabel1" not in adata.obs.columns:
        adata = embed_annotations(adata, METADATA_PATH)
    adata = filter_zero_genes(adata)
    if adata.raw is not None:
        adata.raw = None

    mac_subtypes = subcluster_macrophages(adata)

    cell_type = adata.obs["DeconvolutionLabel1"].astype(str).copy()
    has_subtype = mac_subtypes.notna()
    cell_type.loc[has_subtype] = mac_subtypes.loc[has_subtype].astype(str)
    adata.obs["cell_type_liana"] = cell_type.values

    res, used_labels = run_liana(adata, groupby="cell_type_liana")
    print(f"  LIANA rows: {len(res):,}")

    # Embed macrophage subtype labels so the notebook can re-color the overlay
    mac_lbl_df = (used_labels[used_labels.astype(str).str.startswith("Macrophage")]
                  .rename("subtype")
                  .reset_index()
                  .rename(columns={"index": "barcode"}))
    mac_lbl_df["barcode"] = mac_lbl_df["barcode"].astype(str)
    mac_lbl_df["subtype"] = mac_lbl_df["subtype"].astype(str)

    table = pa.Table.from_pandas(res, preserve_index=False)
    mac_buf = pa.BufferOutputStream()
    pq.write_table(pa.Table.from_pandas(mac_lbl_df, preserve_index=False),
                   mac_buf, compression="zstd")
    mac_blob = mac_buf.getvalue().to_pybytes()

    meta = {
        b"liana_groupby": b"cell_type_liana (DeconvolutionLabel1 + macrophage subclustering)",
        b"liana_resource": b"consensus",
        b"liana_n_perms": b"1000",
        b"liana_n_bins": str(int(adata.n_obs)).encode(),
        b"liana_subset": b"full_tissue",
        b"liana_dataset": b"CRC P1 (de Oliveira et al. 2025, Nat Genet)",
        b"liana_mac_subtype_table_parquet": mac_blob,
    }
    existing = dict(table.schema.metadata or {})
    table = table.replace_schema_metadata({**existing, **meta})
    pq.write_table(table, out_path)

    size_mb = file_mb(out_path)
    print(f"  Saved: {out_path.name} ({size_mb:.2f} MB)")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def generate(n_bins=TARGET_BINS):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Raw checkpoint
    # ------------------------------------------------------------------
    raw_ckpt = CHECKPOINT_DIR / "crc_p1_8um_raw.h5ad"
    if not raw_ckpt.exists():
        print(f"ERROR: Raw checkpoint not found: {raw_ckpt}")
        print("Run the notebook first to generate checkpoints.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Generating precomputed data ({n_bins:,} bins)")
    print(f"{'='*60}")

    print(f"\n--- Raw checkpoint ---")
    adata_raw = sc.read_h5ad(raw_ckpt)
    print(f"  Loaded: {adata_raw.shape[0]:,} x {adata_raw.shape[1]:,}")

    adata_raw = subsample_adata(adata_raw, n_bins)
    adata_raw = filter_zero_genes(adata_raw)
    adata_raw = embed_annotations(adata_raw, METADATA_PATH)

    # Find outs_dir for spatial images
    dataset_dir = DATA_DIR / "Visium_HD_Human_Colon_Cancer_P1"
    for candidate in [dataset_dir, dataset_dir / "outs"]:
        if (candidate / "spatial").exists() or (candidate / "binned_outputs").exists():
            adata_raw = embed_spatial_image(adata_raw, candidate)
            break

    raw_path = OUTPUT_DIR / "crc_p1_raw_50k.h5ad"
    if not save_with_size_check(adata_raw, raw_path):
        return False

    # ------------------------------------------------------------------
    # 2. Annotated checkpoint
    # ------------------------------------------------------------------
    annot_ckpt = CHECKPOINT_DIR / "crc_p1_8um_annotated.h5ad"
    if not annot_ckpt.exists():
        # Fall back to normalized checkpoint
        annot_ckpt = CHECKPOINT_DIR / "crc_p1_8um_normalized.h5ad"
    if not annot_ckpt.exists():
        print(f"ERROR: No annotated/normalized checkpoint found in {CHECKPOINT_DIR}")
        sys.exit(1)

    print(f"\n--- Annotated checkpoint ---")
    adata_annot = sc.read_h5ad(annot_ckpt)
    print(f"  Loaded: {adata_annot.shape[0]:,} x {adata_annot.shape[1]:,}")

    adata_annot = subsample_adata(adata_annot, n_bins)
    adata_annot = filter_zero_genes(adata_annot)

    # Ensure annotations are present
    if "UnsupervisedL1" not in adata_annot.obs.columns:
        adata_annot = embed_annotations(adata_annot, METADATA_PATH)

    # Drop .raw to save space; keep raw counts in layers["counts"]
    if adata_annot.raw is not None:
        if "counts" not in adata_annot.layers:
            # Store raw counts from .raw
            raw_X = adata_annot.raw[:, adata_annot.var_names].X
            adata_annot.layers["counts"] = raw_X
        adata_annot.raw = None
        print("  Dropped .raw, raw counts preserved in layers['counts']")

    # Embed spatial image
    if "spatial_image" not in adata_annot.uns:
        for candidate in [dataset_dir, dataset_dir / "outs"]:
            if (candidate / "spatial").exists() or (candidate / "binned_outputs").exists():
                adata_annot = embed_spatial_image(adata_annot, candidate)
                break

    annot_path = OUTPUT_DIR / "crc_p1_annotated_50k.h5ad"
    if not save_with_size_check(adata_annot, annot_path):
        return False

    # ------------------------------------------------------------------
    # 3. LIANA cell-cell communication artifact
    # ------------------------------------------------------------------
    liana_path = OUTPUT_DIR / "crc_p1_liana.parquet"
    try:
        generate_liana(annot_path, liana_path)
    except ImportError as exc:
        print(f"\n  WARNING: LIANA precompute skipped ({exc}). "
              "Install liana to enable: pip install liana")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"SUCCESS")
    print(f"{'='*60}")
    print(f"  {raw_path.name:40s}  {file_mb(raw_path):6.1f} MB")
    print(f"  {annot_path.name:40s}  {file_mb(annot_path):6.1f} MB")
    if liana_path.exists():
        print(f"  {liana_path.name:40s}  {file_mb(liana_path):6.2f} MB")
    if SPATIAL_OUT.exists():
        for f in sorted(SPATIAL_OUT.iterdir()):
            print(f"  spatial/{f.name:36s}  {file_mb(f):6.1f} MB")
    print(f"\nNext: upload all files in {OUTPUT_DIR} to Figshare article 31937262.")
    return True


if __name__ == "__main__":
    ok = generate(TARGET_BINS)
    if not ok:
        print(f"\nFiles too large at {TARGET_BINS:,} bins. Retrying with {FALLBACK_BINS:,}...")
        ok = generate(FALLBACK_BINS)
        if not ok:
            print("ERROR: Files still too large. Manual intervention needed.")
            sys.exit(1)
