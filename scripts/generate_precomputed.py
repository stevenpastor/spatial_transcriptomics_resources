"""
Generate pre-processed .h5ad files for Colab / Binder distribution.

Run this ONCE on a machine where the full CRC P1 dataset has already been
processed (i.e. the checkpoint files exist in data/checkpoints_crc/).

Usage:
    python scripts/generate_precomputed.py

Outputs:
    data/precomputed/crc_p1_raw_50k.h5ad        (< 100 MB)
    data/precomputed/crc_p1_annotated_50k.h5ad   (< 100 MB)
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
    # Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"SUCCESS")
    print(f"{'='*60}")
    print(f"  {raw_path.name:40s}  {file_mb(raw_path):6.1f} MB")
    print(f"  {annot_path.name:40s}  {file_mb(annot_path):6.1f} MB")
    if SPATIAL_OUT.exists():
        for f in sorted(SPATIAL_OUT.iterdir()):
            print(f"  spatial/{f.name:36s}  {file_mb(f):6.1f} MB")
    return True


if __name__ == "__main__":
    ok = generate(TARGET_BINS)
    if not ok:
        print(f"\nFiles too large at {TARGET_BINS:,} bins. Retrying with {FALLBACK_BINS:,}...")
        ok = generate(FALLBACK_BINS)
        if not ok:
            print("ERROR: Files still too large. Manual intervention needed.")
            sys.exit(1)
