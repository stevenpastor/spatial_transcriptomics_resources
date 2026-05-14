#!/usr/bin/env python3
"""
Precompute LIANA L-R results for the FULL Visium HD CRC P1 tissue.

Self-contained driver for SLURM / headless servers. Goes from a raw 10x
Space Ranger output directory + Oliveira's per-barcode metadata parquet to
the final `crc_p1_liana.parquet` in one job. Skips intermediate h5ad
checkpoints; only the parquet is written.

Place this file next to `generate_precomputed.py` (it imports helpers from
there).

Usage:
    python precompute_liana_full.py \\
        --space-ranger-dir /path/to/Visium_HD_Human_Colon_Cancer_P1 \\
        --metadata-parquet /path/to/P1CRC_Metadata.parquet \\
        --output-parquet  /path/to/crc_p1_liana.parquet

Resource hint: ~30-60 GB RAM, ~1-3 hours wall time at n_perms=1000.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from generate_precomputed import (
    embed_annotations, filter_zero_genes,
    subcluster_macrophages, run_liana,
)


def find_bin8_dir(space_ranger_dir: Path) -> Path:
    """Locate the square_008um output inside a 10x Visium HD bundle."""
    candidates = [
        space_ranger_dir / "binned_outputs" / "square_008um",
        space_ranger_dir / "outs" / "binned_outputs" / "square_008um",
        space_ranger_dir,
    ]
    for c in candidates:
        if (c / "filtered_feature_bc_matrix.h5").exists():
            return c
    raise FileNotFoundError(
        f"No filtered_feature_bc_matrix.h5 found under {space_ranger_dir}. "
        f"Checked: {[str(c) for c in candidates]}"
    )


def load_visium_hd_8um(bin8_dir: Path):
    """Load counts + spatial coords + in_tissue flag for the 8 um bins."""
    h5_path = bin8_dir / "filtered_feature_bc_matrix.h5"
    spatial_dir = bin8_dir / "spatial"
    pos_path = spatial_dir / "tissue_positions.parquet"
    if not pos_path.exists():
        pos_path = spatial_dir / "tissue_positions.csv"

    print(f"Loading counts from {h5_path}")
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    print(f"Loading spatial positions from {pos_path}")
    if pos_path.suffix == ".parquet":
        positions = pd.read_parquet(pos_path)
    else:
        positions = pd.read_csv(pos_path)
    positions = positions.set_index("barcode")
    positions = positions.reindex(adata.obs_names)

    if "in_tissue" in positions.columns:
        adata.obs["in_tissue"] = positions["in_tissue"].values
    coord_cols = ("pxl_col_in_fullres", "pxl_row_in_fullres")
    adata.obsm["spatial"] = positions[list(coord_cols)].to_numpy()
    return adata


def write_liana_parquet(res, used_labels, n_bins, out_path: Path):
    import pyarrow as pa
    import pyarrow.parquet as pq

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
        b"liana_n_bins": str(int(n_bins)).encode(),
        b"liana_subset": b"full_tissue",
        b"liana_dataset": b"CRC P1 (de Oliveira et al. 2025, Nat Genet)",
        b"liana_mac_subtype_table_parquet": mac_blob,
    }
    existing = dict(table.schema.metadata or {})
    table = table.replace_schema_metadata({**existing, **meta})

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, out_path)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--space-ranger-dir", type=Path, required=True,
                   help="Top-level 10x Visium HD Space Ranger output for P1 "
                        "(the dir containing `binned_outputs/` or `outs/`).")
    p.add_argument("--metadata-parquet", type=Path, required=True,
                   help="Path to P1CRC_Metadata.parquet from Oliveira et al.")
    p.add_argument("--output-parquet", type=Path, required=True,
                   help="Where to write crc_p1_liana.parquet.")
    args = p.parse_args()

    bin8_dir = find_bin8_dir(args.space_ranger_dir)
    print(f"Using 8 um bins from: {bin8_dir}")

    adata = load_visium_hd_8um(bin8_dir)
    print(f"Loaded: {adata.shape}")

    if "in_tissue" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["in_tissue"] == 1].copy()
        print(f"in_tissue filter: {before:,} -> {adata.n_obs:,}")

    adata = embed_annotations(adata, args.metadata_parquet)
    if "DeconvolutionLabel1" not in adata.obs.columns:
        raise RuntimeError(
            "DeconvolutionLabel1 missing after embed_annotations. "
            "Check that the metadata parquet's barcodes match this Space "
            "Ranger output (it must be the 8 um bin layout for P1)."
        )

    has_label = adata.obs["DeconvolutionLabel1"].notna()
    print(f"Label coverage: {has_label.sum():,} / {adata.n_obs:,}")
    adata = adata[has_label].copy()

    adata = filter_zero_genes(adata)

    print("Normalizing (normalize_total + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("Subclustering macrophages -> SPP1+ / SELENOP+ / other...")
    mac_subtypes = subcluster_macrophages(adata)
    cell_type = adata.obs["DeconvolutionLabel1"].astype(str).copy()
    has_subtype = mac_subtypes.notna()
    cell_type.loc[has_subtype] = mac_subtypes.loc[has_subtype].astype(str)
    adata.obs["cell_type_liana"] = cell_type.values

    print("Running LIANA rank_aggregate (this is the slow step)...")
    res, used_labels = run_liana(adata, groupby="cell_type_liana")
    print(f"LIANA rows: {len(res):,}")

    write_liana_parquet(res, used_labels, adata.n_obs, args.output_parquet)
    size_mb = args.output_parquet.stat().st_size / 1e6
    print(f"\nDone. Wrote {args.output_parquet} ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()
