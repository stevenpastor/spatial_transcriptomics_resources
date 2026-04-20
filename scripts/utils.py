"""
Utility functions for the Visium HD Spatial Transcriptomics Tutorial.

Functions for data loading, spatial QC, cell type comparison, and visualization.
"""

import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy import sparse
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------------
# Data Loading Helpers
# ---------------------------------------------------------------------------

def load_visium_hd_binned(data_dir: str, bin_size: str = "008um") -> sc.AnnData:
    """Load Visium HD binned output as AnnData.

    Parameters
    ----------
    data_dir : str
        Path to the Space Ranger output directory (the 'outs' folder or
        the parent that contains 'binned_outputs').
    bin_size : str
        Bin size folder name, e.g. '008um', '016um', '002um'.

    Returns
    -------
    AnnData with spatial coordinates in .obsm['spatial'].
    """
    data_path = Path(data_dir)

    # Locate the binned output directory
    candidates = [
        data_path / "binned_outputs" / f"square_{bin_size}",
        data_path / f"square_{bin_size}",
    ]
    bin_dir = None
    for c in candidates:
        if c.exists():
            bin_dir = c
            break
    if bin_dir is None:
        raise FileNotFoundError(
            f"Cannot find binned output for {bin_size} in {data_path}. "
            f"Tried: {[str(c) for c in candidates]}"
        )

    # Load the filtered feature-barcode matrix
    matrix_dir = bin_dir / "filtered_feature_bc_matrix.h5"
    if matrix_dir.exists():
        adata = sc.read_10x_h5(str(matrix_dir))
    else:
        # Try the directory-based format
        matrix_dir = bin_dir / "filtered_feature_bc_matrix"
        adata = sc.read_10x_mtx(str(matrix_dir))

    adata.var_names_make_unique()

    # Load spatial coordinates from tissue_positions.parquet or .csv
    positions_parquet = bin_dir / "spatial" / "tissue_positions.parquet"
    positions_csv = bin_dir / "spatial" / "tissue_positions.csv"

    if positions_parquet.exists():
        positions = pd.read_parquet(positions_parquet)
    elif positions_csv.exists():
        positions = pd.read_csv(positions_csv, header=0,
                                index_col=0)
    else:
        raise FileNotFoundError(
            f"Cannot find tissue_positions in {bin_dir / 'spatial'}"
        )

    # Align positions to adata
    if "barcode" in positions.columns:
        positions = positions.set_index("barcode")

    common = adata.obs_names.intersection(positions.index)
    adata = adata[common].copy()
    positions = positions.loc[common]

    # Store spatial coordinates
    # Typical columns: in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
    coord_cols = [c for c in positions.columns if "pxl" in c.lower()]
    if len(coord_cols) >= 2:
        adata.obsm["spatial"] = positions[coord_cols].values.astype(np.float64)
    else:
        # Fallback to array row/col
        row_col = [c for c in positions.columns if "row" in c.lower() or "col" in c.lower()]
        if len(row_col) >= 2:
            adata.obsm["spatial"] = positions[row_col[-2:]].values.astype(np.float64)

    # Store tissue flag
    if "in_tissue" in positions.columns:
        adata.obs["in_tissue"] = positions["in_tissue"].values
        adata = adata[adata.obs["in_tissue"] == 1].copy()

    adata.obs["bin_size"] = bin_size
    return adata


def load_visium_hd_segmented(data_dir: str) -> sc.AnnData:
    """Load Visium HD segmented (single-cell) output.

    Space Ranger v4.0+ produces segmented outputs via StarDist nuclei
    segmentation. This function loads those results. Spatial coordinates
    are derived from nucleus segmentation polygon centroids in the GeoJSON.

    Parameters
    ----------
    data_dir : str
        Path to the Space Ranger output directory (containing
        'segmented_outputs/' subfolder).

    Returns
    -------
    AnnData with spatial coordinates in .obsm['spatial'].
    """
    import json

    data_path = Path(data_dir)

    # Locate segmented output directory
    seg_dir = None
    for candidate in [data_path / "segmented_outputs", data_path]:
        if (candidate / "filtered_feature_cell_matrix.h5").exists():
            seg_dir = candidate
            break
        if (candidate / "filtered_feature_cell_matrix").exists():
            seg_dir = candidate
            break

    if seg_dir is None:
        raise FileNotFoundError(
            f"Cannot find segmented outputs in {data_path}. "
            "Expected 'segmented_outputs/filtered_feature_cell_matrix.h5'."
        )

    # Load the filtered cell matrix (note: *cell* not *bc* matrix)
    h5_path = seg_dir / "filtered_feature_cell_matrix.h5"
    if h5_path.exists():
        adata = sc.read_10x_h5(str(h5_path))
    else:
        adata = sc.read_10x_mtx(str(seg_dir / "filtered_feature_cell_matrix"))

    adata.var_names_make_unique()

    # Extract spatial coordinates from nucleus segmentation GeoJSON
    # Barcodes are like 'cellid_000000001-1' → cell_id = 1
    geojson_path = seg_dir / "nucleus_segmentations.geojson"
    if geojson_path.exists():
        print(f"  Loading nucleus segmentations from GeoJSON...")
        with open(geojson_path) as f:
            geojson = json.load(f)

        # Compute polygon centroids for each cell_id
        centroids = {}
        for feature in geojson["features"]:
            cell_id = feature["properties"]["cell_id"]
            coords = np.array(feature["geometry"]["coordinates"][0])
            centroid_x = coords[:, 0].mean()
            centroid_y = coords[:, 1].mean()
            centroids[cell_id] = (centroid_y, centroid_x)  # row, col convention

        # Map barcodes to cell_ids and assign coordinates
        # Barcode format: cellid_XXXXXXXXX-1
        cell_ids = []
        for bc in adata.obs_names:
            # Extract numeric cell_id from barcode
            try:
                cid = int(bc.split("_")[1].split("-")[0])
            except (IndexError, ValueError):
                cid = None
            cell_ids.append(cid)

        spatial_coords = np.zeros((adata.n_obs, 2), dtype=np.float64)
        matched = 0
        for i, cid in enumerate(cell_ids):
            if cid is not None and cid in centroids:
                spatial_coords[i] = centroids[cid]
                matched += 1

        adata.obsm["spatial"] = spatial_coords
        print(f"  Matched {matched}/{adata.n_obs} cells to spatial coordinates")
    else:
        print("  Warning: No nucleus_segmentations.geojson found. "
              "Spatial coordinates not available for segmented data.")

    adata.obs["data_type"] = "segmented"
    return adata


# ---------------------------------------------------------------------------
# QC Helpers
# ---------------------------------------------------------------------------

def compute_qc_metrics(adata: sc.AnnData, verbose: bool = True) -> sc.AnnData:
    """Compute comprehensive QC metrics for Visium HD data.

    Adds: total_counts, n_genes_by_counts, pct_counts_mt, pct_counts_ribo,
          log1p_total_counts, log1p_n_genes, complexity.

    Robust to var_names being Ensembl IDs: looks for a gene-symbol column
    in .var (gene_symbols / feature_name / symbol / gene_name) and matches
    against that. Prints a diagnostic of how many genes were flagged so you
    can catch silently-empty matches.
    """
    # Find the most plausible source of gene symbols. read_10x_h5 / read_10x_mtx
    # usually puts symbols in var_names, but with some Visium HD outputs the
    # var_names end up as Ensembl IDs and the symbols live in a var column.
    symbol_cols = ["gene_symbols", "feature_name", "symbol", "gene_name", "Symbol"]
    symbol_source = None
    for col in symbol_cols:
        if col in adata.var.columns:
            symbol_source = col
            break

    if symbol_source is not None:
        symbols = adata.var[symbol_source].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)

    # If var_names look like Ensembl IDs but we somehow didn't find a symbol
    # column, fall back to var_names anyway (matching will just miss).
    looks_ensembl = symbols.str.match(r"^ENS[A-Z]*\d{6,}").mean() > 0.5
    if looks_ensembl and symbol_source is None and verbose:
        print("  [warn] var_names look like Ensembl IDs and no gene-symbol "
              "column was found in .var, mt/ribo matching will likely miss.")

    # Case-insensitive matching catches MT-/mt-, RPS/Rps, RPL/Rpl in one pass.
    mt_mask = symbols.str.match(r"^mt-", case=False, na=False)
    ribo_mask = symbols.str.match(r"^rp[sl]\d", case=False, na=False)

    adata.var["mt"] = mt_mask.values
    adata.var["ribo"] = ribo_mask.values

    n_mt = int(mt_mask.sum())
    n_ribo = int(ribo_mask.sum())

    if verbose:
        src = f".var['{symbol_source}']" if symbol_source else ".var_names"
        print(f"  QC gene matching against {src}: {n_mt} mt, {n_ribo} ribo")
        if n_ribo == 0:
            print("  [info] I removed RPL*/RPS* genes if you are looking for them later on your own.")

    qc_vars = ["mt", "ribo"] if n_ribo > 0 else ["mt"]
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=qc_vars, percent_top=None, log1p=True, inplace=True
    )
    # Keep the column present (as zeros) so downstream code referencing it
    # doesn't KeyError, but it will be uninformative.
    if n_ribo == 0 and "pct_counts_ribo" not in adata.obs.columns:
        adata.obs["pct_counts_ribo"] = 0.0
    adata.uns["qc_has_ribo"] = bool(n_ribo > 0)

    # Complexity: log(genes) / log(counts), measures transcriptomic diversity
    adata.obs["complexity"] = (
        np.log1p(adata.obs["n_genes_by_counts"])
        / np.log1p(adata.obs["total_counts"])
    )
    # Handle edge cases
    adata.obs["complexity"] = adata.obs["complexity"].replace(
        [np.inf, -np.inf], np.nan
    ).fillna(0)

    return adata


def spatial_qc_plots(
    adata: sc.AnnData,
    metrics: list[str] = None,
    figsize: tuple = (16, 4),
    spot_size: float = 1.0,
    save: str = None,
):
    """Plot QC metrics as spatial heatmaps.

    Parameters
    ----------
    adata : AnnData
        Must have .obsm['spatial'] and QC columns in .obs.
    metrics : list of str
        Columns in .obs to plot. Default: total_counts, n_genes, pct_mt.
    figsize : tuple
        Figure size per row.
    spot_size : float
        Size of scatter points.
    save : str or None
        Path to save figure.
    """
    if metrics is None:
        metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]

    available = [m for m in metrics if m in adata.obs.columns]
    n = len(available)
    if n == 0:
        print("No QC metrics found in adata.obs")
        return

    fig, axes = plt.subplots(1, n, figsize=(figsize[0], figsize[1]))
    if n == 1:
        axes = [axes]

    coords = adata.obsm["spatial"]

    for ax, metric in zip(axes, available):
        values = adata.obs[metric].values
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=values, s=spot_size, cmap="viridis",
            edgecolors="none", rasterized=True,
        )
        ax.set_title(metric, fontsize=12)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
        plt.colorbar(sc_plot, ax=ax, shrink=0.6)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()


def spatial_outlier_detection(
    adata: sc.AnnData,
    metric: str = "total_counts",
    n_neighbors: int = 20,
    z_threshold: float = 3.0,
) -> pd.Series:
    """Detect local spatial outliers using neighborhood z-scores.

    For each bin, compute the z-score of a QC metric relative to its
    spatial neighbors. Bins with |z| > threshold are flagged.

    Parameters
    ----------
    adata : AnnData
        Must have spatial neighbors computed (adata.obsp['spatial_connectivities']).
    metric : str
        Column in adata.obs to test.
    n_neighbors : int
        Number of spatial neighbors (used if neighbors not precomputed).
    z_threshold : float
        Absolute z-score threshold for outlier flagging.

    Returns
    -------
    pd.Series of bool, True = outlier.
    """
    import squidpy as sq

    # Ensure spatial neighbors are computed
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic")

    W = adata.obsp["spatial_connectivities"]
    values = adata.obs[metric].values.astype(np.float64)

    # Compute neighborhood mean and std
    # W is sparse, normalize rows to get mean
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # avoid division by zero

    neighbor_mean = np.array(W.dot(values)).flatten() / row_sums

    # Neighborhood variance: E[X^2] - (E[X])^2
    neighbor_sq_mean = np.array(W.dot(values ** 2)).flatten() / row_sums
    neighbor_var = neighbor_sq_mean - neighbor_mean ** 2
    neighbor_var[neighbor_var < 0] = 0
    neighbor_std = np.sqrt(neighbor_var)
    neighbor_std[neighbor_std == 0] = 1  # avoid division by zero

    z_scores = (values - neighbor_mean) / neighbor_std

    outliers = pd.Series(np.abs(z_scores) > z_threshold, index=adata.obs_names)
    adata.obs[f"{metric}_spatial_zscore"] = z_scores
    adata.obs[f"{metric}_spatial_outlier"] = outliers.values

    return outliers


# ---------------------------------------------------------------------------
# Comparison Helpers
# ---------------------------------------------------------------------------

def compare_binned_vs_segmented(
    adata_binned: sc.AnnData,
    adata_segmented: sc.AnnData,
    cell_type_key_binned: str = "cell_type",
    cell_type_key_seg: str = "cell_type",
    figsize: tuple = (16, 5),
    save: str = None,
):
    """Generate side-by-side comparison plots of binned vs segmented data.

    Creates three panels:
    1. Gene detection comparison (violin)
    2. Spatial cell type maps (side-by-side)
    3. Cell type proportion correlation

    Parameters
    ----------
    adata_binned, adata_segmented : AnnData
        Annotated datasets with cell types and spatial coordinates.
    cell_type_key_binned, cell_type_key_seg : str
        Column names for cell type annotations.
    figsize : tuple
        Figure size.
    save : str or None
        Path to save figure.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # Panel 1: Gene detection comparison
    data_to_plot = []
    if "n_genes_by_counts" in adata_binned.obs:
        data_to_plot.append(
            pd.DataFrame({
                "n_genes": adata_binned.obs["n_genes_by_counts"],
                "Data Type": "Binned (8µm)"
            })
        )
    if "n_genes_by_counts" in adata_segmented.obs:
        data_to_plot.append(
            pd.DataFrame({
                "n_genes": adata_segmented.obs["n_genes_by_counts"],
                "Data Type": "Segmented"
            })
        )

    if data_to_plot:
        import seaborn as sns
        combined = pd.concat(data_to_plot, ignore_index=True)
        sns.violinplot(data=combined, x="Data Type", y="n_genes", ax=axes[0],
                       cut=0, inner="box")
        axes[0].set_title("Genes per Observation")
        axes[0].set_ylabel("Number of Genes")

    # Panel 2: Spatial plot colored by total counts
    for ax, adata, label in [
        (axes[1], adata_binned, "Binned (8µm)"),
    ]:
        if "spatial" in adata.obsm:
            coords = adata.obsm["spatial"]
            ax.scatter(
                coords[:, 0], coords[:, 1],
                c=adata.obs["total_counts"], s=0.5,
                cmap="viridis", edgecolors="none", rasterized=True,
            )
            ax.set_title(f"{label}: Total Counts")
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.axis("off")

    if "spatial" in adata_segmented.obsm:
        coords = adata_segmented.obsm["spatial"]
        axes[2].scatter(
            coords[:, 0], coords[:, 1],
            c=adata_segmented.obs["total_counts"], s=0.5,
            cmap="viridis", edgecolors="none", rasterized=True,
        )
        axes[2].set_title("Segmented: Total Counts")
        axes[2].set_aspect("equal")
        axes[2].invert_yaxis()
        axes[2].axis("off")

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()


def cell_type_proportion_comparison(
    adata_binned: sc.AnnData,
    adata_segmented: sc.AnnData,
    ct_key_binned: str = "cell_type",
    ct_key_seg: str = "cell_type",
    figsize: tuple = (8, 6),
    save: str = None,
):
    """Compare cell type proportions between binned and segmented analysis.

    Parameters
    ----------
    adata_binned, adata_segmented : AnnData
        Annotated datasets.
    ct_key_binned, ct_key_seg : str
        Column names for cell type labels.
    figsize : tuple
        Figure size.
    save : str or None
        Path to save figure.
    """
    from scipy.stats import pearsonr

    # Get proportions
    props_bin = adata_binned.obs[ct_key_binned].value_counts(normalize=True)
    props_seg = adata_segmented.obs[ct_key_seg].value_counts(normalize=True)

    # Align cell types
    all_types = sorted(set(props_bin.index) | set(props_seg.index))
    props_bin = props_bin.reindex(all_types, fill_value=0)
    props_seg = props_seg.reindex(all_types, fill_value=0)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(props_bin, props_seg, s=80, edgecolors="black", linewidths=0.5)

    for ct in all_types:
        ax.annotate(ct, (props_bin[ct], props_seg[ct]),
                    fontsize=8, ha="left", va="bottom")

    # Correlation
    r, p = pearsonr(props_bin, props_seg)
    ax.set_xlabel("Binned (8µm) Proportion")
    ax.set_ylabel("Segmented Proportion")
    ax.set_title(f"Cell Type Proportions: Binned vs Segmented\nr={r:.3f}, p={p:.2e}")

    # Diagonal
    lims = [0, max(props_bin.max(), props_seg.max()) * 1.1]
    ax.plot(lims, lims, "k--", alpha=0.3)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    plt.show()
