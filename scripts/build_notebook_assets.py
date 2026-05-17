"""Generate the two images embedded in visium_hd_quickstart.ipynb.

Produces:
  1. A terminal-style PNG of the Space Ranger ``CRC_P1/outs/`` tree.
  2. A PNG raster of the official AnnData schema figure pulled from the
     scverse/anndata docs.

Both are base64-embedded as data URIs into the notebook markdown so the
notebook renders on Colab with no external files. Re-running this script is
idempotent: it replaces the generated cells in place by their stable ids.

Usage:
    python scripts/build_notebook_assets.py
"""

import base64
import io
import json
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from PIL import Image

HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE.parent / "notebook" / "visium_hd_quickstart.ipynb"
ASSETS = HERE.parent / "notebook" / "assets"

ANNDATA_SCHEMA_SVG = (
    "https://raw.githubusercontent.com/scverse/anndata/main/"
    "docs/_static/img/anndata_schema.svg"
)

# Space Ranger outs/ tree. (left, comment) per line; left is rendered light,
# comment muted green, directory names cyan.
TREE_LINES = [
    ("$ tree CRC_P1/outs/", None),
    ("CRC_P1/outs/", None),
    ("├── binned_outputs/", None),
    ("│   ├── square_002um/", "# 2 um bins (native res, sparsest)"),
    ("│   ├── square_008um/", "# 8 um bins (what we load)"),
    ("│   │   ├── filtered_feature_bc_matrix/", "# counts, MEX format"),
    ("│   │   ├── filtered_feature_bc_matrix.h5", "# same counts, HDF5"),
    ("│   │   └── spatial/", None),
    ("│   │       ├── tissue_positions.parquet", "# bin coords on the H&E"),
    ("│   │       ├── scalefactors_json.json", "# pixel-to-um factors"),
    ("│   │       ├── tissue_hires_image.png", "# downscaled H&E"),
    ("│   │       └── tissue_lowres_image.png", "# smaller H&E"),
    ("│   └── square_016um/", "# 16 um bins (quick QC)"),
    ("├── segmented_outputs/", "# v4.0+: one row per cell"),
    ("│   └── filtered_feature_cell_matrix.h5", None),
    ("├── cloupe_008um.cloupe", "# Loupe Browser viewer file"),
    ("├── web_summary.html", "# QC report (read this first)"),
    ("└── metrics_summary.csv", "# QC metrics, machine-readable"),
]

BG = "#1e1e1e"
FG = "#d4d4d4"
GREEN = "#6a9955"
CYAN = "#4ec9b0"
PROMPT = "#9cdcfe"


def _mono_font():
    for name in ("DejaVu Sans Mono", "Liberation Mono", "Courier New"):
        try:
            return font_manager.FontProperties(family=name)
        except Exception:
            continue
    return font_manager.FontProperties(family="monospace")


def make_terminal_png() -> bytes:
    """Render the outs/ tree as a dark terminal screenshot."""
    rows = [left + (("   " + cmt) if cmt else "") for left, cmt in TREE_LINES]
    ncols = max(len(r) for r in rows)
    nrows = len(rows)

    fig_w = 0.085 * ncols + 0.6
    fig_h = 0.30 * nrows + 0.7
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, ncols)
    ax.set_ylim(0, nrows + 1.4)
    ax.invert_yaxis()
    ax.axis("off")
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    # Window chrome dots.
    for i, c in enumerate(("#ff5f56", "#ffbd2e", "#27c93f")):
        ax.scatter(0.6 + i * 1.0, 0.55, s=130, color=c, zorder=5,
                   edgecolors="none", clip_on=False)

    fp = _mono_font()
    char_w = 1.0  # data units per character (xlim == ncols)

    for i, (left, cmt) in enumerate(TREE_LINES):
        y = i + 1.6
        if i == 0:  # prompt line
            ax.text(0, y, "$", color=GREEN, fontproperties=fp, fontsize=11,
                    ha="left", va="center")
            ax.text(2 * char_w, y, left[2:], color=PROMPT, fontproperties=fp,
                    fontsize=11, ha="left", va="center")
            continue
        is_dir = left.rstrip().endswith("/")
        ax.text(0, y, left, color=(CYAN if is_dir else FG),
                fontproperties=fp, fontsize=11, ha="left", va="center")
        if cmt:
            ax.text((len(left) + 3) * char_w, y, cmt, color=GREEN,
                    fontproperties=fp, fontsize=11, ha="left", va="center")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=BG, bbox_inches="tight",
                pad_inches=0.18)
    plt.close(fig)
    return _quantize(buf.getvalue(), colors=32)


def make_anndata_png() -> bytes:
    """Download the official scverse AnnData schema SVG, rasterize to PNG."""
    import cairosvg

    with urllib.request.urlopen(ANNDATA_SCHEMA_SVG, timeout=60) as resp:
        svg = resp.read()
    png = cairosvg.svg2png(bytestring=svg, output_width=1200,
                           background_color="white")
    return _quantize(png, colors=64)


def _quantize(png_bytes: bytes, colors: int) -> bytes:
    """Shrink a flat-color PNG via palette quantization."""
    im = Image.open(io.BytesIO(png_bytes)).convert("RGB")
    im = im.quantize(colors=colors, method=Image.MEDIANCUT)
    out = io.BytesIO()
    im.save(out, format="PNG", optimize=True)
    return out.getvalue()


def _data_uri(png_bytes: bytes) -> str:
    return "data:image/png;base64," + base64.b64encode(png_bytes).decode("ascii")


def build_cells(term_uri: str, anndata_uri: str):
    """Return the four markdown cells keyed by their stable ids."""
    sr_intro = (
        "## From tissue to matrix: Space Ranger\n"
        "\n"
        "Before any of this notebook runs, the raw data has to come out of "
        "**Space Ranger**, the 10x pipeline that turns sequencing reads and "
        "the H&E image into a spatial gene-by-bin count matrix.\n"
        "\n"
        "**What goes in:**\n"
        "\n"
        "- FASTQ files: the raw Illumina reads, carrying the spatial barcode, "
        "UMI, and gene sequence for every captured molecule.\n"
        "- Reference transcriptome: the prebuilt 10x reference for the species "
        "(human GRCh38 here).\n"
        "- Probe set CSV: the Visium HD probe list, matched to the assay "
        "chemistry.\n"
        "- Slide serial number and capture area ID: printed on the slide, used "
        "to look up the fiducial layout.\n"
        "- CytAssist and microscope H&E images: aligned so barcodes land on "
        "the right tissue coordinates.\n"
        "\n"
        "**What comes out (the parts we use):**\n"
        "\n"
        "- `filtered_feature_bc_matrix.h5`: the gene-by-bin counts in one HDF5 "
        "file. This is what we read into AnnData.\n"
        "- `spatial/tissue_positions.parquet`: each bin's pixel coordinates on "
        "the H&E. This is what makes the data spatial.\n"
        "- `spatial/scalefactors_json.json`: pixel-to-um conversion, needed to "
        "overlay bins on the image.\n"
        "- `spatial/tissue_hires_image.png` and `tissue_lowres_image.png`: the "
        "H&E to plot under the bins.\n"
        "\n"
        "Space Ranger also writes `segmented_outputs/` (one row per cell, from "
        "nuclei segmentation), `cloupe` files for the Loupe Browser, and "
        "`web_summary.html` for QC. We don't use those here, but they sit in "
        "the same output folder, shown below."
    )

    sr_tree = (
        f"![Space Ranger CRC_P1/outs/ directory tree]({term_uri})\n"
        "\n"
        "*The Space Ranger output folder. We pull four files from "
        "`square_008um/` and its `spatial/` subfolder.*"
    )

    sr_res = (
        "### Which resolution we use\n"
        "\n"
        "Space Ranger writes the same tissue at three bin sizes. We load "
        "`binned_outputs/square_008um/filtered_feature_bc_matrix.h5`.\n"
        "\n"
        "8 um sits at roughly 1 to 3 cells per bin, which keeps enough counts "
        "per bin to cluster and annotate without needing nuclei segmentation "
        "of the H&E. The 2 um bins are at native resolution but too sparse to "
        "use directly, and the 16 um bins blur too many cell types together "
        "for annotation. When you need true single-cell resolution, the route "
        "is `segmented_outputs/`, not a finer bin."
    )

    anndata_md = (
        "## What's inside the AnnData object\n"
        "\n"
        f"![AnnData data structure]({anndata_uri})\n"
        "\n"
        "*AnnData data structure. Source: anndata documentation (scverse), "
        "https://anndata.readthedocs.io*\n"
        "\n"
        "Everything we just loaded lives in one AnnData object, `adata_8um`. "
        "The pieces:\n"
        "\n"
        "- **`X`**: the counts. A sparse matrix, bins as rows, genes as "
        "columns.\n"
        "- **`obs`**: the per-bin table, one row per bin, barcodes as the "
        "index. The ground-truth `UnsupervisedL1` labels are a column here.\n"
        "- **`var`**: the per-gene table, one row per gene.\n"
        "- **`obsm`**: per-bin matrices, keyed by name. This is the difference "
        "from `obs`: `obs` holds one value per bin per column, `obsm` holds a "
        "whole array per bin. The spatial coordinates are `obsm['spatial']`, "
        "and the UMAP computed later in this notebook lands here as "
        "`obsm['X_umap']`.\n"
        "- **`uns`**: unstructured extras. The H&E image is "
        "`uns['spatial_image']` and the pixel-to-um factors are "
        "`uns['scalefactors']`.\n"
        "- **`layers`** and **`raw`**: alternate copies of the matrix, for "
        "example counts kept before normalization.\n"
        "\n"
        "`obs` and `var` describe the two axes, `obsm` and `varm` attach "
        "matrices to them, and `uns` is the catch-all. The figure shows the "
        "general layout; not every slot is filled in this dataset."
    )

    return {
        "sr-intro-md": sr_intro,
        "sr-tree-img": sr_tree,
        "sr-res-md": sr_res,
        "anndata-struct-md": anndata_md,
    }


def _md_cell(cell_id: str, source: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": source.splitlines(keepends=True),
    }


def patch_notebook(cells: dict):
    nb = json.loads(NOTEBOOK.read_text())
    our_ids = set(cells)

    # Drop any previously inserted copies so re-runs stay idempotent.
    kept = [c for c in nb["cells"] if c.get("id") not in our_ids]

    def find(pred):
        for i, c in enumerate(kept):
            if pred(c):
                return i
        raise RuntimeError("anchor cell not found")

    intro_i = find(
        lambda c: c["cell_type"] == "markdown"
        and "".join(c["source"]).lstrip().startswith(
            "# Visium HD Spatial Transcriptomics Quickstart"
        )
    )
    block1 = [
        _md_cell("sr-intro-md", cells["sr-intro-md"]),
        _md_cell("sr-tree-img", cells["sr-tree-img"]),
        _md_cell("sr-res-md", cells["sr-res-md"]),
    ]
    kept[intro_i + 1 : intro_i + 1] = block1

    he_i = find(
        lambda c: c["cell_type"] == "markdown"
        and "".join(c["source"]).lstrip().startswith("## H&E Image Verification")
    )
    kept[he_i:he_i] = [_md_cell("anndata-struct-md", cells["anndata-struct-md"])]

    nb["cells"] = kept
    NOTEBOOK.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")


def main():
    ASSETS.mkdir(parents=True, exist_ok=True)
    term_png = make_terminal_png()
    anndata_png = make_anndata_png()
    (ASSETS / "spaceranger_outs_tree.png").write_bytes(term_png)
    (ASSETS / "anndata_schema.png").write_bytes(anndata_png)
    print(f"terminal PNG: {len(term_png) / 1024:.1f} KB")
    print(f"anndata PNG:  {len(anndata_png) / 1024:.1f} KB")

    cells = build_cells(_data_uri(term_png), _data_uri(anndata_png))
    patch_notebook(cells)
    print(f"Patched {NOTEBOOK.relative_to(HERE.parent.parent)}")


if __name__ == "__main__":
    main()
