"""Synthetic spatial-transcriptomics game utilities.

Builds a fake Visium-HD-style tissue with retro arcade sprites hidden inside
it. The companion notebook (notebook/visium_hd_game_tutorial.ipynb) reveals a
sprite every time the learner performs an analysis step correctly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# ---------------------------------------------------------------------------
# Sprite pixel art (13x13 each). '#' = filled pixel, '.' = empty.
# ---------------------------------------------------------------------------

SPRITES: dict[str, list[str]] = {
    "star": [
        "......#......",
        "......#......",
        ".....###.....",
        "....#####....",
        "#############",
        ".###########.",
        "..#########..",
        "...#######...",
        "...#######...",
        "..##.###.##..",
        ".##...#...##.",
        "##.........##",
        "#...........#",
    ],
    "heart": [
        ".###...###...",
        "#####.#####..",
        "###########..",
        "###########..",
        "###########..",
        "###########..",
        ".#########...",
        "..#######....",
        "...#####.....",
        "....###......",
        ".....#.......",
        ".............",
        ".............",
    ],
    "pacman": [
        "...#####.....",
        "..#######....",
        ".#########...",
        "##########...",
        "#######......",
        "######.......",
        "#####........",
        "######.......",
        "#######......",
        "##########...",
        ".#########...",
        "..#######....",
        "...#####.....",
    ],
    "ghost": [
        "...#####.....",
        "..#######....",
        ".#########...",
        ".##.###.##...",
        ".##.###.##...",
        ".#########...",
        ".#########...",
        ".#########...",
        ".#########...",
        ".#########...",
        ".#########...",
        ".##.##.##....",
        "..#..#..#....",
    ],
    "invader": [
        "..#.....#....",
        "...#...#.....",
        "..#######....",
        ".##.###.##...",
        "#############",
        "#.#######.#..",
        "#.#.....#.#..",
        "...##.##.....",
        ".............",
        ".............",
        ".............",
        ".............",
        ".............",
    ],
    "mushroom": [
        "...#####.....",
        "..#######....",
        ".##.###.##...",
        "#############",
        "##.#####.##..",
        "#############",
        ".#########...",
        "..#######....",
        "...#####.....",
        "....###......",
        "....###......",
        "....###......",
        ".............",
    ],
    "alien": [
        "....###......",
        "...#####.....",
        "..#######....",
        "..#.###.#....",
        "..#######....",
        "..#.#.#.#....",
        "..#######....",
        "...##.##.....",
        "..##...##....",
        ".##.....##...",
        ".............",
        ".............",
        ".............",
    ],
}


# Sprite metadata: display colour, placement offset on the 100x100 grid, and
# cell-type label used for annotation validation.
SPRITE_INFO: dict[str, dict] = {
    "star":     {"offset": (10, 10),  "color": "#FFD60A", "celltype": "StarCell"},
    "heart":    {"offset": (10, 40),  "color": "#FF3B30", "celltype": "HeartCell"},
    "pacman":   {"offset": (40, 10),  "color": "#FFCC00", "celltype": "PacCell"},
    "ghost":    {"offset": (40, 40),  "color": "#FF2D7A", "celltype": "GhostCell"},
    "invader":  {"offset": (70, 10),  "color": "#30D158", "celltype": "InvaderCell"},
    "mushroom": {"offset": (70, 40),  "color": "#BF5AF2", "celltype": "ShroomCell"},
    "alien":    {"offset": (40, 70),  "color": "#64D2FF", "celltype": "AlienCell"},
}

SPRITE_ORDER: list[str] = list(SPRITE_INFO.keys())


# ---------------------------------------------------------------------------
# Synthetic tissue builder
# ---------------------------------------------------------------------------

def _sprite_mask(name: str) -> np.ndarray:
    rows = SPRITES[name]
    return np.array([[c == "#" for c in row] for row in rows], dtype=bool)


def _build_ground_truth(grid_size: int = 100) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (region_labels[N], xs[N], ys[N]) for all grid points."""
    ys, xs = np.mgrid[0:grid_size, 0:grid_size]
    xs = xs.ravel()
    ys = ys.ravel()
    region = np.full(xs.size, "background", dtype=object)

    for name in SPRITE_ORDER:
        mask = _sprite_mask(name)
        x0, y0 = SPRITE_INFO[name]["offset"]
        h, w = mask.shape
        for dy in range(h):
            for dx in range(w):
                if mask[dy, dx]:
                    gx, gy = x0 + dx, y0 + dy
                    idx = gy * grid_size + gx
                    region[idx] = name
    return region, xs, ys


def make_synthetic_tissue(seed: int = 42, grid_size: int = 100) -> ad.AnnData:
    """Build an AnnData that imitates a Visium HD bin x gene matrix.

    Sprites are stamped on a 100x100 spatial grid. Their pixels have strong
    counts and express sprite-specific marker genes. Background bins have low
    counts, high mitochondrial fraction, and random noise expression. Sprites
    only become visible after QC + normalization.
    """
    rng = np.random.default_rng(seed)

    region, xs, ys = _build_ground_truth(grid_size)
    n_obs = xs.size

    # Gene catalog -----------------------------------------------------------
    marker_genes: dict[str, list[str]] = {
        name: [f"{name.upper()}_MK{i+1}" for i in range(3)] for name in SPRITE_ORDER
    }
    all_markers = [g for gs in marker_genes.values() for g in gs]
    housekeeping = [f"HK_{i:03d}" for i in range(50)]
    mito = [f"MT-GENE{i+1}" for i in range(10)]
    noise = [f"NOISE_{i:03d}" for i in range(20)]
    var_names = all_markers + housekeeping + mito + noise
    n_vars = len(var_names)

    gene_type = (
        ["marker"] * len(all_markers)
        + ["housekeeping"] * len(housekeeping)
        + ["mito"] * len(mito)
        + ["noise"] * len(noise)
    )
    gene_idx = {g: i for i, g in enumerate(var_names)}

    # Counts matrix ----------------------------------------------------------
    X = np.zeros((n_obs, n_vars), dtype=np.float32)

    # Background bins: low counts, high mito fraction, random noise genes.
    bg_mask = region == "background"
    n_bg = bg_mask.sum()
    # housekeeping ~2 counts each
    hk_slice = slice(len(all_markers), len(all_markers) + len(housekeeping))
    X[bg_mask, hk_slice] = rng.poisson(0.3, size=(n_bg, len(housekeeping)))
    # Mito genes: dominant in background → high pct_mt
    mt_slice = slice(len(all_markers) + len(housekeeping),
                     len(all_markers) + len(housekeeping) + len(mito))
    X[bg_mask, mt_slice] = rng.poisson(1.2, size=(n_bg, len(mito)))
    # Noise genes: sparse
    noise_slice = slice(len(all_markers) + len(housekeeping) + len(mito), n_vars)
    X[bg_mask, noise_slice] = rng.poisson(0.4, size=(n_bg, len(noise)))

    # Sprite bins: high total counts, low mito, strong marker expression
    for name in SPRITE_ORDER:
        m = region == name
        n_s = m.sum()
        if n_s == 0:
            continue
        # Housekeeping strong
        X[m, hk_slice] = rng.poisson(3.0, size=(n_s, len(housekeeping)))
        # Mito low
        X[m, mt_slice] = rng.poisson(0.2, size=(n_s, len(mito)))
        # Own markers: very high
        for mg in marker_genes[name]:
            X[m, gene_idx[mg]] = rng.poisson(25.0, size=n_s)
        # Other markers: leak a tiny amount
        for other in SPRITE_ORDER:
            if other == name:
                continue
            for mg in marker_genes[other]:
                X[m, gene_idx[mg]] = rng.poisson(0.2, size=n_s)

    # Assemble AnnData -------------------------------------------------------
    obs = pd.DataFrame(
        {
            "region": pd.Categorical(region, categories=["background"] + SPRITE_ORDER),
            "ground_truth_celltype": pd.Categorical(
                [SPRITE_INFO[r]["celltype"] if r in SPRITE_INFO else "Background"
                 for r in region]
            ),
        },
        index=[f"bin_{i:05d}" for i in range(n_obs)],
    )
    var = pd.DataFrame(
        {"gene_type": gene_type},
        index=pd.Index(var_names, name="gene"),
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = np.column_stack([xs, ys]).astype(np.float32)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.uns["sprite_markers"] = marker_genes
    adata.uns["grid_size"] = grid_size
    return adata


# ---------------------------------------------------------------------------
# Game state + reveal
# ---------------------------------------------------------------------------

class GameState:
    """Tracks unlocked sprites and total score across the notebook."""

    POINTS_PER_SPRITE = 10
    TOTAL_SPRITES = len(SPRITE_ORDER)

    def __init__(self) -> None:
        self.unlocked: list[str] = []
        self.score: int = 0

    def unlock(self, name: str) -> bool:
        if name in self.unlocked:
            return False
        self.unlocked.append(name)
        self.score += self.POINTS_PER_SPRITE
        return True

    def banner(self) -> str:
        bar = "".join("■" if s in self.unlocked else "□" for s in SPRITE_ORDER)
        return f"[{bar}]  {len(self.unlocked)}/{self.TOTAL_SPRITES} sprites   score: {self.score}"

    def final_message(self) -> str:
        if len(self.unlocked) == self.TOTAL_SPRITES:
            return "🏆  LEVEL CLEAR — you mastered the spatial pipeline!"
        return f"Game over. Collected {len(self.unlocked)}/{self.TOTAL_SPRITES}. Retry any step to grab the rest!"


def _grid_image(adata: ad.AnnData, values: np.ndarray) -> np.ndarray:
    n = int(adata.uns["grid_size"])
    img = np.full((n, n), np.nan, dtype=np.float32)
    xs = adata.obsm["spatial"][:, 0].astype(int)
    ys = adata.obsm["spatial"][:, 1].astype(int)
    img[ys, xs] = values
    return img


def reveal_sprite(adata: ad.AnnData, name: str, state: GameState,
                  title_suffix: str = "") -> None:
    """Plot the sprite's marker-expression image and celebrate."""
    markers = adata.uns["sprite_markers"][name]
    expr = np.asarray(adata[:, markers].X).sum(axis=1)
    if hasattr(expr, "A1"):
        expr = expr.A1
    img = _grid_image(adata, expr)

    cmap = ListedColormap(["#111111", SPRITE_INFO[name]["color"]])
    binary = (img > np.nanpercentile(img, 90)).astype(float)

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(binary, cmap=cmap, origin="lower", interpolation="nearest")
    ax.set_title(f"🎉  Sprite unlocked: {name.upper()}  (+{GameState.POINTS_PER_SPRITE})"
                 + (f"\n{title_suffix}" if title_suffix else ""))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    plt.show()
    print(state.banner())


def show_raw_tissue(adata: ad.AnnData) -> None:
    """Plot total counts per bin — raw view, noisy by design."""
    totals = np.asarray(adata.X).sum(axis=1)
    if hasattr(totals, "A1"):
        totals = totals.A1
    img = _grid_image(adata, totals)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.imshow(img, cmap="magma", origin="lower")
    ax.set_title("Raw tissue (total counts per bin)\n"
                 "sprites are buried in background noise")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Step validators
# ---------------------------------------------------------------------------

def _sprite_indices(adata: ad.AnnData, name: str) -> np.ndarray:
    return np.flatnonzero(adata.obs["region"].to_numpy() == name)


def validate_load(adata: ad.AnnData, state: GameState) -> bool:
    if adata.n_obs == 0:
        print("❌  No bins loaded.")
        return False
    print(f"✅  Loaded {adata.n_obs} bins × {adata.n_vars} genes.")
    if state.unlock("star"):
        reveal_sprite(adata, "star", state, "Welcome aboard!")
    return True


def validate_qc(adata_filtered: ad.AnnData, original: ad.AnnData,
                state: GameState) -> bool:
    """Pass if user removed ≥80% of background AND kept ≥80% of sprite bins."""
    kept_ids = set(adata_filtered.obs_names)
    orig_region = original.obs["region"].to_numpy()

    bg_total = (orig_region == "background").sum()
    bg_kept = sum(1 for i, r in enumerate(orig_region)
                  if r == "background" and original.obs_names[i] in kept_ids)
    sprite_total = (orig_region != "background").sum()
    sprite_kept = sum(1 for i, r in enumerate(orig_region)
                      if r != "background" and original.obs_names[i] in kept_ids)

    bg_removed_frac = 1 - bg_kept / max(bg_total, 1)
    sprite_kept_frac = sprite_kept / max(sprite_total, 1)
    print(f"background removed: {bg_removed_frac:.0%}   sprites kept: {sprite_kept_frac:.0%}")

    if bg_removed_frac >= 0.8 and sprite_kept_frac >= 0.8:
        if state.unlock("heart"):
            reveal_sprite(adata_filtered, "heart", state, "Clean QC = happy tissue ❤")
        return True
    print("Hint: aim for a threshold that kills noisy low-count, high-mito bins"
          " without touching the bright spots.")
    return False


def validate_normalization(adata: ad.AnnData, state: GameState) -> bool:
    """Pass if data was log-normalized (max < 20 and floats)."""
    x = np.asarray(adata.X)
    if x.dtype.kind not in {"f"}:
        print("❌  Matrix is still integer — did you normalize + log1p?")
        return False
    if x.max() > 20:
        print("❌  Values look un-normalized (max > 20). Try sc.pp.normalize_total + sc.pp.log1p.")
        return False
    print(f"✅  Normalized matrix (max value = {x.max():.2f}).")
    if state.unlock("mushroom"):
        reveal_sprite(adata, "mushroom", state, "1-UP! Normalized counts.")
    return True


def validate_clusters(adata: ad.AnnData, cluster_key: str,
                      state: GameState) -> bool:
    """Pass if cluster labels align with ground-truth regions (ARI-ish check)."""
    if cluster_key not in adata.obs:
        print(f"❌  obs['{cluster_key}'] not found.")
        return False
    clusters = adata.obs[cluster_key].astype(str).to_numpy()
    regions = adata.obs["region"].astype(str).to_numpy()

    matched_sprites = 0
    for name in SPRITE_ORDER:
        sprite_bins = regions == name
        if sprite_bins.sum() == 0:
            continue
        dominant = pd.Series(clusters[sprite_bins]).mode().iloc[0]
        in_cluster = clusters == dominant
        # precision = fraction of the dominant cluster that IS this sprite
        precision = (in_cluster & sprite_bins).sum() / max(in_cluster.sum(), 1)
        recall = (in_cluster & sprite_bins).sum() / sprite_bins.sum()
        if precision > 0.5 and recall > 0.6:
            matched_sprites += 1

    print(f"clusters matching a sprite region: {matched_sprites}/{len(SPRITE_ORDER)}")
    if matched_sprites >= 5:
        if state.unlock("pacman"):
            reveal_sprite(adata, "pacman", state, "Waka waka — clusters found!")
        return True
    print("Hint: try increasing Leiden resolution so each sprite gets its own cluster.")
    return False


def validate_annotation(adata: ad.AnnData, marker_dict: dict[str, list[str]],
                        state: GameState) -> bool:
    """User provides {celltype: [marker_genes]}. Must correctly tag GhostCell."""
    target = "GhostCell"
    if target not in marker_dict:
        print(f"❌  '{target}' not in your marker dict. Hint: find genes that light up the ghost region.")
        return False
    chosen = marker_dict[target]
    expected = set(adata.uns["sprite_markers"]["ghost"])
    overlap = len(expected.intersection(chosen))
    if overlap < 2:
        print(f"❌  Only {overlap}/3 correct ghost markers. "
              "Look for genes prefixed with 'GHOST_'.")
        return False
    print(f"✅  {overlap}/3 ghost markers correct.")
    if state.unlock("ghost"):
        reveal_sprite(adata, "ghost", state, "Boo! Annotation complete.")
    return True


def validate_neighborhoods(adata: ad.AnnData, state: GameState) -> bool:
    """Pass if a spatial neighbor graph lives in obsp."""
    key = None
    for k in adata.obsp.keys() if hasattr(adata.obsp, "keys") else []:
        if "spatial" in k or "connect" in k:
            key = k
            break
    if key is None and "spatial_connectivities" in adata.obsp:
        key = "spatial_connectivities"
    if key is None:
        print("❌  No spatial neighbor graph in adata.obsp. "
              "Hint: call build_spatial_graph(adata).")
        return False
    graph = adata.obsp[key]
    if graph.nnz == 0:
        print("❌  Neighbor graph is empty.")
        return False
    print(f"✅  Spatial neighbor graph with {graph.nnz} edges.")
    if state.unlock("invader"):
        reveal_sprite(adata, "invader", state, "Neighbors assembled!")
    return True


def validate_svgs(svg_list: list[str], adata: ad.AnnData,
                  state: GameState) -> bool:
    """Pass if user's top SVGs include ≥2 alien markers."""
    expected = set(adata.uns["sprite_markers"]["alien"])
    overlap = len(expected.intersection(set(svg_list)))
    if overlap < 2:
        print(f"❌  Only {overlap}/3 alien markers in your SVG list. "
              "Hint: the alien's genes vary sharply across space.")
        return False
    print(f"✅  {overlap}/3 alien markers present in top SVGs.")
    if state.unlock("alien"):
        reveal_sprite(adata, "alien", state, "Contact made! 👽")
    return True


# ---------------------------------------------------------------------------
# Helpers used by the notebook (non-validating)
# ---------------------------------------------------------------------------

def build_spatial_graph(adata: ad.AnnData, radius: float = 1.5) -> None:
    """Tiny radius-based neighbor graph → adata.obsp['spatial_connectivities']."""
    from scipy.sparse import csr_matrix
    from scipy.spatial import cKDTree

    coords = adata.obsm["spatial"]
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=radius, output_type="ndarray")
    if pairs.size == 0:
        adata.obsp["spatial_connectivities"] = csr_matrix((adata.n_obs, adata.n_obs))
        return
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    data = np.ones_like(rows, dtype=np.float32)
    adata.obsp["spatial_connectivities"] = csr_matrix(
        (data, (rows, cols)), shape=(adata.n_obs, adata.n_obs)
    )


def morans_i(adata: ad.AnnData, top_n: int = 10) -> list[tuple[str, float]]:
    """Minimal Moran's I over the spatial neighbor graph. Returns top_n genes."""
    if "spatial_connectivities" not in adata.obsp:
        build_spatial_graph(adata)
    W = adata.obsp["spatial_connectivities"]
    W_sum = W.sum()
    if W_sum == 0:
        return []
    X = np.asarray(adata.X)
    mean = X.mean(axis=0)
    dev = X - mean
    num_coef = (adata.n_obs / W_sum)
    # For each gene: I = N/W * sum_{ij} W_ij * dev_i * dev_j / sum_i dev_i^2
    Wdev = W @ dev                       # N x G
    numer = (dev * Wdev).sum(axis=0)
    denom = (dev * dev).sum(axis=0) + 1e-12
    I = num_coef * numer / denom
    order = np.argsort(-I)[:top_n]
    return [(adata.var_names[i], float(I[i])) for i in order]
