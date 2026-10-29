#!/usr/bin/env python3
import argparse
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.colors import BoundaryNorm, ListedColormap, to_hex
from matplotlib.patches import Rectangle
from scipy import sparse
from scipy.io import mmread


CLUSTER_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255",
    "#AA4499", "#661100", "#6699CC", "#AA4466",
    "#4477AA", "#228833", "#EE6677", "#BBBBBB",
]

PALETTE = {
    "paper": "#f7f5ef",
    "panel": "#ffffff",
    "text": "#26333c",
    "muted": "#65737c",
    "grid": "#d9dde1",
}


def cluster_colors(n_clusters):
    """Return one deterministic display color per Leiden cluster.

    Use the project-wide named Nature/NPG mapping. More than 32 domains need
    an explicit reviewed extension rather than a position-dependent fallback.
    """
    n_clusters = int(n_clusters)
    if n_clusters <= 0:
        raise ValueError("n_clusters must be positive")
    if n_clusters > len(CLUSTER_COLORS):
        raise ValueError("More than 32 spatial domains require an explicit palette")
    return CLUSTER_COLORS[:n_clusters]


def parse_args():
    parser = argparse.ArgumentParser(description="Run Scanpy expression clustering and UMAP for Spatial SLAM bins.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--matrix-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--bin-size-um", type=float)
    parser.add_argument("--bin-size-pixels", type=int)
    parser.add_argument("--min-bin-umis", type=int, default=20)
    parser.add_argument("--hvg", type=int, default=3000)
    parser.add_argument("--pca-components", type=int, default=25)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--target-clusters", type=int, default=6)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--formats", default="png,pdf")
    parser.add_argument(
        "--assignments",
        default="total",
        help="Comma-separated analyses to run: total,old,new. The reviewed workflow defaults to total RNA only.",
    )
    parser.add_argument("--write-h5ad", action="store_true")
    parser.add_argument(
        "--celatlas-bin-dir", "--celatlas-bin50-dir", dest="celatlas_bin_dir",
        help="Official Celatlas bin directory containing filtered_feature_bc_matrix and spatial/tissue_positions_list.csv",
    )
    parser.add_argument("--leiden-resolution", type=float, help="Fixed Leiden resolution; when set, skip target-cluster auto-selection")
    parser.add_argument("--leiden-random-state", type=int, default=0)
    parser.add_argument("--leiden-n-iterations", type=int, default=2)
    parser.add_argument("--hvg-flavor", default="seurat", choices=["cell_ranger", "seurat", "seurat_v3"])
    parser.add_argument("--umap-point-size", type=float, default=10.0)
    return parser.parse_args()


def import_scanpy():
    import scanpy as sc

    return sc


def save_figure(fig, outdir, stem, formats, dpi):
    paths = []
    for fmt in formats:
        path = outdir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        paths.append(str(path))
    plt.close(fig)
    return paths


def load_matrix_inputs(matrix_dir, sample):
    matrix_dir = Path(matrix_dir)
    feature_meta = pd.read_csv(matrix_dir / "features.tsv", sep="\t", dtype={"gene_id": str})
    barcode_meta = pd.read_csv(matrix_dir / "barcodes.tsv", sep="\t")
    features = feature_meta["gene_id"].astype(str).tolist()
    barcodes = barcode_meta["barcode"].astype(str).tolist()
    barcode_meta["barcode"] = barcode_meta["barcode"].astype(str)
    barcode_meta = barcode_meta.set_index("barcode", drop=False)
    if barcode_meta[["x", "y"]].isna().any().any():
        raise ValueError("barcode metadata is missing x/y for some matrix barcodes")
    return matrix_dir, features, barcodes, barcode_meta


def load_assignment_matrix(matrix_dir, sample, assignment):
    if assignment == "total":
        return (
            mmread(matrix_dir / "old_umi.mtx").tocsr()
            + mmread(matrix_dir / "new_umi.mtx").tocsr()
        )
    return mmread(matrix_dir / f"{assignment}_umi.mtx").tocsr()


def load_official_bin(celatlas_dir):
    """Load a filtered official Celatlas binned matrix and its spatial grid."""
    import anndata as ad

    root = Path(celatlas_dir)
    matrix_dir = root / "filtered_feature_bc_matrix"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    if not matrix_path.exists():
        matrix_path = matrix_dir / "matrix.mtx"
    features_path = matrix_dir / "features.tsv.gz"
    if not features_path.exists():
        features_path = matrix_dir / "features.tsv"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    if not barcodes_path.exists():
        barcodes_path = matrix_dir / "barcodes.tsv"
    opener = gzip.open if str(matrix_path).endswith(".gz") else open
    with opener(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr().T
    feature_meta = pd.read_csv(features_path, sep="\t", header=None, compression="infer", dtype=str)
    features = feature_meta.iloc[:, 0].astype(str).tolist()
    barcode_series = pd.read_csv(barcodes_path, sep="\t", header=None, compression="infer", dtype=str).iloc[:, 0]
    barcodes = barcode_series.astype(str).tolist()
    positions_path = root / "spatial" / "tissue_positions_list.csv"
    if not positions_path.exists():
        raise FileNotFoundError(f"official bin positions missing: {positions_path}")
    positions = pd.read_csv(positions_path, header=None, names=["barcode", "in_tissue", "row", "col", "pxl_row_in_fullres", "pxl_col_in_fullres"], dtype={"barcode": str})
    positions = positions[positions["barcode"].isin(barcodes)].copy()
    positions = positions.set_index("barcode").reindex(barcodes)
    required_coordinates = ["row", "col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if positions[required_coordinates].isna().any().any():
        raise ValueError("official bin positions do not cover all matrix barcodes")
    scalefactors_path = root / "spatial" / "scalefactors_json.json"
    if not scalefactors_path.exists():
        raise FileNotFoundError(f"official bin scalefactors missing: {scalefactors_path}")
    scalefactors = json.loads(scalefactors_path.read_text(encoding="utf-8"))
    spot_diameter = float(scalefactors["spot_diameter_fullres"])
    obs = positions.copy()
    obs.index = pd.Index(barcodes, name="bin_id")
    obs["bin_id"] = barcodes
    obs["bin_x"] = obs["pxl_col_in_fullres"].astype(float)
    obs["bin_y"] = obs["pxl_row_in_fullres"].astype(float)
    obs["display_diameter"] = spot_diameter
    obs["umi_count"] = np.asarray(matrix.sum(axis=1)).ravel().astype(int)
    obs["assignment"] = "total"
    var = pd.DataFrame(index=pd.Index(features, name="gene_id"))
    var["gene_id"] = features
    return ad.AnnData(X=matrix, obs=obs, var=var)


def aggregate_to_bins(matrix_dir, sample, assignment, features, barcodes, barcode_meta, bin_size, min_bin_umis):
    raw = load_assignment_matrix(matrix_dir, sample, assignment)
    meta = barcode_meta.copy()
    meta["bin_x"] = (meta["x"] // bin_size).astype(int) * bin_size + bin_size / 2.0
    meta["bin_y"] = (meta["y"] // bin_size).astype(int) * bin_size + bin_size / 2.0
    meta["bin_id"] = meta["bin_x"].astype(int).astype(str) + "_" + meta["bin_y"].astype(int).astype(str)
    bin_ids = sorted(meta["bin_id"].unique())
    bin_index = {bin_id: idx for idx, bin_id in enumerate(bin_ids)}
    target_bins = meta["bin_id"].map(bin_index).to_numpy(dtype=int)
    barcode_to_bin = sparse.coo_matrix(
        (
            np.ones(len(barcodes), dtype=np.float32),
            (np.arange(len(barcodes)), target_bins),
        ),
        shape=(len(barcodes), len(bin_ids)),
        dtype=np.float32,
    ).tocsr()
    bin_gene = (raw @ barcode_to_bin).T.tocsr()
    totals = np.asarray(bin_gene.sum(axis=1)).ravel()
    keep = totals >= min_bin_umis
    bin_gene = bin_gene[keep]
    kept_bins = [bin_id for bin_id, flag in zip(bin_ids, keep) if flag]
    obs = pd.DataFrame(index=kept_bins)
    obs["bin_id"] = kept_bins
    coords = (
        meta[["bin_id", "bin_x", "bin_y"]]
        .drop_duplicates("bin_id")
        .set_index("bin_id")
        .reindex(kept_bins)
    )
    obs["bin_x"] = coords["bin_x"].astype(float).to_numpy()
    obs["bin_y"] = coords["bin_y"].astype(float).to_numpy()
    obs["assignment"] = assignment
    obs["umi_count"] = totals[keep].astype(int)
    var = pd.DataFrame(index=features)
    var["gene_id"] = features
    return bin_gene, obs, var


def choose_leiden_resolution(sc, adata, target_clusters, seed):
    candidates = [0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.25, 1.5, 2.0]
    best_key = None
    best_score = None
    best_resolution = None
    for resolution in candidates:
        key = f"leiden_{str(resolution).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=resolution, random_state=seed, key_added=key)
        n_clusters = int(adata.obs[key].nunique())
        score = (abs(n_clusters - target_clusters), n_clusters < 2, resolution)
        if best_score is None or score < best_score:
            best_key = key
            best_score = score
            best_resolution = resolution
        if n_clusters == target_clusters:
            break
    labels = adata.obs[best_key].astype(str)
    counts = labels.value_counts()
    ordered = counts.index.tolist()
    remap = {label: str(idx + 1) for idx, label in enumerate(ordered)}
    adata.obs["cluster"] = labels.map(remap).astype("category")
    return best_resolution, int(labels.nunique())


def run_fixed_leiden(sc, adata, args):
    sc.tl.leiden(
        adata,
        resolution=float(args.leiden_resolution),
        random_state=args.leiden_random_state,
        n_iterations=args.leiden_n_iterations,
        key_added="leiden_fixed",
    )
    labels = adata.obs["leiden_fixed"].astype(str)
    counts = labels.value_counts()
    remap = {label: str(idx + 1) for idx, label in enumerate(counts.index.tolist())}
    adata.obs["cluster"] = labels.map(remap).astype("category")
    return float(args.leiden_resolution), int(labels.nunique())


def run_scanpy(adata, args):
    sc = import_scanpy()
    sc.settings.verbosity = 1
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.hvg, flavor=args.hvg_flavor)
    if "highly_variable" in adata.var and int(adata.var["highly_variable"].sum()) >= 50:
        adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    n_comps = max(2, min(args.pca_components, adata.n_obs - 1, adata.n_vars - 1))
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=n_comps, random_state=args.seed)
    if args.leiden_resolution is not None:
        resolution, raw_clusters = run_fixed_leiden(sc, adata, args)
    else:
        resolution, raw_clusters = choose_leiden_resolution(sc, adata, args.target_clusters, args.seed)
    sc.tl.umap(adata, random_state=args.seed)
    return adata, {
        "algorithm": "scanpy_normalize_log1p_hvg_pca_neighbors_leiden_umap",
        "scanpy_version": sc.__version__,
        "n_hvg": int(adata.n_vars),
        "n_pcs": int(n_comps),
        "neighbors": int(args.neighbors),
        "leiden_resolution": float(resolution),
        "raw_leiden_clusters": int(raw_clusters),
        "coord_weight": 0.0,
        "hvg_flavor": args.hvg_flavor,
        "leiden_random_state": int(args.leiden_random_state if args.leiden_resolution is not None else args.seed),
        "leiden_n_iterations": int(args.leiden_n_iterations if args.leiden_resolution is not None else 1),
    }


def output_frame(adata, assignment):
    columns = ["bin_id", "bin_x", "bin_y", "cluster", "umi_count"]
    if "display_diameter" in adata.obs:
        columns.append("display_diameter")
    frame = adata.obs[columns].copy()
    frame["cluster"] = frame["cluster"].astype(int)
    frame["assignment"] = assignment
    frame["UMAP1"] = adata.obsm["X_umap"][:, 0]
    frame["UMAP2"] = adata.obsm["X_umap"][:, 1]
    return frame.reset_index(drop=True)


def common_spatial_style(ax, title, frame):
    ax.set_title(title, fontsize=13, fontweight="bold", color=PALETTE["text"], pad=10)
    ax.set_xlabel("spatial1", fontsize=9, color=PALETTE["muted"])
    ax.set_ylabel("spatial2", fontsize=9, color=PALETTE["muted"])
    ax.set_aspect("equal", adjustable="box")
    pad = 130
    ax.set_xlim(float(frame["bin_x"].min()) - pad, float(frame["bin_x"].max()) + pad)
    ax.set_ylim(float(frame["bin_y"].min()) - pad, float(frame["bin_y"].max()) + pad)
    ax.invert_yaxis()
    ax.tick_params(labelsize=7, colors=PALETTE["muted"], length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)


def assignment_label(assignment):
    return "Total" if assignment == "total" else assignment.capitalize()


def spatial_stem(assignment, bin_size):
    if assignment == "total":
        return f"06_spatial_total_domains_bin{bin_size}"
    return f"06_spatial_{assignment}_clusters_bin{bin_size}"


def umap_stem(assignment):
    if assignment == "total":
        return "07_umap_total_domains"
    return f"07_umap_{assignment}_expression_clusters"


def plot_spatial(frame, assignment, args):
    max_cluster = int(frame["cluster"].max())
    colors = cluster_colors(max_cluster)
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(0.5, max_cluster + 1.5, 1), len(colors))
    fig, ax = plt.subplots(figsize=(7.2, 7.0), facecolor=PALETTE["paper"])
    ax.set_facecolor(PALETTE["panel"])
    display_bin_size = float(frame["display_diameter"].iloc[0]) if "display_diameter" in frame else args.bin_size
    patches = [
        Rectangle((row.bin_x - display_bin_size / 2.0, row.bin_y - display_bin_size / 2.0), display_bin_size, display_bin_size)
        for row in frame.itertuples(index=False)
    ]
    collection = PatchCollection(patches, linewidths=0, edgecolors="none", rasterized=True)
    collection.set_array(frame["cluster"].to_numpy(dtype=float))
    collection.set_cmap(cmap)
    collection.set_norm(norm)
    ax.add_collection(collection)
    common_spatial_style(ax, f"{assignment_label(assignment)} RNA Scanpy Leiden domains (bin{args.bin_size})", frame)
    cbar = fig.colorbar(collection, ax=ax, ticks=np.arange(1, max_cluster + 1), fraction=0.035, pad=0.02)
    cbar.ax.set_yticklabels([f"C{i}" for i in range(1, max_cluster + 1)])
    cbar.set_label("Leiden cluster", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    stem = spatial_stem(assignment, args.bin_size)
    return save_figure(fig, Path(args.outdir), stem, args.formats, args.dpi)


def plot_umap(frame, assignment, args):
    max_cluster = int(frame["cluster"].max())
    colors = cluster_colors(max_cluster)
    cmap = ListedColormap(colors)
    fig, ax = plt.subplots(figsize=(7.4, 6.8), facecolor=PALETTE["paper"])
    ax.set_facecolor(PALETTE["panel"])
    sca = ax.scatter(
        frame["UMAP1"],
        frame["UMAP2"],
        c=frame["cluster"],
        cmap=cmap,
        vmin=0.5,
        vmax=max_cluster + 0.5,
        s=args.umap_point_size,
        linewidths=0,
        alpha=0.86,
        rasterized=True,
    )
    ax.set_title(f"{args.sample} {assignment_label(assignment)} RNA Scanpy UMAP", fontsize=14, fontweight="bold", color=PALETTE["text"], pad=12)
    ax.set_xlabel("UMAP1", fontsize=10, color=PALETTE["muted"])
    ax.set_ylabel("UMAP2", fontsize=10, color=PALETTE["muted"])
    ax.grid(True, color=PALETTE["grid"], linewidth=0.45, alpha=0.55)
    ax.tick_params(labelsize=8, colors=PALETTE["muted"], length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(sca, ax=ax, ticks=np.arange(1, max_cluster + 1), fraction=0.035, pad=0.02)
    cbar.ax.set_yticklabels([f"C{i}" for i in range(1, max_cluster + 1)])
    cbar.set_label("Leiden cluster", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    stem = umap_stem(assignment)
    return save_figure(fig, Path(args.outdir), stem, args.formats, args.dpi)


def main():
    args = parse_args()
    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]
    assignments = [item.strip() for item in args.assignments.split(",") if item.strip()]
    invalid = set(assignments) - {"total", "old", "new"}
    if invalid:
        raise ValueError(f"Unsupported assignments: {', '.join(sorted(invalid))}")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_dir, features, barcodes, barcode_meta = load_matrix_inputs(args.matrix_dir, args.sample)
    summary = {
        "sample": args.sample,
        "bin_size": args.bin_size,
        "min_bin_umis": args.min_bin_umis,
        "method": "scanpy",
        "coordinate_used_for_clustering": False,
        "assignments": {},
        "input_source": f"official_celatlas_bin{args.bin_size}" if args.celatlas_bin_dir else "custom_slam_matrix",
        "bin_size_um": args.bin_size_um if args.bin_size_um is not None else args.bin_size,
        "bin_size_pixels": args.bin_size_pixels if args.bin_size_pixels is not None else args.bin_size,
    }
    for assignment in assignments:
        print(f"[scanpy] aggregate {assignment}", flush=True)
        if args.celatlas_bin_dir:
            if assignment != "total":
                raise ValueError("official Celatlas binned input supports total assignment only")
            adata = load_official_bin(args.celatlas_bin_dir)
        else:
            matrix, obs, var = aggregate_to_bins(
                matrix_dir, args.sample, assignment, features, barcodes, barcode_meta, args.bin_size, args.min_bin_umis
            )
            import anndata as ad
            adata = ad.AnnData(X=matrix, obs=obs, var=var)
        print(f"[scanpy] analyze {assignment}: bins={adata.n_obs} genes={adata.n_vars}", flush=True)
        adata, metrics = run_scanpy(adata, args)
        frame = output_frame(adata, assignment)
        frame.to_csv(outdir / f"{spatial_stem(assignment, args.bin_size)}.tsv", sep="\t", index=False)
        frame.to_csv(outdir / f"{umap_stem(assignment)}.tsv", sep="\t", index=False)
        plot_spatial(frame, assignment, args)
        plot_umap(frame, assignment, args)
        if args.write_h5ad:
            adata.write_h5ad(outdir / f"scanpy_{assignment}_bin{args.bin_size}.h5ad", compression="gzip")
        metrics.update(
            {
                "n_bins": int(adata.n_obs),
                "n_genes_after_filter": int(adata.n_vars),
                "n_clusters": int(frame["cluster"].nunique()),
                "total_umis": int(frame["umi_count"].sum()),
            }
        )
        summary["assignments"][assignment] = metrics

    (outdir / "scanpy_spatial_umap_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print("[scanpy] done", flush=True)


if __name__ == "__main__":
    main()
