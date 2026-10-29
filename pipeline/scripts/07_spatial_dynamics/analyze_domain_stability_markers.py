#!/usr/bin/env python3
"""Validate expression-domain stability and calculate total-RNA domain markers."""
from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

from spatial_slam_scanpy_analysis import import_scanpy, load_official_bin


COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
]
SIGNED = LinearSegmentedColormap.from_list("signed", ["#3C5488", "#F7F7F7", "#E64B35"])
GENERAL = LinearSegmentedColormap.from_list("general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"])
INK = "#25313A"


def natural_domains(values):
    return sorted(pd.Series(values).astype(str).unique(), key=lambda value: int(value) if value.isdigit() else value)


def spatial_edges(bin_ids, positions_path):
    positions = pd.read_csv(
        positions_path, header=None,
        names=["bin_id", "in_tissue", "row", "col", "pxl_row", "pxl_col"],
        dtype={"bin_id": str},
    )
    positions = positions[positions["bin_id"].isin(set(map(str, bin_ids)))].copy()
    index = {str(bin_id): idx for idx, bin_id in enumerate(bin_ids)}
    grid = {(int(row.row), int(row.col)): index[str(row.bin_id)] for row in positions.itertuples(index=False)}
    edges = []
    for (row, col), source in grid.items():
        for neighbor in ((row + 1, col), (row, col + 1)):
            target = grid.get(neighbor)
            if target is not None:
                edges.append((source, target))
    return np.asarray(edges, dtype=int), grid


def fragmentation_metrics(labels, edges):
    labels = np.asarray(labels).astype(str)
    adjacency = [[] for _ in range(len(labels))]
    for source, target in edges:
        adjacency[source].append(target)
        adjacency[target].append(source)
    shared = float(np.mean(labels[edges[:, 0]] == labels[edges[:, 1]])) if len(edges) else np.nan
    rows = []
    for domain in natural_domains(labels):
        members = set(np.flatnonzero(labels == domain).tolist())
        remaining = set(members)
        sizes = []
        while remaining:
            start = remaining.pop()
            queue = deque([start])
            size = 0
            while queue:
                node = queue.popleft()
                size += 1
                for neighbor in adjacency[node]:
                    if neighbor in remaining and labels[neighbor] == domain:
                        remaining.remove(neighbor)
                        queue.append(neighbor)
            sizes.append(size)
        rows.append({
            "domain": domain,
            "n_bins": len(members),
            "connected_components": len(sizes),
            "largest_component_fraction": max(sizes) / len(members) if members else 0.0,
        })
    return shared, pd.DataFrame(rows)


def cluster_metrics(adata, labels, edges, sample_indices):
    shared, fragmentation = fragmentation_metrics(labels, edges)
    sampled_labels = np.asarray(labels)[sample_indices]
    silhouette = np.nan
    if len(np.unique(sampled_labels)) > 1:
        silhouette = float(silhouette_score(adata.obsm["X_pca"][sample_indices], sampled_labels))
    weights = fragmentation["n_bins"].to_numpy(dtype=float)
    weighted_lcc = float(np.average(fragmentation["largest_component_fraction"], weights=weights))
    return {
        "n_clusters": int(len(np.unique(labels))),
        "silhouette_pca": silhouette,
        "local_neighbor_agreement": shared,
        "total_connected_components": int(fragmentation["connected_components"].sum()),
        "weighted_largest_component_fraction": weighted_lcc,
    }, fragmentation


def prepare_expression(bin_dir, domains_path, hvg, n_pcs, neighbors, seed):
    sc = import_scanpy()
    adata = load_official_bin(bin_dir)
    feature_path = Path(bin_dir) / "filtered_feature_bc_matrix/features.tsv.gz"
    if not feature_path.exists():
        feature_path = feature_path.with_suffix("")
    feature_meta = pd.read_csv(feature_path, sep="\t", header=None, compression="infer", dtype=str)
    adata.var["gene_name"] = feature_meta.iloc[:, 1].to_numpy() if feature_meta.shape[1] > 1 else adata.var_names
    domains = pd.read_csv(domains_path, sep="\t", dtype={"bin_id": str, "cluster": str, "domain": str})
    if "domain" not in domains:
        domains = domains.rename(columns={"cluster": "domain"})
    labels = domains.set_index("bin_id")["domain"].reindex(adata.obs_names.astype(str))
    if labels.isna().any():
        raise ValueError("Domain table does not cover all official Celatlas bins")
    adata.obs["domain"] = pd.Categorical(labels.astype(str), categories=natural_domains(labels), ordered=True)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg, flavor="seurat")
    reduced = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(reduced, max_value=10)
    sc.tl.pca(reduced, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(reduced, n_neighbors=neighbors, n_pcs=n_pcs, random_state=seed)
    return sc, adata, reduced


def stability(sc, reduced, domains_path, positions_path, resolutions, seeds, final_resolution, outdir):
    provided = pd.read_csv(domains_path, sep="\t", dtype={"bin_id": str})
    if "domain" not in provided:
        provided = provided.rename(columns={"cluster": "domain"})
    provided = provided.set_index("bin_id").reindex(reduced.obs_names)
    baseline_provided = provided["domain"].astype(str).to_numpy()
    edges, _ = spatial_edges(reduced.obs_names.astype(str), positions_path)
    rng = np.random.default_rng(0)
    sample_indices = np.sort(rng.choice(reduced.n_obs, size=min(10_000, reduced.n_obs), replace=False))

    resolution_rows = []
    fragmentation_frames = []
    for resolution in resolutions:
        key = f"resolution_{resolution:g}"
        sc.tl.leiden(reduced, resolution=resolution, random_state=0, n_iterations=2, key_added=key)
        labels = reduced.obs[key].astype(str).to_numpy()
        metrics, fragmentation = cluster_metrics(reduced, labels, edges, sample_indices)
        resolution_rows.append({"resolution": resolution, "seed": 0, **metrics})
        fragmentation.insert(0, "resolution", resolution)
        fragmentation_frames.append(fragmentation)
    resolution_table = pd.DataFrame(resolution_rows)
    resolution_table.to_csv(outdir / "domain_resolution_sweep.tsv", sep="\t", index=False)
    pd.concat(fragmentation_frames, ignore_index=True).to_csv(
        outdir / "domain_spatial_fragmentation.tsv", sep="\t", index=False
    )

    seed_labels = {}
    seed_rows = []
    seed_size_rows = []
    for seed in seeds:
        key = f"seed_{seed}"
        sc.tl.leiden(reduced, resolution=final_resolution, random_state=seed, n_iterations=2, key_added=key)
        labels = reduced.obs[key].astype(str).to_numpy()
        seed_labels[seed] = labels
        for domain, count in Counter(labels).items():
            seed_size_rows.append({"seed": seed, "resolution": final_resolution, "domain": domain, "n_bins": count})
        metrics, _ = cluster_metrics(reduced, labels, edges, sample_indices)
        seed_rows.append({
            "seed": seed,
            "resolution": final_resolution,
            "ari_vs_seed0": adjusted_rand_score(seed_labels[seeds[0]], labels),
            "nmi_vs_seed0": normalized_mutual_info_score(seed_labels[seeds[0]], labels),
            "ari_vs_formal_domains": adjusted_rand_score(baseline_provided, labels),
            "nmi_vs_formal_domains": normalized_mutual_info_score(baseline_provided, labels),
            **metrics,
        })
    seed_table = pd.DataFrame(seed_rows)
    seed_table.to_csv(outdir / "domain_seed_robustness.tsv", sep="\t", index=False)
    pd.DataFrame(seed_size_rows).to_csv(outdir / "domain_seed_cluster_sizes.tsv", sep="\t", index=False)
    pair_rows = []
    for left in seeds:
        for right in seeds:
            pair_rows.append({
                "seed_a": left, "seed_b": right,
                "ari": adjusted_rand_score(seed_labels[left], seed_labels[right]),
                "nmi": normalized_mutual_info_score(seed_labels[left], seed_labels[right]),
            })
    pd.DataFrame(pair_rows).to_csv(outdir / "domain_seed_pairwise.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.6), constrained_layout=True, facecolor="white")
    axes[0].plot(resolution_table["resolution"], resolution_table["n_clusters"], marker="o", color=COLORS[3])
    axes[0].set_ylabel("Number of domains")
    axes[0].set_title("A  Resolution sweep", fontweight="bold", color=INK)
    axes[1].plot(resolution_table["resolution"], resolution_table["silhouette_pca"], marker="o", color=COLORS[2], label="PCA silhouette")
    axes[1].plot(resolution_table["resolution"], resolution_table["local_neighbor_agreement"], marker="o", color=COLORS[0], label="Spatial neighbor agreement")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].set_title("B  Expression separation and spatial continuity", fontweight="bold", color=INK)
    axes[2].bar(seed_table["seed"].astype(str), seed_table["ari_vs_seed0"], color=COLORS[1])
    axes[2].set_ylim(0, 1.02)
    axes[2].set_ylabel("ARI versus seed 0")
    axes[2].set_title("C  Leiden seed robustness", fontweight="bold", color=INK)
    for axis in axes:
        axis.set_xlabel("Resolution" if axis is not axes[2] else "Random seed")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.55)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"domain_stability.{suffix}", dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return resolution_table, seed_table


def marker_analysis(sc, adata, outdir, top_per_domain):
    adata.var[["gene_name"]].rename_axis("gene_id").reset_index().to_csv(
        outdir / "domain_marker_gene_universe.tsv", sep="\t", index=False
    )
    sc.tl.rank_genes_groups(
        adata, "domain", method="wilcoxon", pts=True,
        n_genes=min(1000, adata.n_vars), use_raw=False,
    )
    pts = adata.uns["rank_genes_groups"]["pts"]
    pts_rest = adata.uns["rank_genes_groups"]["pts_rest"]
    name_map = adata.var["gene_name"].astype(str).to_dict()
    rows = []
    for domain in adata.obs["domain"].cat.categories:
        frame = sc.get.rank_genes_groups_df(adata, group=domain)
        frame = frame.rename(columns={
            "names": "gene_id", "logfoldchanges": "logfoldchange",
            "pvals": "pvalue", "pvals_adj": "adjusted_pvalue",
        })
        frame["domain"] = str(domain)
        frame["gene_name"] = frame["gene_id"].map(name_map).fillna(frame["gene_id"])
        frame["pct_in_domain"] = frame["gene_id"].map(pts[str(domain)])
        frame["pct_outside_domain"] = frame["gene_id"].map(pts_rest[str(domain)])
        frame["pct_difference"] = frame["pct_in_domain"] - frame["pct_outside_domain"]
        frame["rank"] = np.arange(1, len(frame) + 1)
        rows.append(frame[[
            "domain", "rank", "gene_id", "gene_name", "scores", "logfoldchange", "pvalue",
            "adjusted_pvalue", "pct_in_domain", "pct_outside_domain", "pct_difference",
        ]].rename(columns={"scores": "score"}).head(1000))
    markers = pd.concat(rows, ignore_index=True)
    markers.to_csv(outdir / "domain_markers.tsv", sep="\t", index=False)
    eligible = markers[
        (markers["adjusted_pvalue"] < 0.05)
        & (markers["logfoldchange"] > 0.5)
        & (markers["pct_in_domain"] >= 0.05)
        & (markers["pct_difference"] >= 0.05)
    ].copy()
    selected = []
    used = set()
    for domain in adata.obs["domain"].cat.categories:
        subset = eligible[eligible["domain"] == str(domain)].sort_values(["score", "pct_difference"], ascending=False)
        for row in subset.itertuples(index=False):
            if row.gene_id in used:
                continue
            selected.append(row)
            used.add(row.gene_id)
            if sum(item.domain == str(domain) for item in selected) >= top_per_domain:
                break
    selected_frame = pd.DataFrame(selected)
    selected_frame.to_csv(outdir / "domain_marker_display_selection.tsv", sep="\t", index=False)
    plot_marker_panels(adata, selected_frame, outdir)
    return markers, selected_frame


def plot_marker_panels(adata, selected, outdir):
    if selected.empty:
        return
    domains = list(adata.obs["domain"].cat.categories)
    genes = selected["gene_id"].astype(str).tolist()
    labels = selected["gene_name"].astype(str).tolist()
    matrix = adata[:, genes].X
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    means = np.zeros((len(domains), len(genes)), dtype=float)
    fractions = np.zeros_like(means)
    for index, domain in enumerate(domains):
        mask = np.asarray(adata.obs["domain"].astype(str) == str(domain))
        sub = matrix[mask]
        means[index] = np.asarray(sub.mean(axis=0)).ravel()
        fractions[index] = np.asarray((sub > 0).mean(axis=0)).ravel()

    fig_width = max(12.0, 0.42 * len(genes) + 4.0)
    fig, ax = plt.subplots(figsize=(fig_width, 5.8), facecolor="white")
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(len(domains)))
    scatter = ax.scatter(
        x.ravel(), y.ravel(), s=15 + 180 * fractions.ravel(), c=means.ravel(),
        cmap=GENERAL, linewidths=0,
    )
    ax.set_xticks(np.arange(len(genes)), labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(domains)), [f"D{value}" for value in domains])
    ax.invert_yaxis()
    ax.set_title("Domain marker dotplot", fontweight="bold", color=INK)
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("Expression-defined domain")
    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Mean log-normalized expression", fontsize=8)
    for fraction in (0.1, 0.3, 0.6):
        ax.scatter([], [], s=15 + 180 * fraction, color="#8491B4", label=f"{fraction:.0%}")
    ax.legend(title="Detected bins", bbox_to_anchor=(1.08, 1), loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"domain_marker_dotplot.{suffix}", dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    gene_means = means.mean(axis=0)
    gene_sd = means.std(axis=0)
    z = (means - gene_means) / np.where(gene_sd > 0, gene_sd, 1.0)
    limit = max(float(np.nanmax(np.abs(z))), 1.0)
    fig, ax = plt.subplots(figsize=(fig_width, 5.8), facecolor="white")
    image = ax.imshow(z, aspect="auto", cmap=SIGNED, vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(genes)), labels, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(domains)), [f"D{value}" for value in domains])
    ax.set_title("Domain marker heatmap", fontweight="bold", color=INK)
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("Expression-defined domain")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Domain-centered mean expression z-score", fontsize=8)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"domain_marker_heatmap.{suffix}", dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--celatlas-bin-dir", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--resolutions", default="0.4,0.6,0.8,1.0,1.2,1.4")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--final-resolution", type=float, default=1.0)
    parser.add_argument("--hvg", type=int, default=3000)
    parser.add_argument("--pca-components", type=int, default=25)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-markers-per-domain", type=int, default=3)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sc, full, reduced = prepare_expression(
        args.celatlas_bin_dir, args.domains, args.hvg, args.pca_components, args.neighbors, args.seed
    )
    positions = Path(args.celatlas_bin_dir) / "spatial/tissue_positions_list.csv"
    resolution_table, seed_table = stability(
        sc, reduced, args.domains, positions,
        [float(value) for value in args.resolutions.split(",")],
        [int(value) for value in args.seeds.split(",")],
        args.final_resolution, outdir,
    )
    markers, selected = marker_analysis(sc, full, outdir, args.top_markers_per_domain)
    summary = {
        "n_bins": int(full.n_obs),
        "n_genes_tested": int(full.n_vars),
        "marker_method": "scanpy_rank_genes_groups_wilcoxon",
        "marker_rows": int(len(markers)),
        "display_markers": int(len(selected)),
        "resolutions": resolution_table["resolution"].tolist(),
        "seeds": seed_table["seed"].tolist(),
        "mean_seed_ari_vs_seed0": float(seed_table["ari_vs_seed0"].mean()),
        "minimum_seed_ari_vs_seed0": float(seed_table["ari_vs_seed0"].min()),
        "formal_resolution": args.final_resolution,
        "coordinate_used_for_clustering": False,
    }
    (outdir / "domain_stability_marker_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
