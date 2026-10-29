#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


PALETTE = {
    "paper": "#f7f5ef",
    "panel": "#ffffff",
    "text": "#26333c",
    "muted": "#65737c",
    "grid": "#d9dde1",
    "accent": "#b85c38",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Build single-sample Spatial SLAM kinetics summaries from corrected counts.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--corrected-bin-gene", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--label-hours", type=float)
    parser.add_argument("--min-gene-total-umis", type=int, default=100)
    parser.add_argument("--top-genes-per-domain", type=int, default=8)
    parser.add_argument("--smoothing-neighbors", type=int, default=24)
    parser.add_argument("--smoothing-max-radius", type=float)
    parser.add_argument("--official-bin-positions")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--formats", default="png,pdf")
    return parser.parse_args()


def safe_divide(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator, dtype=float), where=denominator > 0)


def add_kinetic_columns(frame, label_hours):
    frame = frame.copy()
    frame["ntr"] = safe_divide(frame["expected_new"], frame["total"])
    if "posterior_variance" in frame:
        frame["ntr_se"] = safe_divide(np.sqrt(np.maximum(frame["posterior_variance"], 0.0)), frame["total"])
        frame["ntr_ci_lower"] = np.maximum(0.0, frame["ntr"] - 1.96 * frame["ntr_se"])
        frame["ntr_ci_upper"] = np.minimum(1.0, frame["ntr"] + 1.96 * frame["ntr_se"])
    if label_hours and label_hours > 0:
        clipped = np.clip(frame["ntr"].to_numpy(dtype=float), 1e-8, 0.999999)
        gamma = -np.log1p(-clipped) / float(label_hours)
        frame["gamma_per_hour"] = gamma
        frame["half_life_hours"] = np.log(2.0) / np.maximum(gamma, 1e-12)
        frame["alpha_relative"] = gamma * np.log1p(frame["total"].to_numpy(dtype=float))
    return frame


def save(fig, outdir, stem, formats, dpi):
    paths = []
    for fmt in formats:
        path = outdir / f"{stem}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        paths.append(str(path))
    plt.close(fig)
    return paths


def load_domains(path):
    domains = pd.read_csv(path, sep="\t")
    if "domain" not in domains.columns:
        domains = domains.rename(columns={"cluster": "domain"})
    domains["bin_id"] = domains["bin_id"].astype(str)
    domains["domain"] = domains["domain"].astype(str)
    return domains[["bin_id", "bin_x", "bin_y", "domain", "UMAP1", "UMAP2"]]


def build_domain_gene(corrected, domains, min_gene_total_umis, label_hours):
    merged = corrected.merge(domains[["bin_id", "domain"]], on="bin_id", how="inner")
    aggregations = {
        "total": ("total", "sum"), "expected_new": ("expected_new", "sum"),
        "expected_old": ("expected_old", "sum"), "n_bins": ("bin_id", "nunique"),
    }
    if "posterior_variance" in merged:
        aggregations["posterior_variance"] = ("posterior_variance", "sum")
    grouped = merged.groupby(["domain", "gene_id", "gene_name"], as_index=False).agg(**aggregations)
    gene_global = (
        grouped.groupby(["gene_id", "gene_name"], as_index=False)
        .agg(global_total=("total", "sum"), global_expected_new=("expected_new", "sum"))
    )
    gene_global["global_ntr"] = safe_divide(gene_global["global_expected_new"], gene_global["global_total"])
    grouped = grouped[grouped["total"] >= min_gene_total_umis].copy()
    grouped = add_kinetic_columns(grouped, label_hours)
    grouped = grouped.merge(gene_global, on=["gene_id", "gene_name"], how="left")
    grouped["delta_ntr"] = grouped["ntr"] - grouped["global_ntr"]
    grouped["dynamic_score"] = grouped["delta_ntr"] * np.log1p(grouped["total"])
    return grouped.sort_values(["domain", "dynamic_score"], ascending=[True, False])


def select_heatmap_genes(domain_gene, top_genes_per_domain):
    chosen = []
    for _, sub in domain_gene.groupby("domain", sort=True):
        ranked = sub.assign(abs_dynamic_score=sub["dynamic_score"].abs()).sort_values("abs_dynamic_score", ascending=False)
        chosen.extend(ranked.head(top_genes_per_domain)["gene_id"].tolist())
    seen = set()
    genes = []
    for gene in chosen:
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return genes[:60]


def select_representative_display_genes(domain_gene):
    """Return one uncertainty-supported display gene per domain.

    The reviewed biology report selects representatives after this stage.  The
    per-bin export must therefore include the same class of candidates, or a
    later report build can produce an empty panel for an otherwise valid gene.
    This deliberately uses only quantities already available in domain_gene;
    marker annotations and story labels remain report-layer decisions.
    """
    frame = domain_gene.copy()
    if frame.empty or "ntr_ci_lower" not in frame or "ntr_ci_upper" not in frame:
        return []
    frame["relative_uncertainty"] = (frame["ntr_ci_upper"] - frame["ntr_ci_lower"]).fillna(1.0)
    frame["display_score"] = (
        frame["dynamic_score"].abs() * np.log1p(frame["n_bins"]) /
        (1.0 + frame["relative_uncertainty"])
    )
    eligible = frame[
        (frame["total"] >= 500) &
        (frame["n_bins"] >= 10) &
        (frame["relative_uncertainty"] <= 0.25)
    ]
    chosen = []
    for _, sub in eligible.sort_values("display_score", ascending=False).groupby("domain", sort=True):
        chosen.append(str(sub.iloc[0]["gene_id"]))
    return chosen


def plot_domain_heatmap(domain_gene, value_col, stem, title, outdir, args, cmap="coolwarm"):
    genes = select_heatmap_genes(domain_gene, args.top_genes_per_domain)
    if not genes:
        return []
    frame = domain_gene[domain_gene["gene_id"].isin(genes)].copy()
    labels = frame.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].to_dict()
    pivot = frame.pivot_table(index="gene_id", columns="domain", values=value_col, aggfunc="mean")
    pivot = pivot.reindex(genes).fillna(0.0)
    fig_h = max(5.0, min(13.0, 0.24 * len(pivot) + 2.0))
    fig, ax = plt.subplots(figsize=(8.8, fig_h), facecolor=PALETTE["paper"])
    ax.set_facecolor(PALETTE["panel"])
    values = pivot.to_numpy(dtype=float)
    if value_col == "delta_ntr":
        vmax = max(float(np.nanmax(np.abs(values))), 0.05)
        im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    else:
        im = ax.imshow(values, aspect="auto", cmap=cmap)
    ax.set_title(title, fontsize=13, fontweight="bold", color=PALETTE["text"], pad=10)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([f"D{d}" for d in pivot.columns], fontsize=8)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels([labels.get(gene, gene) for gene in pivot.index], fontsize=7)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.ax.tick_params(labelsize=7)
    return save(fig, outdir, stem, args.formats, args.dpi)


def add_tissue_components(domains, positions_path):
    result = domains.copy()
    result["tissue_component"] = "all"
    if not positions_path:
        return result
    positions = pd.read_csv(
        positions_path, header=None,
        names=["bin_id", "in_tissue", "row", "col", "pxl_row", "pxl_col"],
        dtype={"bin_id": str},
    )
    positions = positions[positions["bin_id"].isin(result["bin_id"])].copy()
    coordinate_to_bin = {(int(row.row), int(row.col)): str(row.bin_id) for row in positions.itertuples(index=False)}
    component = {}
    component_id = 0
    for coordinate, bin_id in coordinate_to_bin.items():
        if bin_id in component:
            continue
        component_id += 1
        stack = [coordinate]
        component[bin_id] = component_id
        while stack:
            row, col = stack.pop()
            for neighbor in ((row - 1, col), (row + 1, col), (row, col - 1), (row, col + 1)):
                neighbor_bin = coordinate_to_bin.get(neighbor)
                if neighbor_bin is not None and neighbor_bin not in component:
                    component[neighbor_bin] = component_id
                    stack.append(neighbor)
    result["tissue_component"] = result["bin_id"].map(component).fillna(-1).astype(str)
    return result


def median_grid_spacing(positions_path):
    if not positions_path:
        return None
    positions = pd.read_csv(
        positions_path, header=None,
        names=["bin_id", "in_tissue", "row", "col", "pxl_row", "pxl_col"],
    )
    grid = {(int(row.row), int(row.col)): (float(row.pxl_col), float(row.pxl_row))
            for row in positions.itertuples(index=False)}
    distances = []
    for (row, col), (x, y) in grid.items():
        for neighbor in ((row + 1, col), (row, col + 1)):
            if neighbor in grid:
                x2, y2 = grid[neighbor]
                distances.append(float(np.hypot(x - x2, y - y2)))
    return float(np.median(distances)) if distances else None


def smooth_bin_gene(corrected, domains, genes, neighbors, label_hours, max_radius=None):
    selected = corrected[corrected["gene_id"].isin(genes)].merge(
        domains[["bin_id", "domain", "UMAP1", "UMAP2", "tissue_component"]],
        on="bin_id",
        how="inner",
    )
    rows = []
    for (gene_id, tissue_component), sub in selected.groupby(["gene_id", "tissue_component"], sort=False):
        # Smoothing must follow physical tissue neighborhoods, not an expression
        # embedding whose geometry and orientation have no physical meaning.
        coords = sub[["bin_x", "bin_y"]].to_numpy(dtype=float)
        if sub.empty:
            continue
        k = max(1, min(neighbors, len(sub)))
        tree = cKDTree(coords)
        distances, indices = tree.query(coords, k=k, distance_upper_bound=max_radius or np.inf)
        total = sub["total"].to_numpy(dtype=float)
        expected_new = sub["expected_new"].to_numpy(dtype=float)
        for idx, row in enumerate(sub.itertuples(index=False)):
            candidate_indices = np.atleast_1d(indices[idx])
            candidate_distances = np.atleast_1d(distances[idx])
            neigh = candidate_indices[(candidate_indices < len(sub)) & np.isfinite(candidate_distances)]
            if len(neigh) == 0:
                neigh = np.array([idx], dtype=int)
            sm_total = float(total[neigh].sum())
            sm_new = float(expected_new[neigh].sum())
            ntr = sm_new / sm_total if sm_total else 0.0
            record = {
                "bin_id": row.bin_id,
                "bin_x": row.bin_x,
                "bin_y": row.bin_y,
                "domain": row.domain,
                "tissue_component": tissue_component,
                "gene_id": gene_id,
                "gene_name": row.gene_name,
                "smoothed_total": sm_total,
                "smoothed_expected_new": sm_new,
                "smoothed_ntr": ntr,
                "smoothing_neighbors_used": int(len(neigh)),
                "smoothing_max_radius": max_radius,
            }
            if label_hours and label_hours > 0:
                gamma = -np.log1p(-np.clip(ntr, 1e-8, 0.999999)) / float(label_hours)
                record["smoothed_gamma_per_hour"] = gamma
                record["smoothed_half_life_hours"] = np.log(2.0) / max(gamma, 1e-12)
            rows.append(record)
    return pd.DataFrame(rows)


def plot_representative_gene_map(smoothed, outdir, args):
    if smoothed.empty:
        return []
    gene_scores = (
        smoothed.groupby(["gene_id", "gene_name"], as_index=False)
        .agg(mean_ntr=("smoothed_ntr", "mean"), total=("smoothed_total", "sum"))
        .sort_values(["mean_ntr", "total"], ascending=[False, False])
    )
    top = gene_scores.head(4)
    if top.empty:
        return []
    fig, axes = plt.subplots(1, len(top), figsize=(4.2 * len(top), 4.2), facecolor=PALETTE["paper"])
    axes = np.atleast_1d(axes)
    for ax, row in zip(axes, top.itertuples(index=False)):
        sub = smoothed[smoothed["gene_id"] == row.gene_id]
        ax.set_facecolor(PALETTE["panel"])
        sc = ax.scatter(sub["bin_x"], sub["bin_y"], c=sub["smoothed_ntr"], s=5, cmap="magma", vmin=0, vmax=1, linewidths=0)
        ax.set_title(str(row.gene_name), fontsize=11, fontweight="bold", color=PALETTE["text"])
        ax.set_aspect("equal", adjustable="box")
        ax.invert_yaxis()
        ax.tick_params(labelsize=6, colors=PALETTE["muted"], length=2)
        for spine in ax.spines.values():
            spine.set_visible(False)
    cbar = fig.colorbar(sc, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("smoothed NTR", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    return save(fig, outdir, "14_representative_gene_spatial_ntr", args.formats, args.dpi)


def main():
    args = parse_args()
    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    corrected_path = Path(args.corrected_bin_gene)
    corrected = pd.read_pickle(corrected_path) if corrected_path.suffix == ".pkl" else pd.read_csv(corrected_path, sep="\t")
    corrected["bin_id"] = corrected["bin_id"].astype(str)
    domains = add_tissue_components(load_domains(args.domains), args.official_bin_positions)
    grid_spacing = median_grid_spacing(args.official_bin_positions)
    bin_aggregations = {
        "total": ("total", "sum"),
        "expected_new": ("expected_new", "sum"),
        "expected_old": ("expected_old", "sum"),
    }
    if "posterior_variance" in corrected:
        bin_aggregations["posterior_variance"] = ("posterior_variance", "sum")
    bin_summary = corrected.groupby(["bin_id", "bin_x", "bin_y"], as_index=False).agg(**bin_aggregations)
    bin_summary["corrected_ntr"] = safe_divide(bin_summary["expected_new"], bin_summary["total"])
    if "posterior_variance" in bin_summary:
        bin_summary["ntr_se"] = safe_divide(np.sqrt(np.maximum(bin_summary["posterior_variance"], 0.0)), bin_summary["total"])
        bin_summary["ntr_ci_lower"] = np.maximum(0.0, bin_summary["corrected_ntr"] - 1.96 * bin_summary["ntr_se"])
        bin_summary["ntr_ci_upper"] = np.minimum(1.0, bin_summary["corrected_ntr"] + 1.96 * bin_summary["ntr_se"])
    bin_summary.to_csv(
        outdir / "bin_kinetics_summary.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "compresslevel": 1},
    )

    domain_gene = build_domain_gene(corrected, domains, args.min_gene_total_umis, args.label_hours)
    domain_gene.to_csv(outdir / "domain_gene_kinetics.tsv", sep="\t", index=False)
    global_aggregations = {
        "total": ("total", "sum"), "expected_new": ("expected_new", "sum"),
        "expected_old": ("expected_old", "sum"), "n_bins": ("bin_id", "nunique"),
    }
    if "posterior_variance" in corrected:
        global_aggregations["posterior_variance"] = ("posterior_variance", "sum")
    global_gene = corrected.groupby(["gene_id", "gene_name"], as_index=False).agg(**global_aggregations)
    global_gene = add_kinetic_columns(global_gene, args.label_hours)
    global_gene.to_csv(outdir / "global_gene_kinetics.tsv", sep="\t", index=False)
    plot_domain_heatmap(domain_gene, "delta_ntr", "12_domain_ntr_delta_heatmap", "Domain gene NTR enrichment", outdir, args)
    if args.label_hours and args.label_hours > 0 and "half_life_hours" in domain_gene.columns:
        plot_domain_heatmap(
            domain_gene,
            "half_life_hours",
            "13_domain_half_life_heatmap",
            "Domain gene half-life estimates",
            outdir,
            args,
            cmap="viridis",
        )

    smooth_genes = select_heatmap_genes(domain_gene, args.top_genes_per_domain)[:40]
    for gene in select_representative_display_genes(domain_gene):
        if gene not in smooth_genes:
            smooth_genes.append(gene)
    raw_selected = corrected[corrected["gene_id"].isin(smooth_genes)].merge(
        domains[["bin_id", "domain", "tissue_component"]], on="bin_id", how="inner"
    )
    raw_selected.to_csv(
        outdir / "bin_gene_kinetics_raw.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "compresslevel": 1},
    )
    smoothed = smooth_bin_gene(
        corrected, domains, smooth_genes, args.smoothing_neighbors, args.label_hours,
        args.smoothing_max_radius,
    )
    smoothed.to_csv(
        outdir / "bin_gene_kinetics_smoothed.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "compresslevel": 1},
    )
    plot_representative_gene_map(smoothed, outdir, args)

    summary = {
        "sample": args.sample,
        "mode": "single_sample_corrected_ntr",
        "label_hours": args.label_hours,
        "n_domain_gene_rows": int(len(domain_gene)),
        "n_smoothed_rows": int(len(smoothed)),
        "n_bins": int(len(bin_summary)),
        "smoothing_neighbors": args.smoothing_neighbors,
        "smoothing_max_radius": args.smoothing_max_radius,
        "median_four_neighbor_grid_spacing": grid_spacing,
        "smoothing_radius_in_bin_spacings": (
            args.smoothing_max_radius / grid_spacing
            if args.smoothing_max_radius is not None and grid_spacing not in (None, 0)
            else None
        ),
        "tissue_components": int(domains["tissue_component"].nunique()),
        "half_life_reported": bool(args.label_hours and args.label_hours > 0),
        "domains": sorted(domain_gene["domain"].astype(str).unique().tolist()) if not domain_gene.empty else [],
    }
    (outdir / "kinetics_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
