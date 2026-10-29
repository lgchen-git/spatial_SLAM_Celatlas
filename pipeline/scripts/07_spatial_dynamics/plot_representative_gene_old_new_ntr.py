#!/usr/bin/env python3
"""Plot matched old RNA, new RNA, and NTR maps for dynamic genes at one bin scale."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


NPG_GENERAL = LinearSegmentedColormap.from_list(
    "npg_general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"]
)
NPG_ENRICHMENT = LinearSegmentedColormap.from_list(
    "npg_enrichment", ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"]
)
INK = "#25313A"
MUTED = "#64727B"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoothed-bin-gene", required=True)
    parser.add_argument("--domain-gene-kinetics", required=True)
    parser.add_argument("--domains", required=True, help="Total-RNA domain table used as the spatial footprint")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--metadata-dir", help="Optional stage metadata directory for provenance and manifest")
    parser.add_argument("--genes", help="Optional comma-separated gene symbols or IDs")
    parser.add_argument("--n-genes", type=int, default=4)
    parser.add_argument("--min-total", type=float, default=1000.0)
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--background-point-size", type=float, default=8.0)
    parser.add_argument("--point-size", type=float, default=12.0)
    parser.add_argument("--interpretation-mode", choices=("standard", "unlabeled_control"), default="standard")
    parser.add_argument("--bin-label", default="bin50")
    return parser.parse_args()


def select_genes(smoothed, kinetics, genes, n_genes, min_total):
    if genes:
        requested = [item.strip() for item in genes.split(",") if item.strip()]
        selected = []
        for item in requested:
            match = smoothed[(smoothed["gene_id"].astype(str) == item) | (smoothed["gene_name"].astype(str) == item)]
            if not match.empty:
                selected.append(str(match.iloc[0]["gene_id"]))
        if not selected:
            raise ValueError("None of the requested genes are available in the smoothed bin-gene table")
        return selected[:n_genes]
    ranked = (
        kinetics.assign(abs_delta=kinetics["delta_ntr"].abs())
        .groupby(["gene_id", "gene_name"], as_index=False)
        .agg(dynamic_score=("abs_delta", "max"), total=("total", "sum"))
        .sort_values(["dynamic_score", "total"], ascending=False)
    )
    available = set(smoothed["gene_id"].astype(str))
    ranked = ranked[(ranked["total"] >= min_total) & ranked["gene_id"].astype(str).isin(available)]
    return ranked.head(n_genes)["gene_id"].astype(str).tolist()


def style_spatial(ax, row, column):
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    if row == 0:
        ax.set_title(column, fontsize=11, fontweight="bold", color=INK, pad=7)
    ax.tick_params(labelsize=6.5, colors=MUTED, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_map(ax, footprint, frame, values, cmap, vmin, vmax, row, column, background_point_size, point_size):
    ax.scatter(
        footprint["bin_x"], footprint["bin_y"], marker="s", s=background_point_size,
        color="#EEF1F3", linewidths=0, rasterized=True,
    )
    order = np.argsort(values)
    scatter = ax.scatter(
        frame.iloc[order]["bin_x"], frame.iloc[order]["bin_y"], c=values[order],
        marker="s", s=point_size, cmap=cmap, vmin=vmin, vmax=vmax,
        linewidths=0, rasterized=True,
    )
    style_spatial(ax, row, column)
    return scatter


def main():
    args = parse_args()
    control = args.interpretation_mode == "unlabeled_control"
    smoothed = pd.read_csv(args.smoothed_bin_gene, sep="\t", compression="infer")
    kinetics = pd.read_csv(args.domain_gene_kinetics, sep="\t")
    footprint = pd.read_csv(args.domains, sep="\t")[["bin_x", "bin_y"]].drop_duplicates()
    selected = select_genes(smoothed, kinetics, args.genes, args.n_genes, args.min_total)
    if not selected:
        raise ValueError("No genes passed dynamic-score and abundance filters")

    nrows = len(selected)
    fig, axes = plt.subplots(
        nrows, 3, figsize=(13.2, 3.7 * nrows), facecolor="white",
        constrained_layout=True, squeeze=False,
    )
    gene_records = []
    for row, gene_id in enumerate(selected):
        frame = smoothed[smoothed["gene_id"].astype(str) == gene_id].copy()
        gene_name = str(frame.iloc[0]["gene_name"])
        frame["smoothed_expected_old"] = (frame["smoothed_total"] - frame["smoothed_expected_new"]).clip(lower=0)
        old = np.log1p(frame["smoothed_expected_old"].to_numpy(dtype=float))
        new = np.log1p(frame["smoothed_expected_new"].to_numpy(dtype=float))
        shared = np.concatenate([old, new])
        vmin, vmax = np.percentile(shared[np.isfinite(shared)], [2, 98])
        ntr = frame["smoothed_ntr"].to_numpy(dtype=float)
        old_scatter = draw_map(axes[row, 0], footprint, frame, old, NPG_GENERAL, vmin, vmax, row, "Background-old component" if control else "Old RNA", args.background_point_size, args.point_size)
        new_scatter = draw_map(axes[row, 1], footprint, frame, new, NPG_GENERAL, vmin, vmax, row, "Apparent new / T>C-positive" if control else "New RNA", args.background_point_size, args.point_size)
        ntr_scatter = draw_map(axes[row, 2], footprint, frame, ntr, NPG_ENRICHMENT, 0.0, 1.0, row, "Apparent NTR" if control else "New-to-total ratio (NTR)", args.background_point_size, args.point_size)
        axes[row, 0].set_ylabel(
            f"{gene_name}\nSpatial Y (downward)", fontsize=9.5, fontweight="bold", color=INK
        )
        axes[row, 1].set_ylabel("")
        axes[row, 2].set_ylabel("")
        if row == nrows - 1:
            for ax in axes[row]:
                ax.set_xlabel("Spatial X (raw chip units)", fontsize=7.5, color=MUTED)
        expression_cbar = fig.colorbar(new_scatter, ax=axes[row, :2], fraction=0.018, pad=0.012)
        expression_cbar.set_label("log1p smoothed expected molecules\n(shared old/new scale)", fontsize=7)
        expression_cbar.ax.tick_params(labelsize=6.5)
        ntr_cbar = fig.colorbar(ntr_scatter, ax=axes[row, 2], fraction=0.035, pad=0.012)
        ntr_cbar.set_label("Smoothed apparent NTR" if control else "Smoothed NTR", fontsize=7)
        ntr_cbar.ax.tick_params(labelsize=6.5)
        score = kinetics.loc[kinetics["gene_id"].astype(str) == gene_id, "delta_ntr"].abs().max()
        gene_records.append({
            "gene_id": gene_id, "gene_name": gene_name,
            "max_abs_domain_delta_ntr": float(score), "n_detected_bins": int(len(frame)),
            "old_new_display_log1p_limits": [float(vmin), float(vmax)],
        })

    fig.suptitle(
        (f"Representative background-high genes: component and apparent-NTR maps — raw {args.bin_label}"
         if control else f"Representative dynamic genes: matched old RNA, new RNA, and spatial NTR — raw {args.bin_label}"),
        fontsize=14, fontweight="bold", color=INK,
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stem = f"06_representative_dynamic_genes_old_new_ntr_{args.bin_label}"
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"{stem}.{suffix}", dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    metadata = {
        "figure": stem,
        "interpretation_mode": args.interpretation_mode,
        "classification": "expression-derived; no histology or microscopy background",
        "binning": f"raw {args.bin_label}",
        "canonical_orientation": "X increases right; Y increases downward",
        "old_new_comparison": (
            "background-old and apparent-new components use one shared log1p display scale within each gene"
            if control else "old and new use one shared log1p display scale within each gene"
        ),
        "ntr_scale": [0.0, 1.0],
        "smoothing": "physical-neighbor aggregation inherited from spatial_slam_kinetics.py",
        "genes": gene_records,
        "sources": {
            "smoothed_bin_gene": str(Path(args.smoothed_bin_gene).resolve()),
            "domain_gene_kinetics": str(Path(args.domain_gene_kinetics).resolve()),
            "spatial_footprint": str(Path(args.domains).resolve()),
        },
    }
    metadata_dir = Path(args.metadata_dir) if args.metadata_dir else output
    metadata_dir.mkdir(parents=True, exist_ok=True)
    (metadata_dir / f"{stem}.provenance.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    manifest_path = metadata_dir / "figure_manifest.tsv"
    manifest_row = {
        "stem": stem,
        "title": ("Representative background-high genes: components and apparent NTR"
                  if control else "Representative dynamic genes: matched old RNA, new RNA, and NTR"),
        "png": f"02_figures/{stem}.png",
        "pdf": f"02_figures/{stem}.pdf",
        "classification": "expression-derived",
    }
    if manifest_path.is_file():
        manifest = pd.read_csv(manifest_path, sep="\t")
        manifest = manifest[manifest["stem"] != stem]
        manifest = pd.concat([manifest, pd.DataFrame([manifest_row])], ignore_index=True)
    else:
        manifest = pd.DataFrame([manifest_row])
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
