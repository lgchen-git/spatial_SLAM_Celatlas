#!/usr/bin/env python3
"""Build an expression-derived tissue mask from smoothed total-RNA bin counts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from skimage.filters import threshold_multiotsu
from skimage.measure import label
from skimage.morphology import binary_closing, disk, remove_small_objects


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-summary", required=True)
    parser.add_argument("--barcode-metadata", required=True)
    parser.add_argument("--output-dir", required=True, help="Stage 07 result directory")
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--coordinate-um", type=float, default=0.5)
    parser.add_argument("--smooth-sigma-bins", type=float, default=2.0)
    parser.add_argument("--multiotsu-classes", type=int, default=4)
    parser.add_argument("--closing-radius-bins", type=int, default=2)
    parser.add_argument("--min-component-bins", type=int, default=50)
    parser.add_argument("--dpi", type=int, default=450)
    return parser.parse_args()


def grid_from_bins(frame):
    xs = np.sort(frame["bin_x"].unique())
    ys = np.sort(frame["bin_y"].unique())
    x_index = {value: idx for idx, value in enumerate(xs)}
    y_index = {value: idx for idx, value in enumerate(ys)}
    grid = np.full((len(ys), len(xs)), np.nan, dtype=float)
    for row in frame.itertuples(index=False):
        grid[y_index[row.bin_y], x_index[row.bin_x]] = np.log1p(row.total)
    return grid, xs, ys


def build_mask(grid, args):
    occupied = np.isfinite(grid)
    values = np.nan_to_num(grid, nan=0.0)
    weight = gaussian_filter(occupied.astype(float), args.smooth_sigma_bins)
    smoothed = gaussian_filter(values, args.smooth_sigma_bins) / np.maximum(weight, 1e-8)
    thresholds = threshold_multiotsu(smoothed[occupied], classes=args.multiotsu_classes)
    threshold = float(thresholds[-1])
    mask = occupied & (smoothed >= threshold)
    mask = binary_closing(mask, disk(args.closing_radius_bins))
    mask = remove_small_objects(mask, min_size=args.min_component_bins)
    components = label(mask)
    component_sizes = sorted(np.bincount(components.ravel())[1:].tolist(), reverse=True)
    return mask, smoothed, thresholds.tolist(), component_sizes


def save_qc_figure(grid, mask, xs, ys, figures, dpi):
    figures.mkdir(parents=True, exist_ok=True)
    extent = [float(xs.min() - 25), float(xs.max() + 25), float(ys.max() + 25), float(ys.min() - 25)]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.8), facecolor="white", constrained_layout=True)
    image = axes[0].imshow(grid, origin="upper", extent=extent, cmap="viridis", aspect="equal")
    axes[0].set_title("Total RNA density (log1p)", fontweight="bold")
    colorbar = fig.colorbar(image, ax=axes[0], fraction=0.04, pad=0.02)
    colorbar.set_label("log1p molecules per raw bin50")
    axes[1].imshow(mask, origin="upper", extent=extent, cmap="gray_r", vmin=0, vmax=1, aspect="equal")
    axes[1].set_title("Expression-derived tissue mask", fontweight="bold")
    for axis in axes:
        axis.set_xlabel("Spatial X (raw chip units)")
        axis.set_ylabel("Spatial Y (raw chip units; downward)")
        for spine in axis.spines.values():
            spine.set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"01_tissue_mask_expression_derived.{suffix}", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    args = parse_args()
    output = Path(args.output_dir)
    tables, figures, metadata = output / "01_tables", output / "02_figures", output / "05_metadata"
    for path in (tables, figures, metadata):
        path.mkdir(parents=True, exist_ok=True)
    bins = pd.read_csv(args.bin_summary, sep="\t")
    required = {"bin_id", "bin_x", "bin_y", "total"}
    missing = required - set(bins.columns)
    if missing:
        raise ValueError(f"bin summary is missing: {sorted(missing)}")
    grid, xs, ys = grid_from_bins(bins)
    mask, smoothed, thresholds, component_sizes = build_mask(grid, args)
    mask_lookup = {(float(xs[x]), float(ys[y])) for y, x in np.argwhere(mask)}
    bins = bins.copy()
    bins["is_tissue"] = [(float(x), float(y)) in mask_lookup for x, y in zip(bins["bin_x"], bins["bin_y"])]
    bins.to_csv(tables / "01_tissue_bins.tsv", sep="\t", index=False)

    barcodes = pd.read_csv(args.barcode_metadata, sep="\t")
    barcodes["bin_x"] = (barcodes["x"] // args.bin_size).astype(int) * args.bin_size + args.bin_size / 2.0
    barcodes["bin_y"] = (barcodes["y"] // args.bin_size).astype(int) * args.bin_size + args.bin_size / 2.0
    keep = [(float(x), float(y)) in mask_lookup for x, y in zip(barcodes["bin_x"], barcodes["bin_y"])]
    tissue = barcodes.loc[keep, ["barcode", "x", "y", "bin_x", "bin_y"]].copy()
    tissue.to_csv(tables / "01_tissue_barcodes.tsv", sep="\t", index=False)
    save_qc_figure(grid, mask, xs, ys, figures, args.dpi)

    selected = bins[bins["is_tissue"]]
    qc = {
        "method": "expression_derived_smoothed_log_total_multiotsu",
        "image_classification": "expression-derived; no histology or microscopy input",
        "source_bin_summary": str(Path(args.bin_summary).resolve()),
        "source_barcode_metadata": str(Path(args.barcode_metadata).resolve()),
        "bin_size_raw": args.bin_size,
        "coordinate_um_per_raw_unit": args.coordinate_um,
        "bin_size_um": args.bin_size * args.coordinate_um,
        "canonical_orientation": "X increases right; Y increases downward",
        "smoothing_sigma_bins": args.smooth_sigma_bins,
        "multiotsu_classes": args.multiotsu_classes,
        "multiotsu_thresholds_log1p": thresholds,
        "selected_threshold_log1p": thresholds[-1],
        "closing_radius_bins": args.closing_radius_bins,
        "min_component_bins": args.min_component_bins,
        "component_sizes_bins": component_sizes,
        "n_input_bins": int(len(bins)),
        "n_tissue_bins": int(mask.sum()),
        "tissue_bin_fraction": float(mask.sum() / max(np.isfinite(grid).sum(), 1)),
        "n_input_barcodes": int(len(barcodes)),
        "n_tissue_barcodes": int(len(tissue)),
        "tissue_barcode_fraction": float(len(tissue) / max(len(barcodes), 1)),
        "tissue_total_molecules": float(selected["total"].sum()),
        "interpretive_limit": "The mask is inferred from RNA density and is not histology or staining evidence. A manual lasso GeoJSON can override it.",
    }
    (metadata / "tissue_mask_qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
