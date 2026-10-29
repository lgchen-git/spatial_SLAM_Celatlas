#!/usr/bin/env python3
"""Create tissue barcode sets from automatic counts or TissUUmaps GeoJSON."""
import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def point_in_ring(x, y, ring):
    inside = False
    for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1 + 1e-12) + x1:
            inside = not inside
    return inside


def inside_geojson(x, y, data):
    features = data.get("features", [data]) if isinstance(data, dict) else []
    for feature in features:
        geom = feature.get("geometry", feature)
        if geom.get("type") == "Polygon":
            rings = geom["coordinates"]
            if point_in_ring(x, y, rings[0]) and not any(point_in_ring(x, y, hole) for hole in rings[1:]):
                return True
        elif geom.get("type") == "MultiPolygon":
            for polygon in geom["coordinates"]:
                if point_in_ring(x, y, polygon[0]) and not any(point_in_ring(x, y, hole) for hole in polygon[1:]):
                    return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-molecules", type=int, default=3)
    parser.add_argument("--celatlas-dir", help="Celatlas binSegment output directory")
    parser.add_argument("--roi", help="TissUUmaps GeoJSON; overrides automatic mask")
    parser.add_argument("--roi-coordinate-um", type=float, default=1.0, help="Multiply raw chip x/y before testing a TissUUmaps ROI")
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    molecule_counts = Counter()
    coordinates = {}
    for chunk in pd.read_csv(
        args.molecule_calls,
        sep="\t",
        compression="gzip",
        usecols=["barcode", "x", "y"],
        chunksize=args.chunk_size,
    ):
        grouped = chunk.groupby(["barcode", "x", "y"]).size()
        for (barcode, x, y), value in grouped.items():
            barcode = str(barcode)
            coordinates.setdefault(barcode, (x, y))
            molecule_counts[barcode] += int(value)
    spots = pd.DataFrame(
        [(barcode, *coordinates[barcode], molecule_counts[barcode]) for barcode in sorted(molecule_counts)],
        columns=["barcode", "x", "y", "molecules"],
    )
    if args.roi:
        roi = json.loads(Path(args.roi).read_text(encoding="utf-8"))
        keep = np.array([inside_geojson(row.x * args.roi_coordinate_um, row.y * args.roi_coordinate_um, roi) for row in spots.itertuples(index=False)])
        method = "tissuumaps_roi"
    elif args.celatlas_dir:
        candidates = sorted(Path(args.celatlas_dir).rglob("*tissue*.csv"))
        recovered = None
        recovered_barcodes = None
        for candidate in candidates:
            try:
                table = pd.read_csv(candidate)
                if not {"barcode", "in_tissue"}.issubset(table.columns) and table.shape[1] >= 6:
                    table = pd.read_csv(
                        candidate,
                        header=None,
                        names=[
                            "barcode",
                            "in_tissue",
                            "array_row",
                            "array_col",
                            "pxl_row_in_fullres",
                            "pxl_col_in_fullres",
                        ],
                    )
                if {"barcode", "in_tissue"}.issubset(table.columns):
                    recovered_barcodes = set(table.loc[table["in_tissue"].astype(int) == 1, "barcode"].astype(str))
                    break
                if {"x", "y"}.issubset(table.columns):
                    recovered = table[["x", "y"]].drop_duplicates()
                    break
            except (OSError, pd.errors.ParserError):
                continue
        if recovered_barcodes is not None:
            keep = spots["barcode"].astype(str).isin(recovered_barcodes).to_numpy()
            method = "celatlas_gene_expression_segmentation_barcode"
        elif recovered is not None:
            keep_xy = set(map(tuple, recovered.to_numpy()))
            keep = np.array([(row.x, row.y) in keep_xy for row in spots.itertuples(index=False)])
            method = "celatlas_gene_expression_segmentation"
        else:
            keep = spots["molecules"].to_numpy() >= args.min_molecules
            method = "celatlas_output_without_barcode_table_depth_fallback"
    else:
        # Celatlas auto-tissue output can replace this file through --roi. This conservative
        # depth filter is the deterministic fallback when an image/count segmentation is absent.
        keep = spots["molecules"].to_numpy() >= args.min_molecules
        method = "celatlas_auto_depth_fallback"
    tissue = spots.loc[keep].copy()
    tissue.to_csv(outdir / "01_tissue_barcodes.tsv", sep="\t", index=False)
    total_molecules = int(spots["molecules"].sum())
    tissue_total_molecules = int(tissue["molecules"].sum())
    (outdir / "tissue_mask_qc.json").write_text(json.dumps({
        "method": method,
        "n_spots": int(len(spots)),
        "n_tissue_spots": int(len(tissue)),
        "tissue_spot_fraction": float(len(tissue) / max(len(spots), 1)),
        "total_molecules": total_molecules,
        "tissue_total_molecules": tissue_total_molecules,
        "tissue_molecule_fraction": float(tissue_total_molecules / max(total_molecules, 1)),
        "min_molecules": args.min_molecules,
        "roi": args.roi,
        "roi_coordinate_um": args.roi_coordinate_um,
        "celatlas_dir": args.celatlas_dir,
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
