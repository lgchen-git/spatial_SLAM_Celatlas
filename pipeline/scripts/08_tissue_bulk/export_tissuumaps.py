#!/usr/bin/env python3
"""Export Spatial SLAM bin markers and a portable TissUUmaps project."""
import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial-domains", required=True)
    parser.add_argument("--corrected-bin-gene", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--coordinate-um", type=float, required=True)
    args = parser.parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    domains = pd.read_csv(args.spatial_domains, sep="\t")
    counts_path = Path(args.corrected_bin_gene)
    counts = pd.read_pickle(counts_path) if counts_path.suffix == ".pkl" else pd.read_csv(counts_path, sep="\t", compression="gzip")
    total_column = "total_expected" if "total_expected" in counts.columns else "total"
    new_column = "new_expected" if "new_expected" in counts.columns else "expected_new"
    missing = {"bin_id", total_column, new_column} - set(counts.columns)
    if missing:
        raise ValueError(f"Corrected bin-gene table is missing required columns: {sorted(missing)}")
    summary = counts.groupby("bin_id", as_index=False).agg(total_expected=(total_column, "sum"), new_expected=(new_column, "sum"))
    summary["NTR"] = summary["new_expected"] / summary["total_expected"].clip(lower=1e-12)
    markers = domains.merge(summary, on="bin_id", how="left")
    markers["X"] = markers["bin_x"] * args.coordinate_um
    markers["Y"] = markers["bin_y"] * args.coordinate_um
    markers = markers.rename(columns={"cluster": "domain", "umi_count": "total_umi"})
    markers.to_csv(outdir / "spatial_slam_bins.csv", index=False)
    tmap = {
        "schemaVersion": "1.3",
        "filename": "spatial_slam.tmap",
        "mpp": 1.0,
        "markerFiles": [{
            "title": "Spatial SLAM bins",
            "name": "NTR and domains",
            "path": "spatial_slam_bins.csv",
            "autoLoad": True,
            "expectedHeader": {"X": "X", "Y": "Y", "gb_col": "domain", "cb_col": "NTR", "scale_col": "total_umi", "tooltip_fmt": "{domain}; NTR={NTR}; total={total_umi}"},
        }],
        "regionFiles": [],
    }
    (outdir / "spatial_slam.tmap").write_text(json.dumps(tmap, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
