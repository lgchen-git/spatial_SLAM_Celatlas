#!/usr/bin/env python3
"""Validate an existing Celatlas 1.8 result before read-only reuse."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def required_artifacts(source: Path, chip: str) -> dict[str, Path]:
    return {
        "metrics": source / ".metrics.json",
        "count_detail": source / f"05.count/{chip}_count_detail.txt",
        "bin50_matrix": source / f"06.segment/01.binsegment/square_bin/{chip}_bin50/filtered_feature_bc_matrix/matrix.mtx.gz",
        "bin50_features": source / f"06.segment/01.binsegment/square_bin/{chip}_bin50/filtered_feature_bc_matrix/features.tsv.gz",
        "bin50_barcodes": source / f"06.segment/01.binsegment/square_bin/{chip}_bin50/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "bin50_positions": source / f"06.segment/01.binsegment/square_bin/{chip}_bin50/spatial/tissue_positions_list.csv",
        "raw_positions": source / f"06.segment/01.binsegment/square_bin/{chip}_Raw/spatial/tissue_positions_list.csv",
        "scalefactors": source / f"06.segment/01.binsegment/square_bin/{chip}_Raw/spatial/scalefactors_json.json",
        "registration_image": source / f"06.segment/01.binsegment/images/{chip}_regist.tif",
        "hires_image": source / "06.segment/01.binsegment/images/tissue_hires_image.png",
        "lowres_image": source / "06.segment/01.binsegment/images/tissue_lowres_image.png",
        "bin50_h5ad": source / f"07.outs/binned_outputs/h5ad/{chip}_bin50.h5ad",
        "report": source / f"{chip}_spatial_analysis_report.html",
    }


def validate(source: Path, sample: str, chip: str, chemistry: str, species: str) -> dict:
    source = source.resolve()
    required = required_artifacts(source, chip)
    missing = [str(path) for path in required.values() if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise ValueError("Celatlas reuse source is incomplete:\n" + "\n".join(missing))
    data = json.loads((source / ".data.json").read_text(encoding="utf-8"))
    sample_summary = {
        item.get("name"): item.get("value")
        for item in data.get("sample_summary", {}).get("metric_list", [])
    }
    star_summary = {
        item.get("name"): item.get("value")
        for item in data.get("star_summary", {}).get("metric_list", [])
    }
    checks = {
        "chip": [sample_summary.get("Chip Number"), chip],
        "chemistry": [sample_summary.get("Chemistry"), chemistry],
        "species": [star_summary.get("Genome"), species],
    }
    bad = {name: values for name, values in checks.items() if values[0] != values[1]}
    if bad:
        raise ValueError("Celatlas reuse metadata mismatch: " + json.dumps(bad, ensure_ascii=False))
    return {
        "schema_version": 1,
        "ok": True,
        "reuse_mode": "validated_existing_result",
        "sample": sample,
        "chip": chip,
        "chemistry": chemistry,
        "species": species,
        "source_result_dir": str(source),
        "artifacts": {
            name: {"path": str(path), "bytes": path.stat().st_size}
            for name, path in required.items()
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chip", required=True)
    parser.add_argument("--chemistry", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        payload = validate(args.source, args.sample, args.chip, args.chemistry, args.species)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
