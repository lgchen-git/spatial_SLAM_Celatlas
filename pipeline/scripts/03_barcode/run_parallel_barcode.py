#!/usr/bin/env python3
"""Run barcode extraction per FASTQ lane, then atomically merge lane outputs."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def split_paths(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def run_lane(index: int, r1: str, r2: str, args, extractor: Path) -> dict:
    lane = f"lane{index:03d}"
    lane_sample = f"{args.sample}.{lane}"
    lane_root = args.tmp_dir / lane
    filtered = lane_root / f"{lane_sample}.filtered_r2.fastq.gz"
    table = lane_root / "01_tables" / f"{lane_sample}.read_barcode_umi.tsv.gz"
    qc_path = lane_root / "05_metadata" / f"{lane_sample}.barcode_qc.json"
    log = args.log_dir / f"barcode_{lane}.log"
    lane_root.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable, str(extractor), "--r1", r1, "--r2", r2,
        "--whitelist", str(args.whitelist), "--sample", lane_sample,
        "--outdir", str(lane_root), "--filtered-r2-output", str(filtered),
        "--pattern", args.pattern, "--min-polyt", str(args.min_polyt),
        "--gzip-threads", str(args.gzip_threads),
    ]
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    return {"lane": lane, "filtered": filtered, "table": table, "qc": qc, "log": log}


def merge_gzip_members(inputs: list[Path], output: Path) -> None:
    partial = output.with_name(f".{output.name}.partial.{os.getpid()}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with partial.open("wb") as target:
        for source in inputs:
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, target, length=8 * 1024 * 1024)
    os.replace(partial, output)


def merge_tables(inputs: list[Path], output: Path) -> None:
    partial = output.with_name(f".{output.name}.partial.{os.getpid()}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(partial, "wt") as target:
        expected_header = None
        for source in inputs:
            with gzip.open(source, "rt") as handle:
                header = handle.readline()
                if expected_header is None:
                    expected_header = header
                    target.write(header)
                elif header != expected_header:
                    raise ValueError(f"lane table header mismatch: {source}")
                shutil.copyfileobj(handle, target, length=8 * 1024 * 1024)
    os.replace(partial, output)


def aggregate_qc(parts: list[dict], args, table: Path, filtered: Path) -> dict:
    qcs = [part["qc"] for part in parts]
    total = sum(int(qc["total_reads"]) for qc in qcs)
    kept = sum(int(qc["kept_reads"]) for qc in qcs)
    result = {
        "sample": args.sample, "pattern": args.pattern, "input_fastq_pairs": len(qcs),
        "total_reads": total, "kept_reads": kept,
        "kept_fraction": kept / total if total else 0.0,
    }
    for key in ("failed_pattern", "barcode_not_in_whitelist", "failed_poly_t"):
        result[key] = sum(int(qc[key]) for qc in qcs)
    for key in ("mean_r1_quality", "mean_r2_quality"):
        result[key] = round(sum(float(qc[key]) * int(qc["kept_reads"]) for qc in qcs) / kept, 3) if kept else 0.0
    for key in ("median_umi_length", "median_polyt_length"):
        values = {int(qc[key]) for qc in qcs}
        if len(values) != 1:
            raise ValueError(f"lane QC disagreement for {key}: {sorted(values)}")
        result[key] = values.pop()
    result.update({
        "barcode_table": str(table), "filtered_r2": str(filtered),
        "execution": "parallel_by_fastq_pair", "workers": min(args.workers, len(parts)),
        "lane_logs": [str(part["log"]) for part in parts],
    })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--whitelist", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--filtered-r2-output", type=Path, required=True)
    parser.add_argument("--tmp-dir", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--pattern", default="C4L15C4L15C4U10T18")
    parser.add_argument("--min-polyt", type=int, default=12)
    parser.add_argument("--gzip-threads", type=int, default=2)
    args = parser.parse_args()
    r1, r2 = split_paths(args.r1), split_paths(args.r2)
    if len(r1) != len(r2) or not r1:
        raise ValueError(f"R1/R2 FASTQ count mismatch: {len(r1)} != {len(r2)}")
    extractor = Path(__file__).with_name("extract_spatial_barcodes.py")
    with ThreadPoolExecutor(max_workers=min(args.workers, len(r1))) as pool:
        futures = [pool.submit(run_lane, index, left, right, args, extractor) for index, (left, right) in enumerate(zip(r1, r2), 1)]
        parts = [future.result() for future in futures]
    table = args.outdir / "01_tables" / f"{args.sample}.read_barcode_umi.tsv.gz"
    qc_path = args.outdir / "05_metadata" / f"{args.sample}.barcode_qc.json"
    merge_gzip_members([part["filtered"] for part in parts], args.filtered_r2_output)
    merge_tables([part["table"] for part in parts], table)
    qc = aggregate_qc(parts, args, table, args.filtered_r2_output)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    partial_qc = qc_path.with_name(f".{qc_path.name}.partial.{os.getpid()}")
    partial_qc.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    os.replace(partial_qc, qc_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
