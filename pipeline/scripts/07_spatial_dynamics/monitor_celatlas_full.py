#!/usr/bin/env python3
"""Sample and summarize resources for each substep of a Celatlas full run."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import subprocess
import time
from collections import defaultdict
from pathlib import Path


STAGE_PATTERNS = (
    ("00_sample", (" rna sample ",)),
    ("01_barcode", (" rna barcode ",)),
    ("02_cutadapt", (" rna cutadapt ", "cutadapt ")),
    ("03_star", (" rna star ", "STAR ", "STAR --")),
    ("04_featureCounts", (" rna featureCounts ", "featureCounts ")),
    ("05_count", (" rna count ",)),
    ("06_binSegment", (" rna binSegment ",)),
    ("07_analysis", (" rna analysis ",)),
    ("09_report", ("celatlas_spatial_report", "spatial_report_generator")),
)


def stage_for(command: str) -> str | None:
    padded = f" {command} "
    for stage, patterns in STAGE_PATTERNS:
        if any(pattern in padded for pattern in patterns):
            return stage
    return None


def process_rows(workspace: str) -> list[tuple[str, float, int]]:
    output = subprocess.check_output(
        ["ps", "-eo", "pid=,%cpu=,rss=,args="], text=True, errors="replace"
    )
    rows: list[tuple[str, float, int]] = []
    own_pid = os.getpid()
    for line in output.splitlines():
        fields = line.strip().split(None, 3)
        if len(fields) != 4:
            continue
        pid_text, cpu_text, rss_text, command = fields
        if int(pid_text) == own_pid or workspace not in command:
            continue
        stage = stage_for(command)
        if stage is not None:
            rows.append((stage, float(cpu_text), int(rss_text)))
    return rows


def summarize(samples_path: Path, output_path: Path, interval: float) -> None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    if samples_path.exists():
        with samples_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                grouped[row["substage"]].append(row)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        fields = (
            "substage",
            "first_sample_utc",
            "last_sample_utc",
            "observed_wall_seconds",
            "peak_cpu_percent",
            "peak_rss_gb",
            "sample_count",
            "measurement_note",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for stage, _ in STAGE_PATTERNS:
            rows = grouped.get(stage, [])
            if not rows:
                continue
            writer.writerow(
                {
                    "substage": stage,
                    "first_sample_utc": rows[0]["timestamp_utc"],
                    "last_sample_utc": rows[-1]["timestamp_utc"],
                    "observed_wall_seconds": round((len(rows) - 1) * interval + interval, 1),
                    "peak_cpu_percent": round(max(float(r["cpu_percent"]) for r in rows), 1),
                    "peak_rss_gb": round(max(float(r["rss_mb"]) for r in rows) / 1024, 3),
                    "sample_count": len(rows),
                    "measurement_note": f"process sampling every {interval:g}s; short setup gaps may be absent",
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--exit-file", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--interval", type=float, default=5.0)
    args = parser.parse_args()
    if args.interval <= 0:
        parser.error("--interval must be positive")

    args.samples.parent.mkdir(parents=True, exist_ok=True)
    new_file = not args.samples.exists() or args.samples.stat().st_size == 0
    with args.samples.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        if new_file:
            writer.writerow(("timestamp_utc", "substage", "process_count", "cpu_percent", "rss_mb"))
            handle.flush()
        empty_samples = 0
        while not args.exit_file.exists():
            rows = process_rows(args.workspace)
            if rows:
                empty_samples = 0
                grouped: dict[str, list[tuple[float, int]]] = defaultdict(list)
                for stage, cpu, rss in rows:
                    grouped[stage].append((cpu, rss))
                timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
                for stage, values in grouped.items():
                    writer.writerow(
                        (
                            timestamp,
                            stage,
                            len(values),
                            round(sum(v[0] for v in values), 1),
                            round(sum(v[1] for v in values) / 1024, 1),
                        )
                    )
                handle.flush()
            else:
                empty_samples += 1
            # Exit if the target disappeared before it could write an exit file.
            if empty_samples >= max(12, int(300 / args.interval)):
                break
            time.sleep(args.interval)
    summarize(args.samples, args.summary, args.interval)


if __name__ == "__main__":
    main()
