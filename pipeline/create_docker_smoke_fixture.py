#!/usr/bin/env python3
"""Create an isolated synchronized real-FASTQ Docker smoke project."""
from __future__ import annotations

import argparse
import csv
import gzip
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def normalized_read_id(header: str) -> str:
    return header.split()[0].removesuffix("/1").removesuffix("/2")


def subset_pair(source_r1: Path, source_r2: Path, target_r1: Path, target_r2: Path, pairs: int) -> None:
    target_r1.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source_r1, "rt", encoding="utf-8", errors="strict") as r1, \
         gzip.open(source_r2, "rt", encoding="utf-8", errors="strict") as r2, \
         gzip.open(target_r1, "wt", encoding="utf-8") as out1, \
         gzip.open(target_r2, "wt", encoding="utf-8") as out2:
        for pair_index in range(pairs):
            record1 = [r1.readline() for _ in range(4)]
            record2 = [r2.readline() for _ in range(4)]
            if not all(record1) or not all(record2):
                raise RuntimeError(f"source contains fewer than {pairs} read pairs: {source_r1}, {source_r2}")
            if normalized_read_id(record1[0]) != normalized_read_id(record2[0]):
                raise RuntimeError(f"R1/R2 read ID mismatch at pair {pair_index + 1}")
            out1.writelines(record1)
            out2.writelines(record2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--pairs-per-library", type=int, default=500000)
    args = parser.parse_args()
    target = args.target.resolve()
    if target.exists() and any(target.iterdir()):
        raise SystemExit(f"target must be absent or empty: {target}")
    for name in ("config", "inputs/fastq", "inputs/celatlas_reuse", "log", "outputs", "tmp"):
        (target / name).mkdir(parents=True, exist_ok=True)
    os.symlink(ROOT / "pipeline", target / "pipeline", target_is_directory=True)
    os.symlink(ROOT / "inputs/ST_chip", target / "inputs/ST_chip", target_is_directory=True)
    os.symlink(ROOT / "inputs/celatlas_reuse/example_sample", target / "inputs/celatlas_reuse/example_sample", target_is_directory=True)
    os.symlink(ROOT / "reference", target / "reference", target_is_directory=True)

    libraries = []
    with (ROOT / "config/libraries_example_run.tsv").open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in source_rows:
        output_row = {"sample_id": row["sample_id"], "library_id": row["library_id"]}
        source_r1 = (ROOT / row["fastq_r1"]).resolve()
        source_r2 = (ROOT / row["fastq_r2"]).resolve()
        destination_r1 = target / "inputs/fastq" / source_r1.name
        destination_r2 = target / "inputs/fastq" / source_r2.name
        subset_pair(source_r1, source_r2, destination_r1, destination_r2, args.pairs_per_library)
        output_row["fastq_r1"] = f"inputs/fastq/{source_r1.name}"
        output_row["fastq_r2"] = f"inputs/fastq/{source_r2.name}"
        libraries.append(output_row)
    with (target / "config/libraries_example_run.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sample_id", "library_id", "fastq_r1", "fastq_r2"), delimiter="\t")
        writer.writeheader(); writer.writerows(libraries)
    (target / "config/samples_example_run.tsv").write_text(
        (ROOT / "config/samples_example_run.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    config_text = (ROOT / "config/run_example_run.yaml").read_text(encoding="utf-8")
    config_text = config_text.replace("threads: 32\n", "threads: 8\n", 1)
    config_text = config_text.replace("memory_mb: 98304\n", "memory_mb: 32768\n", 1)
    config_text += "smoke: true\n"
    (target / "config/run_example_run.yaml").write_text(config_text, encoding="utf-8")
    print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
