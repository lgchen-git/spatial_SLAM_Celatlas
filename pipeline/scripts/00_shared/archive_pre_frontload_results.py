#!/usr/bin/env python3
"""Archive pre-frontload results/runtime trees with a reversible move manifest."""
from __future__ import annotations

import argparse
import csv
import os
from datetime import datetime, timezone
from pathlib import Path

import yaml


FIELDS = ("source", "destination", "kind", "file_count", "bytes", "status")


def contained(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
        return True
    except ValueError:
        return False


def tree_stats(path: Path) -> tuple[int, int]:
    files = 0
    total_bytes = 0
    for current_root, _, names in os.walk(path):
        for name in names:
            item = Path(current_root) / name
            if item.is_symlink():
                continue
            files += 1
            total_bytes += item.stat().st_size
    return files, total_bytes


def write_manifest(path: Path, rows: list[dict]) -> None:
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(partial, path)


def planned_rows(config: dict, samples_path: Path, stamp: str) -> tuple[Path, list[dict]]:
    results_root = Path(config["results_root"]).resolve()
    tmp_root = Path(config["tmp_root"]).resolve()
    project_id = str(config["project_id"])
    with samples_path.open(encoding="utf-8", newline="") as handle:
        sample_ids = [
            row["sample_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["project_id"] == project_id
        ]
    if not sample_ids:
        raise ValueError(f"no samples found for project {project_id}")

    archive_root = tmp_root / project_id / "results_legacy" / f"pre_celatlas_frontload_{stamp}"
    candidates: list[tuple[Path, Path, str]] = []
    for sample_id in sample_ids:
        candidates.append((results_root / sample_id, archive_root / "results" / sample_id, "result"))
        candidates.append((tmp_root / project_id / sample_id, archive_root / "runtime" / sample_id, "runtime"))
    candidates.append((results_root / "_project", archive_root / "results" / "_project", "project_qc"))

    rows = []
    for source, destination, kind in candidates:
        if not source.exists():
            continue
        if source.is_symlink() or not source.is_dir():
            raise ValueError(f"archive source must be a real directory: {source}")
        if not contained(source, results_root) and not contained(source, tmp_root):
            raise ValueError(f"source outside configured roots: {source}")
        if not contained(destination, tmp_root):
            raise ValueError(f"destination outside configured tmp root: {destination}")
        file_count, total_bytes = tree_stats(source)
        rows.append(
            {
                "source": str(source),
                "destination": str(destination),
                "kind": kind,
                "file_count": file_count,
                "bytes": total_bytes,
                "status": "planned",
            }
        )
    return archive_root, rows


def apply_archive(archive_root: Path, rows: list[dict]) -> Path:
    if archive_root.exists():
        raise FileExistsError(f"archive destination already exists: {archive_root}")
    for row in rows:
        destination = Path(row["destination"])
        if destination.exists():
            raise FileExistsError(destination)
        if destination.parent.exists() and destination.parent.is_symlink():
            raise ValueError(f"refusing symlink destination parent: {destination.parent}")

    archive_root.mkdir(parents=True)
    manifest = archive_root / "migration_manifest.tsv"
    write_manifest(manifest, rows)
    for row in rows:
        source = Path(row["source"])
        destination = Path(row["destination"])
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.rename(source, destination)
        row["status"] = "archived"
        write_manifest(manifest, rows)
    return manifest


def rollback(manifest: Path) -> None:
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    active = [row for row in rows if row["status"] == "archived"]
    for row in active:
        source = Path(row["source"])
        destination = Path(row["destination"])
        if source.exists():
            raise FileExistsError(f"cannot restore over existing source: {source}")
        if not destination.is_dir() or destination.is_symlink():
            raise FileNotFoundError(f"archived directory unavailable: {destination}")
    for row in reversed(active):
        source = Path(row["source"])
        destination = Path(row["destination"])
        source.parent.mkdir(parents=True, exist_ok=True)
        os.rename(destination, source)
        row["status"] = "restored"
        write_manifest(manifest, rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/analysis_layout.yaml"))
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument("--stamp", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    parser.add_argument("--apply", action="store_true", help="Apply the planned directory moves")
    parser.add_argument("--rollback-manifest", type=Path, help="Restore an applied archive manifest")
    args = parser.parse_args()

    if args.rollback_manifest:
        if args.apply:
            parser.error("--apply and --rollback-manifest are mutually exclusive")
        rollback(args.rollback_manifest)
        print(args.rollback_manifest)
        return 0

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    archive_root, rows = planned_rows(config, args.samples, args.stamp)
    print(f"archive_root\t{archive_root}")
    print(f"directories\t{len(rows)}")
    print(f"files\t{sum(int(row['file_count']) for row in rows)}")
    print(f"bytes\t{sum(int(row['bytes']) for row in rows)}")
    for row in rows:
        print(f"{row['kind']}\t{row['source']}\t{row['destination']}")
    if not args.apply:
        print("dry_run\ttrue")
        return 0
    manifest = apply_archive(archive_root, rows)
    print(f"manifest\t{manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
