#!/usr/bin/env python3
"""Write and validate stage input/output manifests and completion metadata."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

EXCLUDED = {"input_manifest.tsv", "stage_output_manifest.tsv", "stage_complete.json"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def input_rows(paths: list[Path]) -> list[dict[str, str | int]]:
    rows = []
    for path in paths:
        resolved = path.resolve()
        if not resolved.is_file():
            rows.append({"path": str(resolved), "bytes": -1, "mtime_ns": -1, "sha256": "missing"})
            continue
        stat = resolved.stat()
        rows.append({
            "path": str(resolved), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": digest(resolved) if stat.st_size <= 64 * 1024 * 1024 else "size_mtime_signature",
        })
    return rows


def signature(rows: list[dict[str, str | int]]) -> str:
    payload = json.dumps(rows, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def output_rows(stage: Path) -> list[dict[str, str | int]]:
    rows = []
    for path in sorted(
        item for item in stage.rglob("*")
        if item.is_file() and item.name not in EXCLUDED and "progress_runtime" not in item.parts
    ):
        stat = path.stat()
        rows.append({
            "relative_path": str(path.relative_to(stage)), "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": digest(path) if stat.st_size <= 64 * 1024 * 1024 else "size_mtime_signature",
        })
    return rows


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def check(stage: Path, current_signature: str) -> bool:
    marker = stage / "05_metadata/stage_complete.json"
    manifest = stage / "05_metadata/stage_output_manifest.tsv"
    if not marker.is_file() or not manifest.is_file():
        return False
    metadata = json.loads(marker.read_text(encoding="utf-8"))
    if metadata.get("input_signature") != current_signature or metadata.get("exit_code") != 0:
        return False
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            path = stage / row["relative_path"]
            if not path.is_file() or path.stat().st_size != int(row["bytes"]):
                return False
            if row["sha256"] == "size_mtime_signature":
                if path.stat().st_mtime_ns != int(row["mtime_ns"]):
                    return False
            elif digest(path) != row["sha256"]:
                return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--stage-id", required=True)
    parser.add_argument("--input", action="append", type=Path, default=[])
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--import-existing", action="store_true",
        help="Mark manifests as adopted from an existing successful result rather than a new execution",
    )
    args = parser.parse_args()
    stage = args.stage_dir.resolve()
    metadata = stage / "05_metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    rows = input_rows(args.input)
    current_signature = signature(rows)
    if args.check:
        return 0 if check(stage, current_signature) else 1
    write_tsv(metadata / "input_manifest.tsv", rows, ["path", "bytes", "mtime_ns", "sha256"])
    outputs = output_rows(stage)
    write_tsv(metadata / "stage_output_manifest.tsv", outputs, ["relative_path", "bytes", "mtime_ns", "sha256"])
    marker = {
        "schema_version": 1, "stage_id": args.stage_id, "exit_code": 0,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_signature": current_signature, "output_count": len(outputs),
        "completion_mode": "imported_existing_outputs" if args.import_existing else "executed",
    }
    (metadata / "stage_complete.json").write_text(json.dumps(marker, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
