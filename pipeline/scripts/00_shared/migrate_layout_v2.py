#!/usr/bin/env python3
"""Plan, apply, or roll back the Celatlas Spatial SLAM v2 layout migration."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

MANIFEST_FIELDS = ("source", "destination", "bytes", "sha256", "reason", "status")


@dataclass(frozen=True)
class Move:
    source: Path
    destination: Path
    reason: str


def scalar(value: str):
    value = value.strip()
    try:
        return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value


def load_config(path: Path) -> dict:
    try:
        import yaml  # type: ignore
    except ModuleNotFoundError:
        data = {}
        for raw in path.read_text(encoding="utf-8").splitlines():
            line = raw.split("#", 1)[0].rstrip()
            if line and not line.startswith(" "):
                key, value = line.split(":", 1)
                if value.strip():
                    data[key.strip()] = scalar(value)
        return data
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def same_file_content(left: Path, right: Path) -> bool:
    return left.stat().st_size == right.stat().st_size and sha256(left) == sha256(right)


def add_files(moves: list[Move], source_dir: Path, destination_dir: Path, reason: str) -> None:
    if source_dir.is_dir():
        for source in sorted(path for path in source_dir.rglob("*") if path.is_file() or path.is_symlink()):
            moves.append(Move(source, destination_dir / source.relative_to(source_dir), reason))


def legacy_path(legacy: Path, source: Path, sample: Path) -> Path:
    return legacy / source.relative_to(sample)


def build_moves(sample: Path, tmp_root: Path, stamp: str) -> list[Move]:
    """Central, reviewable mapping from the legacy sample layout to v2."""
    legacy = tmp_root / "results_legacy" / sample.name / stamp
    runtime_tmp = tmp_root / sample.parent.name / sample.name
    moves: list[Move] = []

    def move(relative: str, destination: Path, reason: str) -> None:
        source = sample / relative
        if source.exists() or source.is_symlink():
            moves.append(Move(source, destination, reason))

    # Stage 02: retained tables/QC plus a non-deliverable large FASTQ.
    for source in sorted((sample / "02_barcode").glob("*.read_barcode_umi.tsv.gz")):
        moves.append(Move(source, sample / "02_barcode/01_tables" / source.name, "stage02_table"))
    for source in sorted((sample / "02_barcode").glob("*.barcode_qc.json")):
        moves.append(Move(source, sample / "02_barcode/05_metadata" / source.name, "stage02_qc"))
    for source in sorted((sample / "02_barcode").glob("*.filtered_r2.fastq.gz")):
        moves.append(Move(source, runtime_tmp / "02_barcode" / source.name, "large_intermediate"))

    # Stage 03: Celatlas/STAR output, logs, tables, and large cutadapt FASTQ.
    stage03 = sample / "03_trim_align"
    for pattern in ("*.bam", "*_report.html"):
        for source in sorted(stage03.glob(pattern)):
            moves.append(Move(source, stage03 / "04_tool_output/celatlas" / source.name, "vendor_output"))
    for pattern in ("*_Log*.out", "*.log"):
        for source in sorted(stage03.glob(pattern)):
            moves.append(Move(source, stage03 / "03_logs" / source.name, "stage03_log"))
    for pattern in ("*_SJ.out.tab", "stat.txt"):
        for source in sorted(stage03.glob(pattern)):
            moves.append(Move(source, stage03 / "01_tables" / source.name, "stage03_table"))
    for source in sorted(stage03.glob(".*.json")):
        moves.append(Move(source, stage03 / "05_metadata" / source.name.lstrip("."), "stage03_metadata"))
    cutadapt = stage03 / "cutadapt"
    for source in sorted(cutadapt.glob("*")) if cutadapt.is_dir() else []:
        if not source.is_file():
            continue
        destination = (runtime_tmp / "03_trim_align/cutadapt" / source.name) if source.name.endswith((".fq.gz", ".fastq.gz")) else (stage03 / "03_logs" / f"cutadapt_{source.name}")
        moves.append(Move(source, destination, "cutadapt_intermediate" if destination.is_relative_to(tmp_root) else "cutadapt_log"))

    # Stage 04: molecule tables, QC, and retained annotation tool output.
    stage04 = sample / "04_molecule_call"
    for pattern in ("*.summary.tsv", "*.molecule_calls.tsv.gz"):
        for source in sorted(stage04.glob(pattern)):
            moves.append(Move(source, stage04 / "01_tables" / source.name, "stage04_table"))
    for source in sorted(stage04.glob("*.qc.json")):
        moves.append(Move(source, stage04 / "05_metadata" / source.name, "stage04_qc"))
    add_files(moves, stage04 / "annotation", stage04 / "04_tool_output/annotation", "annotation_output")

    # Stage 05: one canonical matrix schema under retained tool output.
    stage05 = sample / "05_matrix"
    canonical = {"old_umi.mtx", "new_umi.mtx", "features.tsv", "barcodes.tsv"}
    for source in sorted(path for path in stage05.glob("*") if path.is_file()):
        stripped = source.name.removeprefix(sample.name + ".")
        if stripped in canonical and source.name != stripped:
            destination = stage05 / "04_tool_output/matrix" / stripped
            existing = stage05 / stripped
            if existing.exists():
                destination = legacy_path(legacy, source, sample)
            moves.append(Move(source, destination, "canonical_matrix" if destination.is_relative_to(stage05) else "duplicate_matrix"))
    for name in canonical:
        move(f"05_matrix/{name}", stage05 / "04_tool_output/matrix" / name, "canonical_matrix")

    # Stage 06: flatten the nested stage; retire two redundant figures and renumber the retained combined view.
    showcase = sample / "06_spatial_dynamics/bin50_showcase"
    for subdir in ("01_tables", "04_tool_output", "05_metadata"):
        source_dir = showcase / subdir
        if source_dir.is_dir():
            for source in sorted(path for path in source_dir.rglob("*") if path.is_file() or path.is_symlink()):
                destination = sample / "06_spatial_dynamics" / subdir / source.relative_to(source_dir)
                if destination.exists() or destination.is_symlink():
                    # Existing stage-level metadata is authoritative; retain the nested historical copy.
                    destination = legacy_path(legacy, source, sample)
                    reason = "duplicate_bin50_metadata"
                else:
                    reason = "flatten_bin50"
                moves.append(Move(source, destination, reason))
    if (showcase / "03_logs").is_dir():
        for source in sorted(path for path in (showcase / "03_logs").rglob("*") if path.is_file()):
            if source.name != "tissuumaps_export.time.txt":
                moves.append(Move(source, sample / "06_spatial_dynamics/03_logs" / source.relative_to(showcase / "03_logs"), "flatten_bin50"))
    figures = showcase / "02_figures"
    if figures.is_dir():
        for source in sorted(path for path in figures.iterdir() if path.is_file()):
            if source.name.startswith(("06_representative_dynamic_genes_spatial_ntr", "07_spatial_slam_seq_capability_overview")):
                destination = legacy_path(legacy, source, sample)
                reason = "redundant_spatial_figure"
            elif source.name.startswith("08_representative_dynamic_genes_old_new_ntr"):
                destination = sample / "06_spatial_dynamics/02_figures" / source.name.replace("08_", "06_", 1)
                reason = "renumber_spatial_figure"
            else:
                destination = sample / "06_spatial_dynamics/02_figures" / source.name
                reason = "flatten_bin50"
            moves.append(Move(source, destination, reason))
    move("06_spatial_dynamics/bin50_showcase/03_logs/tissuumaps_export.time.txt", sample / "07_tissue_bulk/03_logs/tissuumaps_export.time.txt", "stage07_log")

    # Stage 07: one tool-output location and one ordered, non-duplicated figure series.
    add_files(moves, sample / "07_tissue_bulk/tissuumaps", sample / "07_tissue_bulk/04_tool_output/tissuumaps", "tissuumaps_tool_output")
    stage07_tables = sample / "07_tissue_bulk/01_tables"
    table_renames = {
        "tissue_barcodes.tsv": "01_tissue_barcodes.tsv",
        "tissue_bins.tsv": "01_tissue_bins.tsv",
        "utr.bed": "02_utr_regions.bed",
        "04_complete_12_channel_mismatch_spectrum.tsv": "05_complete_12_channel_mismatch_spectrum.tsv",
    }
    for old, new in table_renames.items():
        move(f"07_tissue_bulk/01_tables/{old}", stage07_tables / new, "renumber_bulk_table")
    stage07_fig = sample / "07_tissue_bulk/02_figures"
    renames = {
        "01_expression_derived_tissue_mask": "01_tissue_mask_expression_derived",
        "01_slamdunk_conversion_qc": "02_slamdunk_conversion_qc",
        "02_grand_slam_gene_ntr": "03_grand_slam_gene_ntr",
        "03_grand_slam_model_qc": "04_grand_slam_model_qc",
        "04_tc_mismatch_enrichment_over_background": "05_tc_mismatch_enrichment_over_background",
    }
    for source in sorted(stage07_fig.glob("*")) if stage07_fig.is_dir() else []:
        if not source.is_file():
            continue
        if source.name.startswith("03_grand_slam_mismatch_barplot"):
            moves.append(Move(source, legacy_path(legacy, source, sample), "redundant_mismatch_figure"))
            continue
        for old, new in renames.items():
            if source.name.startswith(old):
                moves.append(Move(source, stage07_fig / source.name.replace(old, new, 1), "renumber_bulk_figure"))
                break
    move("07_tissue_bulk/01_tables/grand_slam_mismatch_barplot_significance.tsv", legacy / "07_tissue_bulk/01_tables/grand_slam_mismatch_barplot_significance.tsv", "obsolete_significance_table")

    # Stage 08 old report-only tables/figures and sample-root recovery logs are not deliverables.
    for relative in ("08_report/01_tables/02_qc_comparison_spatialrna_bin50.tsv", "08_report/01_tables/03_total_new_spatial_difference_genes.tsv"):
        move(relative, legacy / relative, "obsolete_report_table")
    add_files(moves, sample / "08_report/02_figures", legacy / "08_report/02_figures", "obsolete_report_figure")
    for source in sorted(path for path in sample.iterdir() if path.is_file()):
        if source.name.endswith((".log", ".exit")) or ".recovery." in source.name:
            moves.append(Move(source, legacy_path(legacy, source, sample), "sample_root_log"))
    add_files(moves, sample / ".runtime", legacy / ".runtime", "sample_root_runtime")
    return moves


def validate_move(move: Move, sample: Path, tmp_root: Path) -> None:
    source = move.source.absolute()
    destination = move.destination.absolute()
    source.relative_to(sample.absolute())
    if not (destination.is_relative_to(sample.absolute()) or destination.is_relative_to(tmp_root.absolute())):
        raise ValueError(f"destination outside configured result/tmp roots: {destination}")
    source_parts = source.relative_to(sample.absolute()).parts
    if "04_tool_output" in source_parts and "bin50_showcase" not in source_parts:
        raise ValueError(f"refusing to migrate inside vendor output: {source}")
    if source.is_symlink():
        raise ValueError(f"refusing to migrate symlink: {source}")
    if source.is_dir():
        raise ValueError(f"only files may be migrated: {source}")


def nearest_existing_parent(path: Path) -> Path:
    candidate = path
    while not candidate.exists():
        if candidate == candidate.parent:
            raise FileNotFoundError(path)
        candidate = candidate.parent
    return candidate


def execute(moves: list[Move], sample: Path, tmp_root: Path, apply: bool) -> list[dict[str, str]]:
    rows = []
    seen_destinations: set[Path] = set()
    # Complete preflight before the first mutation, so conflicts cannot produce a partial migration.
    for move in moves:
        validate_move(move, sample, tmp_root)
        if move.destination in seen_destinations:
            raise RuntimeError(f"multiple sources map to {move.destination}")
        seen_destinations.add(move.destination)
        size = move.source.stat().st_size
        digest = sha256(move.source) if apply else ""
        status = "planned"
        if move.destination.exists() or move.destination.is_symlink():
            if not move.destination.is_file() or not same_file_content(move.source, move.destination):
                raise FileExistsError(f"different destination already exists: {move.destination}")
            # Preserve both copies reversibly: the central map must explicitly route duplicates to legacy.
            raise FileExistsError(f"identical destination already exists; map duplicate to legacy: {move.destination}")
        if apply and move.source.stat().st_dev != nearest_existing_parent(move.destination.parent).stat().st_dev:
            raise OSError(f"migration requires same-filesystem rename: {move.source} -> {move.destination}")
        rows.append({"source": str(move.source), "destination": str(move.destination), "bytes": str(size), "sha256": digest, "reason": move.reason, "status": status})
    if apply:
        for move, row in zip(moves, rows):
            move.destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(move.source, move.destination)
            if move.destination.stat().st_size != int(row["bytes"]) or sha256(move.destination) != row["sha256"]:
                raise OSError(f"post-move verification failed: {move.destination}")
            row["status"] = "moved"
    return rows


def write_manifest(path: Path, rows: list[dict[str, str]], exclusive: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x" if exclusive else "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def remove_empty_legacy_containers(sample: Path) -> None:
    """Remove only now-empty wrapper directories superseded by the standard layout."""
    roots = (
        sample / "03_trim_align/cutadapt",
        sample / "04_molecule_call/annotation",
        sample / "06_spatial_dynamics/bin50_showcase",
        sample / "07_tissue_bulk/tissuumaps",
        sample / ".runtime",
    )
    for root in roots:
        if not root.is_dir():
            continue
        for directory in sorted((path for path in root.rglob("*") if path.is_dir()), key=lambda path: len(path.parts), reverse=True):
            try:
                directory.rmdir()
            except OSError:
                pass
        try:
            root.rmdir()
        except OSError:
            pass


def rollback(path: Path) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in reversed(rows):
        source, destination = Path(row["source"]), Path(row["destination"])
        if row.get("status") != "moved":
            continue
        if not destination.is_file():
            raise FileNotFoundError(f"rollback source missing: {destination}")
        if source.exists() or source.is_symlink():
            raise FileExistsError(f"rollback destination exists: {source}")
        if destination.stat().st_size != int(row["bytes"]) or sha256(destination) != row["sha256"]:
            raise OSError(f"rollback checksum mismatch: {destination}")
        source.parent.mkdir(parents=True, exist_ok=True)
        os.replace(destination, source)
        if source.stat().st_size != int(row["bytes"]) or sha256(source) != row["sha256"]:
            raise OSError(f"rollback verification failed: {source}")
    print(json.dumps({"rolled_back": str(path), "moves": sum(row.get("status") == "moved" for row in rows)}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/analysis_layout.yaml"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--sample", help="sample ID (required unless rolling back)")
    parser.add_argument("--timestamp", help="stable migration stamp; defaults to current UTC")
    parser.add_argument("--apply", action="store_true", help="perform moves; default is a JSON dry-run")
    parser.add_argument("--rollback-manifest", type=Path)
    args = parser.parse_args()
    if args.rollback_manifest:
        if args.apply or args.sample:
            parser.error("--rollback-manifest cannot be combined with --apply or --sample")
        return rollback(args.rollback_manifest.resolve())
    if not args.sample:
        parser.error("--sample is required")
    repo = args.repo_root.resolve()
    config_path = args.config if args.config.is_absolute() else repo / args.config
    config = load_config(config_path)
    results_root = (args.results_root or Path(config["results_root"])).resolve()
    tmp_root = Path(config["tmp_root"]).resolve()
    sample = results_root / args.sample
    if not sample.is_dir():
        raise FileNotFoundError(sample)
    stamp = args.timestamp or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    moves = build_moves(sample, tmp_root, stamp)
    manifest = sample / "08_report/05_metadata/layout_migration_v2.tsv"
    if args.apply and (manifest.exists() or manifest.is_symlink()):
        raise FileExistsError(f"refusing to overwrite existing migration manifest: {manifest}")
    rows = execute(moves, sample, tmp_root, args.apply)
    if args.apply:
        write_manifest(manifest, rows)
        remove_empty_legacy_containers(sample)
        print(json.dumps({"mode": "apply", "moves": len(rows), "manifest": str(manifest)}, indent=2))
    else:
        print(json.dumps({"mode": "dry-run", "moves": rows}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, OSError, ValueError, RuntimeError) as exc:
        print(f"migration error: {exc}", file=sys.stderr)
        raise SystemExit(2)
