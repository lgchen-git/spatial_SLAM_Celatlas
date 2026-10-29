#!/usr/bin/env python3
"""Run a full-read A/B spatial-chip assignment cross-check before Celatlas."""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SCHEMA_VERSION = 1
METHOD = "full_r1_exact_primary_hamming1_directional_cross_whitelist"


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial.{os.getpid()}")
    with partial.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(partial, path)


def compile_scanner(source: Path, binary: Path) -> None:
    rustc = shutil.which("rustc")
    if rustc is None:
        raise RuntimeError("rustc is required for full-read chip-assignment QC")
    linker = shutil.which("cc") or shutil.which("x86_64-conda-linux-gnu-cc")
    if linker is None:
        raise RuntimeError("a C linker is required for the Rust chip-assignment scanner")
    if binary.is_file() and binary.stat().st_mtime_ns >= source.stat().st_mtime_ns:
        return
    binary.parent.mkdir(parents=True, exist_ok=True)
    partial = binary.with_name(f".{binary.name}.partial.{os.getpid()}")
    subprocess.run([rustc, str(source), "-C", f"linker={linker}", "-O", "-o", str(partial)], check=True)
    os.replace(partial, binary)


def run_scan(
    binary: Path,
    sample_id: str,
    candidates: list[dict[str, str]],
    r1_paths: list[str],
    threads: int,
) -> list[dict[str, str]]:
    command = [str(binary), sample_id, f"--targets={len(candidates)}"]
    for candidate in candidates:
        whitelist = Path(candidate["mask_dir"]) / f"{candidate['chip_id']}_FilterBarcodes.csv"
        if not whitelist.is_file():
            raise FileNotFoundError(whitelist)
        command.append(f"{candidate['chip_id']}={whitelist}")
    command.extend(r1_paths)
    environment = os.environ.copy()
    environment["CHIP_QC_PIGZ_THREADS"] = str(threads)
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        capture_output=True,
        env=environment,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"chip assignment scanner failed with {completed.returncode}: {completed.stderr.strip()}"
        )
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    if len(rows) != len(candidates):
        raise RuntimeError(f"scanner returned {len(rows)} rows for {sample_id}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--libraries", type=Path, required=True)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tmp-dir", type=Path, required=True)
    parser.add_argument("--minimum-fold", type=float, default=2.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.minimum_fold <= 1:
        raise ValueError("--minimum-fold must be greater than 1")
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    table_path = args.output_dir / "01_tables" / "chip_assignment_crosscheck.tsv"
    json_path = args.output_dir / "05_metadata" / "chip_assignment_qc.json"
    if table_path.is_file() and json_path.is_file() and not args.force:
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        if payload.get("method") == METHOD:
            print(json.dumps(payload, indent=2))
            return 0 if payload.get("passed") else 3

    samples = [row for row in read_table(args.samples) if row["project_id"] == args.project_id]
    if len(samples) not in (1, 2):
        raise ValueError(
            f"chip cross-check supports one or two project samples; found {len(samples)}"
        )
    if len(samples) == 2 and len({row["chip_id"] for row in samples}) != 2:
        raise ValueError("chip IDs must be unique within the two-sample project")
    libraries = read_table(args.libraries)
    r1_by_sample = {
        row["sample_id"]: [
            library["fastq_r1"]
            for library in libraries
            if library["sample_id"] == row["sample_id"]
        ]
        for row in samples
    }
    for sample_id, paths in r1_by_sample.items():
        if not paths:
            raise ValueError(f"no R1 FASTQ configured for {sample_id}")
        for path in paths:
            if not Path(path).is_file():
                raise FileNotFoundError(path)

    if table_path.is_file() and not args.force:
        raw_rows = read_table(table_path)
    else:
        source = Path(__file__).with_name("chip_assignment_scan.rs")
        binary = args.tmp_dir / "bin" / "chip_assignment_scan"
        compile_scanner(source, binary)
        raw_rows: list[dict[str, str]] = []
        for sample in samples:
            raw_rows.extend(
                run_scan(
                    binary,
                    sample["sample_id"],
                    samples,
                    r1_by_sample[sample["sample_id"]],
                    args.threads,
                )
            )

    declared_by_sample = {row["sample_id"]: row["chip_id"] for row in samples}
    result_rows: list[dict] = []
    exact_rates: dict[tuple[str, str], float] = {}
    corrected_rates: dict[tuple[str, str], float] = {}
    for row in raw_rows:
        total = int(row["total_reads"])
        corrected_poly_t = int(row["corrected_poly_t_reads"])
        rate = corrected_poly_t / total if total else 0.0
        exact_poly_t = int(row["exact_poly_t_reads"])
        exact_rate = exact_poly_t / total if total else 0.0
        exact_rates[(row["sample_id"], row["tested_chip_id"])] = exact_rate
        corrected_rates[(row["sample_id"], row["tested_chip_id"])] = rate
        result_rows.append(
            {
                **row,
                "one_mismatch_reads": int(row["corrected_reads"]) - int(row["exact_reads"]),
                "one_mismatch_poly_t_reads": corrected_poly_t - exact_poly_t,
                "exact_poly_t_rate": f"{exact_rate:.10f}",
                "corrected_poly_t_rate": f"{rate:.10f}",
                "declared_chip": str(row["tested_chip_id"] == declared_by_sample[row["sample_id"]]).lower(),
            }
        )

    columns = [
        "sample_id",
        "tested_chip_id",
        "whitelist_barcodes",
        "total_reads",
        "poly_t_reads",
        "exact_reads",
        "corrected_reads",
        "exact_poly_t_reads",
        "corrected_poly_t_reads",
        "one_mismatch_reads",
        "one_mismatch_poly_t_reads",
        "exact_poly_t_rate",
        "corrected_poly_t_rate",
        "declared_chip",
    ]

    if len(samples) == 1:
        sample = samples[0]
        sample_id = sample["sample_id"]
        own = next(row for row in raw_rows if row["sample_id"] == sample_id and row["tested_chip_id"] == sample["chip_id"])
        total = int(own["total_reads"])
        exact = int(own["exact_reads"])
        exact_poly_t = int(own["exact_poly_t_reads"])
        single_sample_passed = total > 0 and exact > 0 and exact_poly_t > 0
        sample_decisions = [
            {
                "sample_id": sample_id,
                "declared_chip_id": sample["chip_id"],
                "alternate_chip_id": "NA",
                "declared_exact_poly_t_rate": exact_rates[(sample_id, sample["chip_id"])],
                "alternate_exact_poly_t_rate": None,
                "declared_over_alternate_exact_fold": None,
                "declared_corrected_poly_t_rate": corrected_rates[(sample_id, sample["chip_id"])],
                "alternate_corrected_poly_t_rate": None,
                "declared_over_alternate_corrected_fold": None,
                "decision": "declared_chip_not_cross_validated_single_sample",
            }
        ]
        payload = {
            "schema_version": SCHEMA_VERSION,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "project_id": args.project_id,
            "method": METHOD,
            "mode": "single_sample_structural_qc",
            "barcode_pattern": "C4L15C4L15C4U10T18",
            "minimum_poly_t_bases": 12,
            "minimum_fold": args.minimum_fold,
            "primary_identity_metric": "exact_poly_t_rate",
            "cross_sample_assignment": "not_applicable",
            "input_gate_passed": single_sample_passed,
            "diagonal_exact_poly_t_rate_sum": exact_rates[(sample_id, sample["chip_id"])],
            "off_diagonal_exact_poly_t_rate_sum": None,
            "diagonal_over_off_diagonal_exact_fold": None,
            "diagonal_corrected_poly_t_rate_sum": corrected_rates[(sample_id, sample["chip_id"])],
            "off_diagonal_corrected_poly_t_rate_sum": None,
            "diagonal_over_off_diagonal_corrected_fold": None,
            "sample_decisions": sample_decisions,
            "passed": single_sample_passed,
            "failure_policy": "hard_stop_without_automatic_sample_sheet_mutation",
            "interpretation_note": "单样本仅完成声明芯片的结构性输入QC，不能验证芯片交换。",
        }
        write_tsv(table_path, result_rows, columns)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        partial_json = json_path.with_name(f".{json_path.name}.partial.{os.getpid()}")
        partial_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(partial_json, json_path)
        print(json.dumps(payload, indent=2))
        return 0 if single_sample_passed else 3

    sample_decisions = []
    diagonal_exact_sum = 0.0
    off_diagonal_exact_sum = 0.0
    diagonal_corrected_sum = 0.0
    off_diagonal_corrected_sum = 0.0
    for sample in samples:
        sample_id = sample["sample_id"]
        declared = sample["chip_id"]
        alternate = next(candidate["chip_id"] for candidate in samples if candidate["chip_id"] != declared)
        own_exact = exact_rates[(sample_id, declared)]
        other_exact = exact_rates[(sample_id, alternate)]
        own_corrected = corrected_rates[(sample_id, declared)]
        other_corrected = corrected_rates[(sample_id, alternate)]
        diagonal_exact_sum += own_exact
        off_diagonal_exact_sum += other_exact
        diagonal_corrected_sum += own_corrected
        off_diagonal_corrected_sum += other_corrected
        if (
            own_exact >= args.minimum_fold * other_exact
            and own_corrected >= other_corrected
        ):
            decision = "declared_mapping_supported"
        elif (
            other_exact >= args.minimum_fold * own_exact
            and other_corrected >= own_corrected
        ):
            decision = "swapped_mapping_supported"
        else:
            decision = "ambiguous"
        sample_decisions.append(
            {
                "sample_id": sample_id,
                "declared_chip_id": declared,
                "alternate_chip_id": alternate,
                "declared_exact_poly_t_rate": own_exact,
                "alternate_exact_poly_t_rate": other_exact,
                "declared_over_alternate_exact_fold": (
                    own_exact / other_exact if other_exact else None
                ),
                "declared_corrected_poly_t_rate": own_corrected,
                "alternate_corrected_poly_t_rate": other_corrected,
                "declared_over_alternate_corrected_fold": (
                    own_corrected / other_corrected if other_corrected else None
                ),
                "decision": decision,
            }
        )
    project_exact_fold = (
        diagonal_exact_sum / off_diagonal_exact_sum if off_diagonal_exact_sum else None
    )
    project_corrected_fold = (
        diagonal_corrected_sum / off_diagonal_corrected_sum
        if off_diagonal_corrected_sum
        else None
    )
    passed = (
        all(item["decision"] == "declared_mapping_supported" for item in sample_decisions)
        and (project_exact_fold is None or project_exact_fold >= args.minimum_fold)
        and diagonal_corrected_sum >= off_diagonal_corrected_sum
    )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_id": args.project_id,
        "method": METHOD,
        "barcode_pattern": "C4L15C4L15C4U10T18",
        "minimum_poly_t_bases": 12,
        "minimum_fold": args.minimum_fold,
        "primary_identity_metric": "exact_poly_t_rate",
        "hamming1_policy": "diagnostic_and_must_not_contradict_exact_direction",
        "diagonal_exact_poly_t_rate_sum": diagonal_exact_sum,
        "off_diagonal_exact_poly_t_rate_sum": off_diagonal_exact_sum,
        "diagonal_over_off_diagonal_exact_fold": project_exact_fold,
        "diagonal_corrected_poly_t_rate_sum": diagonal_corrected_sum,
        "off_diagonal_corrected_poly_t_rate_sum": off_diagonal_corrected_sum,
        "diagonal_over_off_diagonal_corrected_fold": project_corrected_fold,
        "sample_decisions": sample_decisions,
        "passed": passed,
        "failure_policy": "hard_stop_without_automatic_sample_sheet_mutation",
    }
    write_tsv(table_path, result_rows, columns)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    partial_json = json_path.with_name(f".{json_path.name}.partial.{os.getpid()}")
    partial_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(partial_json, json_path)
    print(json.dumps(payload, indent=2))
    return 0 if passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
