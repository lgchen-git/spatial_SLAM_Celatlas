#!/usr/bin/env python3
"""Validate the generic sample and library sheets before a run."""
import argparse
import csv
from pathlib import Path


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--libraries", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    samples = [row for row in read_tsv(args.samples) if row["sample_id"] == args.sample]
    if len(samples) != 1:
        raise SystemExit(f"Expected exactly one sample row for {args.sample}, found {len(samples)}")
    sample = samples[0]
    libraries = [row for row in read_tsv(args.libraries) if row["sample_id"] == args.sample]
    if not libraries:
        raise SystemExit(f"No libraries found for {args.sample}")
    required = [
        Path(sample["mask_dir"]) / f'{sample["chip_id"]}.barcodeToPos.h5',
        Path(sample["mask_dir"]) / f'{sample["chip_id"]}_FilterBarcodes.csv',
        Path(sample["mask_dir"]) / f'{sample["chip_id"]}_tissue_bbox.csv',
        Path(sample["reference_dir"]) / "genome.fa",
        Path(sample["reference_dir"]) / "genes.gtf",
        Path(sample["reference_dir"]) / "Genome",
    ]
    for row in libraries:
        required.extend([Path(row["fastq_r1"]), Path(row["fastq_r2"])])
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        raise SystemExit("Missing required files:\n" + "\n".join(absent))
    if not sample["coordinate_um_per_unit"]:
        raise SystemExit("coordinate_um_per_unit must be set")
    labeling_status = (sample.get("labeling_status") or "unknown").strip()
    allowed_statuses = {"unknown", "labeled", "unlabeled_control"}
    if labeling_status not in allowed_statuses:
        raise SystemExit(
            f"labeling_status must be one of {sorted(allowed_statuses)}, got {labeling_status!r}"
        )
    celatlas_result_dir = (sample.get("celatlas_result_dir") or "").strip()
    if celatlas_result_dir and not Path(celatlas_result_dir).is_dir():
        raise SystemExit(f"celatlas_result_dir is not a directory: {celatlas_result_dir}")
    print(
        f"validated sample={args.sample} libraries={len(libraries)} "
        f"labeling_status={labeling_status} celatlas_reuse={'yes' if celatlas_result_dir else 'no'}"
    )


if __name__ == "__main__":
    main()
