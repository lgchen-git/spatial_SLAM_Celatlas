#!/usr/bin/env python3
"""Retain ordered barcode-extracted R2 reads whose barcode belongs to tissue."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def read_tissue_barcodes(tissue_barcodes):
    with open(tissue_barcodes, encoding="utf-8", newline="") as handle:
        return {row["barcode"] for row in csv.DictReader(handle, delimiter="\t")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--barcode-table", required=True)
    parser.add_argument("--tissue-barcodes", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", help="Optional JSON summary written after a complete validated stream")
    args = parser.parse_args()
    tissue = read_tissue_barcodes(args.tissue_barcodes)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    total = retained = 0
    # The barcode table and filtered FASTQ are emitted in the same read order by
    # stage 02. Streaming them in lockstep avoids retaining tens of millions of
    # read names in memory. A mismatch is fatal rather than silently misfiltering.
    with gzip.open(args.barcode_table, "rt") as barcode_handle, gzip.open(args.fastq, "rt") as source, gzip.open(output, "wt", compresslevel=1) as target:
        rows = csv.DictReader(barcode_handle, delimiter="\t")
        for row in rows:
            header = source.readline()
            if not header:
                raise SystemExit("FASTQ ended before barcode table")
            seq, plus, qual = source.readline(), source.readline(), source.readline()
            if not qual:
                raise SystemExit("Truncated FASTQ")
            read_name = header.split()[0][1:]
            if read_name != row["read_name"]:
                raise SystemExit(f"Barcode/FASTQ order mismatch: {row['read_name']} != {read_name}")
            total += 1
            if row["barcode"] in tissue:
                target.write(header + seq + plus + qual)
                retained += 1
        if source.readline():
            raise SystemExit("FASTQ contains more reads than barcode table")
    print(f"tissue_fastq_reads={retained} total_barcode_reads={total}")
    if args.summary:
        summary = Path(args.summary)
        summary.parent.mkdir(parents=True, exist_ok=True)
        summary.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "total_barcode_reads": total,
                    "tissue_fastq_reads": retained,
                    "retained_fraction": retained / total if total else 0.0,
                    "tissue_barcode_count": len(tissue),
                    "stream_order_validated": True,
                    "input_fastq": str(Path(args.fastq).resolve()),
                    "input_barcode_table": str(Path(args.barcode_table).resolve()),
                    "input_tissue_barcodes": str(Path(args.tissue_barcodes).resolve()),
                    "output_fastq": str(output.resolve()),
                    "output_fastq_bytes": output.stat().st_size,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
