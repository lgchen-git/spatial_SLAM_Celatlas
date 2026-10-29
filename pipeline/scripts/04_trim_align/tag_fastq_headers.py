#!/usr/bin/env python3
"""Attach spatial barcode metadata to synchronized FASTQ read names."""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def gzip_reader(path: Path, threads: int):
    if threads <= 1 or shutil.which("pigz") is None:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
        return
    process = subprocess.Popen(
        ["pigz", "-dc", "-p", str(threads), str(path)], stdout=subprocess.PIPE
    )
    assert process.stdout is not None
    handle = io.TextIOWrapper(process.stdout, encoding="utf-8")
    try:
        yield handle
    finally:
        handle.close()
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)


@contextmanager
def gzip_writer(path: Path, threads: int):
    if threads <= 1 or shutil.which("pigz") is None:
        with gzip.open(path, "wt", encoding="utf-8", compresslevel=1) as handle:
            yield handle
        return
    with path.open("wb") as raw:
        process = subprocess.Popen(
            ["pigz", "-c", "-p", str(threads), "-1"], stdin=subprocess.PIPE, stdout=raw
        )
        assert process.stdin is not None
        handle = io.TextIOWrapper(process.stdin, encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()
            if process.wait() != 0:
                raise subprocess.CalledProcessError(process.returncode, process.args)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--barcode-table", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--progress-every", type=int, default=10_000_000)
    args = parser.parse_args()

    source_fastq = Path(args.fastq)
    source_table = Path(args.barcode_table)
    output = Path(args.output)
    stats_path = Path(args.stats)
    output.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(f".{output.name}.partial.{os.getpid()}")
    gzip_threads = max(1, min(4, args.threads // 3))
    records = 0

    try:
        with gzip_reader(source_fastq, gzip_threads) as fastq, gzip_reader(
            source_table, gzip_threads
        ) as table, gzip_writer(partial, gzip_threads) as target:
            columns = table.readline().rstrip("\n").split("\t")
            expected = ["read_name", "barcode", "umi", "x", "y"]
            if columns[:5] != expected:
                raise ValueError(f"unexpected barcode table header: {columns}")
            for row in table:
                fields = row.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(f"malformed barcode row after record {records}")
                header = fastq.readline()
                sequence = fastq.readline()
                plus = fastq.readline()
                quality = fastq.readline()
                if not quality:
                    raise ValueError("FASTQ ended before barcode table")
                header_text = header.rstrip("\n")
                first, separator, remainder = header_text.partition(" ")
                read_name = first[1:] if first.startswith("@") else first
                if read_name != fields[0]:
                    raise ValueError(
                        f"FASTQ/barcode order mismatch at record {records + 1}: "
                        f"{read_name!r} != {fields[0]!r}"
                    )
                tagged = (
                    f"@{read_name}|CB={fields[1]}|UB={fields[2]}|X={fields[3]}|Y={fields[4]}"
                )
                if separator:
                    tagged += f" {remainder}"
                target.write(tagged + "\n")
                target.write(sequence)
                target.write(plus)
                target.write(quality)
                records += 1
                if args.progress_every and records % args.progress_every == 0:
                    print(f"[tag_fastq_headers] records={records}", flush=True)
            if fastq.readline():
                raise ValueError("FASTQ contains records absent from barcode table")
        partial.replace(output)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise

    stats_path.write_text(
        json.dumps(
            {
                "records": records,
                "metadata_source": str(source_table.resolve()),
                "fastq_source": str(source_fastq.resolve()),
                "tag_format": "read_name|CB=<barcode>|UB=<umi>|X=<x>|Y=<y>",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
