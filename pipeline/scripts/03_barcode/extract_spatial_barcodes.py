#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import io
import json
import shutil
import subprocess
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

import pandas as pd


SEGMENT_LENGTHS = {"C": 0, "L": 0, "U": 0, "T": 0}


def parse_pattern(pattern: str):
    segments = []
    i = 0
    while i < len(pattern):
        token = pattern[i]
        if token not in SEGMENT_LENGTHS:
            raise ValueError(f"unsupported chemistry token: {token}")
        i += 1
        j = i
        while j < len(pattern) and pattern[j].isdigit():
            j += 1
        if j == i:
            raise ValueError(f"missing length after token: {token}")
        length = int(pattern[i:j])
        segments.append((token, length))
        i = j
    return segments


def read_fastq(handle):
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            raise ValueError("truncated fastq record")
        yield header.rstrip(), seq.rstrip(), plus.rstrip(), qual.rstrip()


def parse_fastq_list(value: str):
    paths = [item.strip() for item in value.split(",") if item.strip()]
    if not paths:
        raise ValueError("empty FASTQ list")
    return paths


def extract_barcode_umi(seq: str, segments):
    cursor = 0
    barcode_parts = []
    umi_parts = []
    poly_t = ""
    for token, length in segments:
        sub = seq[cursor : cursor + length]
        if len(sub) != length:
            return None, None, None
        if token == "C":
            barcode_parts.append(sub)
        elif token == "U":
            umi_parts.append(sub)
        elif token == "T":
            poly_t += sub
        cursor += length
    return "".join(barcode_parts), "".join(umi_parts), poly_t


def mean_quality(qual: str):
    if not qual:
        return 0.0
    return sum(ord(ch) - 33 for ch in qual) / len(qual)


def integer_median_from_counts(counts: Counter) -> int:
    total = sum(counts.values())
    if not total:
        return 0
    left_rank, right_rank = (total - 1) // 2, total // 2
    cumulative = 0
    left = right = 0
    for value in sorted(counts):
        cumulative += counts[value]
        if cumulative > left_rank and left == 0:
            left = value
        if cumulative > right_rank:
            right = value
            break
    return int((left + right) / 2)


@contextmanager
def gzip_text_reader(path: str | Path, threads: int):
    if threads <= 1 or shutil.which("pigz") is None:
        with gzip.open(path, "rt") as handle:
            yield handle
        return
    process = subprocess.Popen(["pigz", "-dc", "-p", str(threads), str(path)], stdout=subprocess.PIPE)
    assert process.stdout is not None
    handle = io.TextIOWrapper(process.stdout)
    try:
        yield handle
    finally:
        handle.close()
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)


@contextmanager
def gzip_text_writer(path: str | Path, threads: int):
    if threads <= 1 or shutil.which("pigz") is None:
        with gzip.open(path, "wt") as handle:
            yield handle
        return
    with Path(path).open("wb") as raw:
        process = subprocess.Popen(["pigz", "-c", "-p", str(threads)], stdin=subprocess.PIPE, stdout=raw)
        assert process.stdin is not None
        handle = io.TextIOWrapper(process.stdin)
        try:
            yield handle
        finally:
            handle.close()
            if process.wait() != 0:
                raise subprocess.CalledProcessError(process.returncode, process.args)


def main():
    parser = argparse.ArgumentParser(description="Extract spatial barcode/UMI and write filtered R2.")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--whitelist", required=True, help="CSV with x,y,barcode")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True, help="Stage 02 result root")
    parser.add_argument(
        "--filtered-r2-output",
        help="Explicit path for the large filtered R2 FASTQ (normally under the configured tmp root)",
    )
    parser.add_argument("--pattern", default="C4L15C4L15C4U10T18")
    parser.add_argument("--min-polyt", type=int, default=12)
    parser.add_argument("--gzip-threads", type=int, default=1)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    table_dir = outdir / "01_tables"
    metadata_dir = outdir / "05_metadata"
    table_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    segments = parse_pattern(args.pattern)
    r1_paths = parse_fastq_list(args.r1)
    r2_paths = parse_fastq_list(args.r2)
    if len(r1_paths) != len(r2_paths):
        raise ValueError(f"R1/R2 FASTQ count mismatch: {len(r1_paths)} != {len(r2_paths)}")

    whitelist = pd.read_csv(
        args.whitelist,
        header=None,
        names=["x", "y", "barcode"],
        dtype={"x": int, "y": int, "barcode": str},
    )
    whitelist["barcode"] = whitelist["barcode"].str.strip()
    barcode_to_pos = whitelist.set_index("barcode")[["x", "y"]].to_dict("index")
    whitelist_set = set(barcode_to_pos)

    filtered_r2 = (
        Path(args.filtered_r2_output)
        if args.filtered_r2_output
        else outdir / "04_tool_output" / f"{args.sample}.filtered_r2.fastq.gz"
    )
    filtered_r2.parent.mkdir(parents=True, exist_ok=True)
    barcode_table = table_dir / f"{args.sample}.read_barcode_umi.tsv.gz"
    qc_json = metadata_dir / f"{args.sample}.barcode_qc.json"

    counters = Counter()
    quality_sums = Counter()
    umi_lengths = Counter()
    polyt_lengths = Counter()

    with gzip_text_writer(filtered_r2, args.gzip_threads) as out_fq, gzip_text_writer(barcode_table, args.gzip_threads) as out_tsv:
        out_tsv.write("read_name\tbarcode\tumi\tx\ty\tmean_r1_qual\tmean_r2_qual\tpolyt_len\n")
        for r1_path, r2_path in zip(r1_paths, r2_paths):
            counters["input_fastq_pairs"] += 1
            with gzip_text_reader(r1_path, args.gzip_threads) as r1h, gzip_text_reader(r2_path, args.gzip_threads) as r2h:
                for r1, r2 in zip(read_fastq(r1h), read_fastq(r2h)):
                    counters["total_reads"] += 1
                    h1, s1, _, q1 = r1
                    h2, s2, p2, q2 = r2
                    read_name_1 = h1.split()[0][1:]
                    read_name_2 = h2.split()[0][1:]
                    if read_name_1 != read_name_2:
                        raise ValueError(f"read name mismatch: {read_name_1} != {read_name_2}")
                    barcode, umi, poly_t = extract_barcode_umi(s1, segments)
                    if not barcode or not umi:
                        counters["failed_pattern"] += 1
                        continue
                    if barcode not in whitelist_set:
                        counters["barcode_not_in_whitelist"] += 1
                        continue
                    if poly_t.count("T") < args.min_polyt:
                        counters["failed_poly_t"] += 1
                        continue
                    pos = barcode_to_pos[barcode]
                    counters["kept_reads"] += 1
                    qual_r1 = mean_quality(q1)
                    qual_r2 = mean_quality(q2)
                    quality_sums["r1"] += qual_r1
                    quality_sums["r2"] += qual_r2
                    umi_lengths[len(umi)] += 1
                    polyt_lengths[len(poly_t)] += 1
                    out_tsv.write(
                        f"{read_name_1}\t{barcode}\t{umi}\t{pos['x']}\t{pos['y']}\t{qual_r1:.2f}\t{qual_r2:.2f}\t{len(poly_t)}\n"
                    )
                    out_fq.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")

    kept_fraction = counters["kept_reads"] / counters["total_reads"] if counters["total_reads"] else 0.0
    qc = {
        "sample": args.sample,
        "pattern": args.pattern,
        "input_fastq_pairs": counters["input_fastq_pairs"],
        "total_reads": counters["total_reads"],
        "kept_reads": counters["kept_reads"],
        "kept_fraction": kept_fraction,
        "failed_pattern": counters["failed_pattern"],
        "barcode_not_in_whitelist": counters["barcode_not_in_whitelist"],
        "failed_poly_t": counters["failed_poly_t"],
        "mean_r1_quality": round(quality_sums["r1"] / counters["kept_reads"], 3) if counters["kept_reads"] else 0.0,
        "mean_r2_quality": round(quality_sums["r2"] / counters["kept_reads"], 3) if counters["kept_reads"] else 0.0,
        "median_umi_length": integer_median_from_counts(umi_lengths),
        "median_polyt_length": integer_median_from_counts(polyt_lengths),
        "barcode_table": str(barcode_table),
        "filtered_r2": str(filtered_r2),
    }
    qc_json.write_text(json.dumps(qc, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
