#!/usr/bin/env python3
"""Call old/new molecules from metadata-tagged alignments with bounded memory."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import multiprocessing as mp
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam


BIN_SIZE = 20_000
OUTPUT_COLUMNS = (
    "barcode", "x", "y", "gene_id", "gene_name", "umi", "supporting_reads",
    "tc_events", "t_covered", "feature_class", "assignment",
)
ANNOTATIONS = None
GENE_NAMES = None


def load_annotations(annotation_dir):
    exon_df = pd.read_csv(Path(annotation_dir) / "exons.tsv", sep="\t", dtype={"chrom": str})
    intron_df = pd.read_csv(Path(annotation_dir) / "introns.tsv", sep="\t", dtype={"chrom": str})
    annotations = defaultdict(lambda: defaultdict(list))
    gene_names = {}
    for frame in (exon_df, intron_df):
        for row in frame.itertuples(index=False):
            gene_id = sys.intern(str(row.gene_id))
            gene_name = sys.intern(str(row.gene_name))
            gene_names.setdefault(gene_id, gene_name)
            interval = {
                "start": int(row.start), "end": int(row.end), "gene_id": gene_id,
                "strand": row.strand, "feature_class": row.feature_class,
            }
            for bin_id in range(interval["start"] // BIN_SIZE, (interval["end"] - 1) // BIN_SIZE + 1):
                annotations[str(row.chrom)][bin_id].append(interval)
    return annotations, gene_names


def overlap_fraction(read_blocks, interval):
    covered = 0
    for start, end in read_blocks:
        ov_start = max(start, interval["start"])
        ov_end = min(end, interval["end"])
        if ov_end > ov_start:
            covered += ov_end - ov_start
    return covered


def assign_gene_feature(chrom, read_blocks, annotations):
    overlaps = defaultdict(lambda: [0, 0, ""])
    candidates = {}
    chrom_bins = annotations.get(str(chrom))
    if not chrom_bins:
        return None
    for start, end in read_blocks:
        for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
            for interval in chrom_bins.get(bin_id, []):
                candidates[id(interval)] = interval
    for interval in candidates.values():
        score = overlap_fraction(read_blocks, interval)
        if score <= 0:
            continue
        entry = overlaps[interval["gene_id"]]
        entry[0] += score
        entry[1] |= 1 if interval["feature_class"] == "exon" else 2
        entry[2] = interval["strand"]
    if not overlaps:
        return None
    sorted_hits = sorted(overlaps.items(), key=lambda item: item[1][0], reverse=True)
    top_gene, top = sorted_hits[0]
    if len(sorted_hits) > 1 and sorted_hits[1][1][0] == top[0]:
        return {"gene_id": None, "feature_class": "ambiguous_gene"}
    feature_class = {1: "exon", 2: "intron"}.get(top[1], "ambiguous_feature")
    return {"gene_id": top_gene, "strand": top[2], "feature_class": feature_class}


def block_reference_cache(read, fasta):
    return [(start, end, fasta.fetch(read.reference_name, start, end).upper()) for start, end in read.get_blocks()]


def count_tc_events(read, fasta, strand):
    tc_events = 0
    t_covered = 0
    block_cache = block_reference_cache(read, fasta)
    block_index = 0
    query_sequence = (read.query_sequence or "").upper()
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is None or rpos is None:
            continue
        while block_index < len(block_cache) and rpos >= block_cache[block_index][1]:
            block_index += 1
        if block_index >= len(block_cache):
            break
        block_start, block_end, ref_seq = block_cache[block_index]
        if rpos < block_start or rpos >= block_end or qpos >= len(query_sequence):
            continue
        ref_base = ref_seq[rpos - block_start]
        read_base = query_sequence[qpos]
        if strand == "+" and ref_base == "T":
            t_covered += 1
            tc_events += read_base == "C"
        elif strand == "-" and ref_base == "A":
            t_covered += 1
            tc_events += read_base == "G"
    return int(tc_events), t_covered


def parse_tagged_metadata(query_name, barcode_cache):
    fields = query_name.rsplit("|", 4)
    if len(fields) != 5:
        return None
    tags = {}
    for field in fields[1:]:
        key, separator, value = field.partition("=")
        if not separator:
            return None
        tags[key] = value
    if set(tags) != {"CB", "UB", "X", "Y"}:
        return None
    barcode = sys.intern(tags["CB"])
    coordinates = barcode_cache.get(barcode)
    if coordinates is None:
        coordinates = (int(tags["X"]), int(tags["Y"]))
        barcode_cache[barcode] = coordinates
    return barcode, sys.intern(tags["UB"]), coordinates[0], coordinates[1]


def process_reference(task):
    index, reference, bam_path, fasta_path, part_dir, threshold, progress_every = task
    assert ANNOTATIONS is not None and GENE_NAMES is not None
    part_path = Path(part_dir) / f"{index:04d}.molecules.tsv.gz"
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)
    qc = Counter()
    feature_stats = defaultdict(lambda: {"tc_events": 0, "t_covered": 0, "reads": 0})
    molecule_reads = {}
    barcode_cache = {}
    print(f"[call_new_old] reference={reference} start", file=sys.stderr, flush=True)

    for read in bam.fetch(reference):
        qc["total_bam_reads"] += 1
        if progress_every and qc["total_bam_reads"] % progress_every == 0:
            print(
                f"[call_new_old] reference={reference} reads={qc['total_bam_reads']} molecules={len(molecule_reads)}",
                file=sys.stderr, flush=True,
            )
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            qc["skipped_non_primary"] += 1
            continue
        meta = parse_tagged_metadata(read.query_name, barcode_cache)
        if meta is None:
            qc["missing_barcode_metadata"] += 1
            continue
        assign = assign_gene_feature(read.reference_name, read.get_blocks(), ANNOTATIONS)
        if assign is None:
            qc["unassigned_gene"] += 1
            continue
        if assign.get("gene_id") is None:
            qc["ambiguous_gene"] += 1
            continue
        tc_events, t_covered = count_tc_events(read, fasta, assign["strand"])
        feature_class = assign["feature_class"]
        barcode, umi, x, y = meta
        key = (barcode, x, y, assign["gene_id"], umi)
        entry = molecule_reads.get(key)
        if entry is None:
            entry = [0, 0, 0, 0]
            molecule_reads[key] = entry
        entry[0] += tc_events
        entry[1] += t_covered
        entry[2] += 1
        entry[3] |= {"exon": 1, "intron": 2}.get(feature_class, 4)
        if feature_class in {"exon", "intron"}:
            feature_stats[feature_class]["tc_events"] += tc_events
            feature_stats[feature_class]["t_covered"] += t_covered
            feature_stats[feature_class]["reads"] += 1
        else:
            qc["ambiguous_feature_reads"] += 1

    with gzip.open(part_path, "wt", encoding="utf-8", compresslevel=1, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for key, record in molecule_reads.items():
            barcode, x, y, gene_id, umi = key
            feature_class = {1: "exon", 2: "intron"}.get(record[3], "ambiguous_feature")
            assignment = "new" if record[0] >= threshold else "old"
            writer.writerow((
                barcode, x, y, gene_id, GENE_NAMES.get(gene_id, gene_id), umi,
                record[2], record[0], record[1], feature_class, assignment,
            ))
            qc[f"molecules_{assignment}"] += 1
    qc["molecules_total"] = len(molecule_reads)
    bam.close()
    fasta.close()
    print(
        f"[call_new_old] reference={reference} done reads={qc['total_bam_reads']} molecules={qc['molecules_total']}",
        file=sys.stderr, flush=True,
    )
    return index, str(part_path), dict(qc), {key: dict(value) for key, value in feature_stats.items()}


def main():
    global ANNOTATIONS, GENE_NAMES
    parser = argparse.ArgumentParser(description="Call old/new molecules and exon/intron conversion.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--barcode-table", required=True, help="Provenance input; metadata is carried in QNAME")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--new-threshold", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    table_dir = outdir / "01_tables"
    metadata_dir = outdir / "05_metadata"
    part_dir = Path(args.tmp_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    if part_dir.exists():
        shutil.rmtree(part_dir)
    part_dir.mkdir(parents=True)

    print("[call_new_old] loading annotations", file=sys.stderr, flush=True)
    ANNOTATIONS, GENE_NAMES = load_annotations(args.annotation_dir)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        references = [stat.contig for stat in bam.get_index_statistics() if stat.mapped > 0]
        unmapped = int(bam.unmapped)
    workers = max(1, min(args.threads, len(references)))
    tasks = [
        (index, reference, args.bam, args.fasta, str(part_dir), args.new_threshold, args.progress_every)
        for index, reference in enumerate(references)
    ]
    print(
        f"[call_new_old] references={len(references)} workers={workers} metadata=qname",
        file=sys.stderr, flush=True,
    )
    context = mp.get_context("fork")
    with context.Pool(processes=workers, maxtasksperchild=1) as pool:
        results = list(pool.imap_unordered(process_reference, tasks))
    results.sort(key=lambda item: item[0])

    qc = Counter({"unmapped_bam_reads": unmapped})
    feature_stats = defaultdict(lambda: {"tc_events": 0, "t_covered": 0, "reads": 0})
    molecule_path = table_dir / f"{args.sample}.molecule_calls.tsv.gz"
    partial = molecule_path.with_name(f".{molecule_path.name}.partial")
    with gzip.open(partial, "wt", encoding="utf-8", compresslevel=1, newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(OUTPUT_COLUMNS)
    with partial.open("ab") as target:
        for _, part_path, part_qc, part_features in results:
            qc.update(part_qc)
            for feature_class, values in part_features.items():
                for key, value in values.items():
                    feature_stats[feature_class][key] += value
            with Path(part_path).open("rb") as source:
                shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    partial.replace(molecule_path)

    rate_path = table_dir / f"{args.sample}.exon_intron_conversion.summary.tsv"
    with rate_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample", "feature_class", "reads", "t_covered", "tc_events", "tc_rate"),
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        for feature_class in ("exon", "intron"):
            values = feature_stats[feature_class]
            writer.writerow({
                "sample": args.sample, "feature_class": feature_class, **values,
                "tc_rate": values["tc_events"] / values["t_covered"] if values["t_covered"] else 0.0,
            })
    qc["metadata_source_qname"] = 1
    qc["worker_processes"] = workers
    (metadata_dir / f"{args.sample}.qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    shutil.rmtree(part_dir)


if __name__ == "__main__":
    main()
