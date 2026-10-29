#!/usr/bin/env python3
"""Call Spatial SLAM molecules with unique genomic-site and base-quality control.

The BAM is coordinate sorted. Molecules are retained only while their assigned
gene can still receive another alignment, then flushed when the stream passes
the gene end. This bounds memory without weakening the molecule key.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import json
import multiprocessing as mp
import re
import shutil
import struct
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import pysam


BIN_SIZE = 20_000
BASES = "ACGT"
BASE_CODE = {base: index for index, base in enumerate("ACGTN")}
CODE_BASE = "ACGTN"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
SITE_STRUCT = struct.Struct("<IBBBB")  # genomic position, BQ, base, read position, lane
OUTPUT_COLUMNS = (
    "barcode", "x", "y", "gene_id", "gene_name", "umi", "supporting_reads",
    "raw_tc_events", "raw_t_covered", "unique_tc_events_bq0", "unique_t_covered_bq0",
    "duplicate_t_evaluations", "ambiguous_t_sites", "masked_snp_sites",
    "tc_events", "t_covered", "feature_class", "assignment",
)

ANNOTATIONS = None
GENE_NAMES = None
GENE_ENDS = None
SNP_MASK = None
SENSITIVITY_BQ = None
DEFAULT_BQ = None


def load_annotations(annotation_dir):
    exon_df = pd.read_csv(Path(annotation_dir) / "exons.tsv", sep="\t", dtype={"chrom": str})
    intron_df = pd.read_csv(Path(annotation_dir) / "introns.tsv", sep="\t", dtype={"chrom": str})
    annotations = defaultdict(lambda: defaultdict(list))
    gene_names = {}
    gene_ends = defaultdict(dict)
    for frame in (exon_df, intron_df):
        for row in frame.itertuples(index=False):
            chrom = str(row.chrom)
            gene_id = sys.intern(str(row.gene_id))
            gene_name = sys.intern(str(row.gene_name))
            gene_names.setdefault(gene_id, gene_name)
            gene_ends[chrom][gene_id] = max(gene_ends[chrom].get(gene_id, 0), int(row.end))
            interval = {
                "start": int(row.start), "end": int(row.end), "gene_id": gene_id,
                "strand": row.strand, "feature_class": row.feature_class,
            }
            for bin_id in range(interval["start"] // BIN_SIZE, (interval["end"] - 1) // BIN_SIZE + 1):
                annotations[chrom][bin_id].append(interval)
    return annotations, gene_names, gene_ends


def overlap_length(read_blocks, interval):
    return sum(
        max(0, min(end, interval["end"]) - max(start, interval["start"]))
        for start, end in read_blocks
    )


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
        score = overlap_length(read_blocks, interval)
        if score <= 0:
            continue
        entry = overlaps[interval["gene_id"]]
        entry[0] += score
        entry[1] |= 1 if interval["feature_class"] == "exon" else 2
        entry[2] = interval["strand"]
    if not overlaps:
        return None
    hits = sorted(overlaps.items(), key=lambda item: item[1][0], reverse=True)
    top_gene, top = hits[0]
    if len(hits) > 1 and hits[1][1][0] == top[0]:
        return {"gene_id": None, "feature_class": "ambiguous_gene"}
    return {
        "gene_id": top_gene,
        "strand": top[2],
        "feature_class": {1: "exon", 2: "intron"}.get(top[1], "ambiguous_feature"),
    }


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


def lane_code(query_name):
    original = query_name.split("|", 1)[0]
    for token in original.split(":"):
        match = re.fullmatch(r"L(\d{1,3})", token)
        if match:
            return min(int(match.group(1)), 254)
    return 255


def lane_label(code):
    return "unknown" if int(code) == 255 else f"L{int(code):02d}"


def load_snp_vcf(path):
    mask = defaultdict(set)
    records = []
    if not path:
        return mask, records
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[6] not in {"PASS", "."}:
                continue
            chrom, pos, _, ref, alt, _, filt, info, fmt, sample = fields[:10]
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
                continue
            values = dict(zip(fmt.split(":"), sample.split(":")))
            depth = int(values.get("DP", values.get("SDP", "0")) or 0)
            alt_depth = int(values.get("AD", "0") or 0)
            frequency = values.get("FREQ", "0").rstrip("%")
            try:
                frequency_value = float(frequency) / 100.0
            except ValueError:
                frequency_value = alt_depth / depth if depth else 0.0
            position0 = int(pos) - 1
            mask[str(chrom)].add(position0)
            records.append({
                "chrom": str(chrom), "position_1based": int(pos), "ref": ref, "alt": alt,
                "filter": filt, "depth_bq15": depth, "ref_depth": int(values.get("RD", "0") or 0),
                "alt_depth": alt_depth, "alt_fraction": frequency_value,
                "alt_forward": int(values.get("ADF", "0") or 0),
                "alt_reverse": int(values.get("ADR", "0") or 0),
                "strand_supported_both": int(values.get("ADF", "0") or 0) > 0 and int(values.get("ADR", "0") or 0) > 0,
            })
    return mask, records


def block_reference_cache(read, fasta):
    return [(start, end, fasta.fetch(read.reference_name, start, end).upper()) for start, end in read.get_blocks()]


def collect_aligned_observations(read, fasta, strand, base_counter, target_pos_counter, target_bq_counter, lane_counter):
    """Return compact transcript-oriented T-site observations for one read."""
    blocks = block_reference_cache(read, fasta)
    block_index = 0
    sequence = (read.query_sequence or "").upper()
    qualities = read.query_qualities
    lane = lane_code(read.query_name)
    packed = bytearray()
    raw_tc = 0
    raw_t = 0
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is None or rpos is None or qpos >= len(sequence):
            continue
        while block_index < len(blocks) and rpos >= blocks[block_index][1]:
            block_index += 1
        if block_index >= len(blocks):
            break
        start, end, ref_sequence = blocks[block_index]
        if rpos < start or rpos >= end:
            continue
        ref_base = ref_sequence[rpos - start]
        read_base = sequence[qpos]
        if strand == "-":
            ref_base = ref_base.translate(COMPLEMENT)
            read_base = read_base.translate(COMPLEMENT)
        if ref_base not in BASES or read_base not in BASES:
            continue
        quality = int(qualities[qpos]) if qualities is not None and qpos < len(qualities) else 0
        quality = min(quality, 255)
        read_position = min(int(qpos) + 1, 255)
        base_counter[(lane, ref_base, read_base, quality)] += 1
        if ref_base != "T":
            continue
        converted = read_base == "C"
        raw_t += 1
        raw_tc += int(converted)
        target_pos_counter[(lane, read_position, quality, int(converted))] += 1
        target_bq_counter[(lane, quality, int(converted))] += 1
        lane_counter[(lane, int(converted))] += 1
        packed.extend(SITE_STRUCT.pack(int(rpos), quality, BASE_CODE.get(read_base, 4), read_position, lane))
    return raw_tc, raw_t, bytes(packed)


def iter_sites(payload):
    for offset in range(0, len(payload), SITE_STRUCT.size):
        yield SITE_STRUCT.unpack_from(payload, offset)


def summarize_sites(payload, supporting_reads, chrom, thresholds, snp_mask):
    """Collapse repeated molecule observations by genomic position."""
    if supporting_reads == 1:
        selected = []
        for position, quality, base_code, read_position, lane in iter_sites(payload):
            selected.append((position, quality, base_code, read_position, lane, False))
        distinct_positions = len(selected)
    else:
        best = {}
        for position, quality, base_code, read_position, lane in iter_sites(payload):
            current = best.get(position)
            if current is None or quality > current[0]:
                best[position] = [quality, base_code, read_position, lane, False]
            elif quality == current[0] and base_code != current[1]:
                current[4] = True
        selected = [(position, *values) for position, values in best.items()]
        distinct_positions = len(best)
    ambiguous = sum(int(row[5]) for row in selected)
    nonambiguous = [row for row in selected if not row[5]]
    unique_bq0_t = len(nonambiguous)
    unique_bq0_tc = sum(int(CODE_BASE[row[2]] == "C") for row in nonambiguous)
    masked_positions = snp_mask.get(str(chrom), set())
    masked = sum(int(row[0] in masked_positions) for row in nonambiguous)
    by_threshold = {}
    for threshold in thresholds:
        kept = [row for row in nonambiguous if row[1] >= threshold and row[0] not in masked_positions]
        by_threshold[threshold] = {
            "t": len(kept),
            "tc": sum(int(CODE_BASE[row[2]] == "C") for row in kept),
        }
    return {
        "unique_bq0_t": unique_bq0_t,
        "unique_bq0_tc": unique_bq0_tc,
        "distinct_positions": distinct_positions,
        "ambiguous": ambiguous,
        "masked": masked,
        "by_threshold": by_threshold,
    }


def read_class(count):
    if count == 1:
        return "1"
    if count == 2:
        return "2"
    return ">=3"


def process_reference(task):
    index, reference, bam_path, fasta_path, part_dir, new_threshold, progress_every = task
    assert all(value is not None for value in (ANNOTATIONS, GENE_NAMES, GENE_ENDS, SNP_MASK, SENSITIVITY_BQ, DEFAULT_BQ))
    part_path = Path(part_dir) / f"{index:04d}.molecules.tsv.gz"
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)
    qc = Counter()
    duplication = defaultdict(Counter)
    k_distribution = Counter()
    t_distribution = Counter()
    base_counter = Counter()
    target_pos_counter = Counter()
    target_bq_counter = Counter()
    lane_counter = Counter()
    active_genes = {}
    active_heap = []
    barcode_cache = {}
    thresholds = sorted(set(SENSITIVITY_BQ) | {DEFAULT_BQ})
    chrom_gene_ends = GENE_ENDS.get(str(reference), {})

    def write_gene(gene_id, writer):
        molecules = active_genes.pop(gene_id, None)
        if not molecules:
            return
        for (barcode, x, y, umi), record in molecules.items():
            supporting_reads, raw_tc, raw_t, feature_mask, payload = record
            site = summarize_sites(payload, supporting_reads, reference, thresholds, SNP_MASK)
            selected = site["by_threshold"][DEFAULT_BQ]
            feature_class = {1: "exon", 2: "intron"}.get(feature_mask, "ambiguous_feature")
            assignment = "new" if selected["tc"] >= new_threshold else "old"
            row = [
                barcode, x, y, gene_id, GENE_NAMES.get(gene_id, gene_id), umi, supporting_reads,
                raw_tc, raw_t, site["unique_bq0_tc"], site["unique_bq0_t"],
                max(raw_t - site["distinct_positions"], 0), site["ambiguous"], site["masked"],
                selected["tc"], selected["t"], feature_class, assignment,
            ]
            for threshold in SENSITIVITY_BQ:
                row.extend([site["by_threshold"][threshold]["tc"], site["by_threshold"][threshold]["t"]])
            writer.writerow(row)
            category = read_class(supporting_reads)
            stats = duplication[category]
            stats["molecules"] += 1
            stats["raw_t_coverage"] += raw_t
            stats["raw_tc_events"] += raw_tc
            stats["unique_t_coverage_bq0"] += site["unique_bq0_t"]
            stats["unique_tc_events_bq0"] += site["unique_bq0_tc"]
            stats["duplicated_t_evaluations"] += max(raw_t - site["distinct_positions"], 0)
            stats["ambiguous_t_sites"] += site["ambiguous"]
            stats["masked_snp_sites"] += site["masked"]
            if raw_t > site["distinct_positions"]:
                stats["molecules_with_duplicated_sites"] += 1
            for threshold in thresholds:
                stats[f"t_bq{threshold}"] += site["by_threshold"][threshold]["t"]
                stats[f"tc_bq{threshold}"] += site["by_threshold"][threshold]["tc"]
            k_distribution[("raw_sum", raw_tc)] += 1
            k_distribution[(f"unique_bq{DEFAULT_BQ}_snp_masked", selected["tc"])] += 1
            t_distribution[("raw_sum", raw_t)] += 1
            t_distribution[(f"unique_bq{DEFAULT_BQ}_snp_masked", selected["t"])] += 1
            qc[f"molecules_{assignment}"] += 1
            qc["molecules_total"] += 1

    print(f"[reviewed_molecule_call] reference={reference} start", file=sys.stderr, flush=True)
    with gzip.open(part_path, "wt", encoding="utf-8", compresslevel=1, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for read in bam.fetch(reference):
            qc["total_bam_reads"] += 1
            if progress_every and qc["total_bam_reads"] % progress_every == 0:
                active_molecules = sum(len(value) for value in active_genes.values())
                print(
                    f"[reviewed_molecule_call] reference={reference} reads={qc['total_bam_reads']} active_molecules={active_molecules}",
                    file=sys.stderr, flush=True,
                )
            current_start = max(int(read.reference_start), 0)
            while active_heap and active_heap[0][0] <= current_start:
                _, finished_gene = heapq.heappop(active_heap)
                write_gene(finished_gene, writer)
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
            gene_id = assign["gene_id"]
            if gene_id not in active_genes:
                active_genes[gene_id] = {}
                heapq.heappush(active_heap, (chrom_gene_ends[gene_id], gene_id))
            raw_tc, raw_t, payload = collect_aligned_observations(
                read, fasta, assign["strand"], base_counter, target_pos_counter, target_bq_counter, lane_counter
            )
            barcode, umi, x, y = meta
            key = (barcode, x, y, umi)
            entry = active_genes[gene_id].get(key)
            if entry is None:
                entry = [0, 0, 0, 0, bytearray()]
                active_genes[gene_id][key] = entry
            entry[0] += 1
            entry[1] += raw_tc
            entry[2] += raw_t
            entry[3] |= {"exon": 1, "intron": 2}.get(assign["feature_class"], 4)
            entry[4].extend(payload)
            if assign["feature_class"] == "ambiguous_feature":
                qc["ambiguous_feature_reads"] += 1
        while active_heap:
            _, finished_gene = heapq.heappop(active_heap)
            write_gene(finished_gene, writer)

    bam.close()
    fasta.close()
    print(
        f"[reviewed_molecule_call] reference={reference} done reads={qc['total_bam_reads']} molecules={qc['molecules_total']}",
        file=sys.stderr, flush=True,
    )
    return (
        index, str(part_path), dict(qc), {key: dict(value) for key, value in duplication.items()},
        dict(k_distribution), dict(t_distribution), dict(base_counter), dict(target_pos_counter),
        dict(target_bq_counter), dict(lane_counter),
    )


def aggregate_counter_results(results, position):
    combined = Counter()
    for result in results:
        combined.update(result[position])
    return combined


def write_substitution_outputs(metadata_dir, results, thresholds):
    base_counts = aggregate_counter_results(results, 6)
    position_counts = aggregate_counter_results(results, 7)
    quality_counts = aggregate_counter_results(results, 8)
    lane_counts = aggregate_counter_results(results, 9)
    rows = []
    lanes = sorted({key[0] for key in base_counts})
    for threshold in thresholds:
        for scope, lane_values in [("global", lanes), *[(lane_label(lane), [lane]) for lane in lanes]]:
            for ref in BASES:
                denominator = sum(
                    count for (lane, observed_ref, _read, quality), count in base_counts.items()
                    if lane in lane_values and observed_ref == ref and quality >= threshold
                )
                for read_base in BASES:
                    if read_base == ref:
                        continue
                    mismatches = sum(
                        count for (lane, observed_ref, observed_read, quality), count in base_counts.items()
                        if lane in lane_values and observed_ref == ref and observed_read == read_base and quality >= threshold
                    )
                    rows.append({
                        "scope": scope, "min_base_quality": threshold, "channel": f"{ref}>{read_base}",
                        "reference_base_evaluations": denominator, "mismatches": mismatches,
                        "rate": mismatches / denominator if denominator else 0.0,
                        "orientation": "transcript_oriented",
                    })
    pd.DataFrame(rows).to_csv(metadata_dir / "spatial_substitution_spectrum.tsv", sep="\t", index=False)

    position_rows = []
    for threshold in thresholds:
        for position in sorted({key[1] for key in position_counts}):
            evaluations = sum(
                count for (lane, pos, quality, converted), count in position_counts.items()
                if pos == position and quality >= threshold
            )
            conversions = sum(
                count for (lane, pos, quality, converted), count in position_counts.items()
                if pos == position and quality >= threshold and converted == 1
            )
            position_rows.append({
                "min_base_quality": threshold, "read_position_1based": position,
                "t_evaluations": evaluations, "tc_events": conversions,
                "tc_rate": conversions / evaluations if evaluations else 0.0,
            })
    pd.DataFrame(position_rows).to_csv(metadata_dir / "tc_by_read_position.tsv", sep="\t", index=False)

    quality_rows = []
    for quality in sorted({key[1] for key in quality_counts}):
        evaluations = sum(count for (_lane, q, _converted), count in quality_counts.items() if q == quality)
        conversions = sum(count for (_lane, q, converted), count in quality_counts.items() if q == quality and converted == 1)
        quality_rows.append({
            "base_quality": quality, "t_evaluations": evaluations, "tc_events": conversions,
            "tc_rate": conversions / evaluations if evaluations else 0.0,
        })
    pd.DataFrame(quality_rows).to_csv(metadata_dir / "tc_by_base_quality.tsv", sep="\t", index=False)

    lane_rows = []
    for threshold in thresholds:
        for lane in sorted({key[0] for key in lane_counts}):
            evaluations = sum(
                count for (observed_lane, _converted), count in lane_counts.items()
                if observed_lane == lane
            )
            conversions = sum(
                count for (observed_lane, converted), count in lane_counts.items()
                if observed_lane == lane and converted == 1
            )
            # Exact-quality filtering is derived from the richer base counter.
            evaluations_bq = sum(
                count for (observed_lane, ref, _read, quality), count in base_counts.items()
                if observed_lane == lane and ref == "T" and quality >= threshold
            )
            conversions_bq = sum(
                count for (observed_lane, ref, read, quality), count in base_counts.items()
                if observed_lane == lane and ref == "T" and read == "C" and quality >= threshold
            )
            lane_rows.append({
                "lane": lane_label(lane), "min_base_quality": threshold,
                "t_evaluations": evaluations_bq, "tc_events": conversions_bq,
                "tc_rate": conversions_bq / evaluations_bq if evaluations_bq else 0.0,
                "raw_t_evaluations_bq0": evaluations, "raw_tc_events_bq0": conversions,
            })
    pd.DataFrame(lane_rows).to_csv(metadata_dir / "lane_tc_background.tsv", sep="\t", index=False)


def main():
    global ANNOTATIONS, GENE_NAMES, GENE_ENDS, SNP_MASK, SENSITIVITY_BQ, DEFAULT_BQ
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--annotation-dir", required=True)
    parser.add_argument("--barcode-table", required=True, help="Provenance input; metadata is carried in QNAME")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--new-threshold", type=int, default=1)
    parser.add_argument("--min-base-quality", type=int, default=20)
    parser.add_argument("--sensitivity-base-qualities", default="15,20,25,30")
    parser.add_argument("--snp-vcf")
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()

    SENSITIVITY_BQ = sorted({int(value) for value in args.sensitivity_base_qualities.split(",") if value.strip()})
    DEFAULT_BQ = int(args.min_base_quality)
    if DEFAULT_BQ not in SENSITIVITY_BQ:
        SENSITIVITY_BQ.append(DEFAULT_BQ)
        SENSITIVITY_BQ.sort()
    outdir = Path(args.outdir)
    table_dir = outdir / "01_tables"
    metadata_dir = outdir / "05_metadata"
    part_dir = Path(args.tmp_dir)
    table_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    if part_dir.exists():
        shutil.rmtree(part_dir)
    part_dir.mkdir(parents=True)

    print("[reviewed_molecule_call] loading annotations and SNP mask", file=sys.stderr, flush=True)
    ANNOTATIONS, GENE_NAMES, GENE_ENDS = load_annotations(args.annotation_dir)
    SNP_MASK, snp_records = load_snp_vcf(args.snp_vcf)
    if snp_records:
        pd.DataFrame(snp_records).to_csv(metadata_dir / "recurrent_variant_candidates.tsv", sep="\t", index=False)
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        references = [stat.contig for stat in bam.get_index_statistics() if stat.mapped > 0]
        unmapped = int(bam.unmapped)
    workers = max(1, min(args.threads, len(references)))
    tasks = [
        (index, reference, args.bam, args.fasta, str(part_dir), args.new_threshold, args.progress_every)
        for index, reference in enumerate(references)
    ]
    context = mp.get_context("fork")
    with context.Pool(processes=workers, maxtasksperchild=1) as pool:
        results = list(pool.imap_unordered(process_reference, tasks))
    results.sort(key=lambda item: item[0])

    molecule_path = table_dir / f"{args.sample}.molecule_calls.tsv.gz"
    partial = molecule_path.with_name(f".{molecule_path.name}.partial")
    sensitivity_columns = tuple(
        column for threshold in SENSITIVITY_BQ
        for column in (f"tc_events_bq{threshold}", f"t_covered_bq{threshold}")
    )
    with gzip.open(partial, "wt", encoding="utf-8", compresslevel=1, newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow((*OUTPUT_COLUMNS, *sensitivity_columns))
    with partial.open("ab") as target:
        for _, part_path, *_rest in results:
            with Path(part_path).open("rb") as source:
                shutil.copyfileobj(source, target, length=16 * 1024 * 1024)
    partial.replace(molecule_path)

    qc = Counter({"unmapped_bam_reads": unmapped})
    duplication = defaultdict(Counter)
    k_distribution = Counter()
    t_distribution = Counter()
    for result in results:
        qc.update(result[2])
        for category, values in result[3].items():
            duplication[category].update(values)
        k_distribution.update(result[4])
        t_distribution.update(result[5])

    total_molecules = qc["molecules_total"]
    duplication_rows = []
    for category in ("1", "2", ">=3"):
        values = duplication[category]
        duplication_rows.append({
            "supporting_reads_class": category,
            "molecules": values["molecules"],
            "molecule_fraction": values["molecules"] / total_molecules if total_molecules else 0.0,
            **{key: value for key, value in sorted(values.items()) if key != "molecules"},
        })
    pd.DataFrame(duplication_rows).to_csv(table_dir / "molecule_site_duplication_qc.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"method": method, "tc_events": value, "molecules": count} for (method, value), count in sorted(k_distribution.items())]
    ).to_csv(table_dir / "molecule_k_distribution.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"method": method, "t_covered": value, "molecules": count} for (method, value), count in sorted(t_distribution.items())]
    ).to_csv(table_dir / "molecule_t_distribution.tsv", sep="\t", index=False)
    write_substitution_outputs(metadata_dir, results, SENSITIVITY_BQ)

    totals = Counter()
    for values in duplication.values():
        totals.update(values)
    summary = {
        "sample": args.sample,
        "counting_method": "unique_genomic_site_highest_base_quality",
        "molecule_key": ["barcode", "x", "y", "gene_id", "umi"],
        "min_base_quality": DEFAULT_BQ,
        "sensitivity_base_qualities": SENSITIVITY_BQ,
        "snp_mask_source": str(Path(args.snp_vcf).resolve()) if args.snp_vcf else None,
        "snp_mask_records": len(snp_records),
        "molecules": total_molecules,
        "one_read_molecules": duplication["1"]["molecules"],
        "one_read_molecule_fraction": duplication["1"]["molecules"] / total_molecules if total_molecules else 0.0,
        "two_read_molecules": duplication["2"]["molecules"],
        "two_read_molecule_fraction": duplication["2"]["molecules"] / total_molecules if total_molecules else 0.0,
        "three_or_more_read_molecules": duplication[">=3"]["molecules"],
        "three_or_more_read_molecule_fraction": duplication[">=3"]["molecules"] / total_molecules if total_molecules else 0.0,
        "molecules_with_duplicated_sites": totals["molecules_with_duplicated_sites"],
        "molecules_with_duplicated_sites_fraction": totals["molecules_with_duplicated_sites"] / total_molecules if total_molecules else 0.0,
        "raw_t_coverage": totals["raw_t_coverage"],
        "unique_site_t_coverage_bq0": totals["unique_t_coverage_bq0"],
        "raw_tc_events": totals["raw_tc_events"],
        "unique_site_tc_events_bq0": totals["unique_tc_events_bq0"],
        "duplicated_t_evaluations": totals["duplicated_t_evaluations"],
        "ambiguous_t_sites": totals["ambiguous_t_sites"],
        "masked_snp_sites": totals["masked_snp_sites"],
        "selected_t_coverage": totals[f"t_bq{DEFAULT_BQ}"],
        "selected_tc_events": totals[f"tc_bq{DEFAULT_BQ}"],
        "worker_processes": workers,
        "streaming_strategy": "coordinate_sorted_bam_flush_after_gene_end",
    }
    (metadata_dir / "molecule_site_duplication_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    qc.update({
        "metadata_source_qname": 1,
        "worker_processes": workers,
        "min_base_quality": DEFAULT_BQ,
        "snp_mask_records": len(snp_records),
        "unique_site_counting": 1,
    })
    (metadata_dir / f"{args.sample}.qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    shutil.rmtree(part_dir)


if __name__ == "__main__":
    main()
