#!/usr/bin/env python3
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_attributes(field: str):
    attrs = {}
    for item in field.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def subtract_intervals(gene_span, exons):
    remaining = []
    cursor = gene_span[0]
    for start, end in exons:
        if cursor < start:
            remaining.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < gene_span[1]:
        remaining.append((cursor, gene_span[1]))
    return [(s, e) for s, e in remaining if e > s]


def main():
    parser = argparse.ArgumentParser(description="Build flattened exon/intron annotations from GTF.")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = {}
    exons_by_gene = defaultdict(list)

    with open(args.gtf, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            attr_map = parse_attributes(attrs)
            gene_id = attr_map.get("gene_id")
            gene_name = attr_map.get("gene_name", gene_id)
            if not gene_id:
                continue
            start0 = int(start) - 1
            end1 = int(end)
            if feature == "gene":
                genes[gene_id] = {
                    "chrom": chrom,
                    "start": start0,
                    "end": end1,
                    "strand": strand,
                    "gene_name": gene_name,
                }
            elif feature == "exon":
                exons_by_gene[gene_id].append((start0, end1))
                if gene_id not in genes:
                    genes[gene_id] = {
                        "chrom": chrom,
                        "start": start0,
                        "end": end1,
                        "strand": strand,
                        "gene_name": gene_name,
                    }
                else:
                    genes[gene_id]["start"] = min(genes[gene_id]["start"], start0)
                    genes[gene_id]["end"] = max(genes[gene_id]["end"], end1)

    exon_rows = []
    intron_rows = []
    gene_rows = []
    for gene_id, meta in genes.items():
        exons = merge_intervals(exons_by_gene.get(gene_id, []))
        introns = subtract_intervals((meta["start"], meta["end"]), exons)
        gene_rows.append(
            {
                "chrom": meta["chrom"],
                "start": meta["start"],
                "end": meta["end"],
                "gene_id": gene_id,
                "gene_name": meta["gene_name"],
                "strand": meta["strand"],
            }
        )
        for start, end in exons:
            exon_rows.append(
                {
                    "chrom": meta["chrom"],
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "gene_name": meta["gene_name"],
                    "strand": meta["strand"],
                    "feature_class": "exon",
                }
            )
        for start, end in introns:
            intron_rows.append(
                {
                    "chrom": meta["chrom"],
                    "start": start,
                    "end": end,
                    "gene_id": gene_id,
                    "gene_name": meta["gene_name"],
                    "strand": meta["strand"],
                    "feature_class": "intron",
                }
            )

    exon_df = pd.DataFrame(exon_rows).sort_values(["chrom", "start", "end", "gene_id"])
    intron_df = pd.DataFrame(intron_rows).sort_values(["chrom", "start", "end", "gene_id"])
    gene_df = pd.DataFrame(gene_rows).sort_values(["chrom", "start", "end", "gene_id"])

    exon_df.to_csv(outdir / "exons.tsv", sep="\t", index=False)
    intron_df.to_csv(outdir / "introns.tsv", sep="\t", index=False)
    gene_df.to_csv(outdir / "genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
