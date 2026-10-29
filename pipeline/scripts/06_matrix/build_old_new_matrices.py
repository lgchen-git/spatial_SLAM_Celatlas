#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def write_matrix(counts, out_prefix, features, barcodes):
    feature_index = {gene: idx for idx, gene in enumerate(features)}
    barcode_index = {barcode: idx for idx, barcode in enumerate(barcodes)}

    if not counts:
        matrix = sparse.coo_matrix((len(features), len(barcodes))).tocsr()
    else:
        rows = []
        cols = []
        data = []
        for (gene_id, barcode), value in counts.items():
            rows.append(feature_index[gene_id])
            cols.append(barcode_index[barcode])
            data.append(value)
        matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(features), len(barcodes))).tocsr()
    mmwrite(f"{out_prefix}.mtx", matrix)


def main():
    parser = argparse.ArgumentParser(description="Build old/new sparse matrices from molecule calls.")
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--sample", required=True, help="Sample identifier retained for CLI compatibility")
    parser.add_argument("--outdir", required=True, help="Stage 05 result root")
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    args = parser.parse_args()

    matrix_dir = Path(args.outdir) / "04_tool_output" / "matrix"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    feature_names = {}
    barcode_positions = {}
    counts = {"new": Counter(), "old": Counter()}
    usecols = ["assignment", "gene_id", "gene_name", "barcode", "x", "y"]
    for chunk in pd.read_csv(args.molecule_calls, sep="\t", usecols=usecols, chunksize=args.chunk_size):
        chunk = chunk.dropna(subset=["gene_id", "barcode"])
        chunk["gene_id"] = chunk["gene_id"].astype(str)
        chunk["barcode"] = chunk["barcode"].astype(str)
        for row in chunk[["gene_id", "gene_name"]].drop_duplicates("gene_id").itertuples(index=False):
            feature_names.setdefault(row.gene_id, row.gene_name)
        for row in chunk[["barcode", "x", "y"]].drop_duplicates("barcode").itertuples(index=False):
            barcode_positions.setdefault(row.barcode, (row.x, row.y))
        grouped = chunk.groupby(["assignment", "gene_id", "barcode"]).size()
        for (assignment, gene_id, barcode), value in grouped.items():
            counts[assignment][(gene_id, barcode)] += int(value)

    features = sorted(feature_names)
    barcodes = sorted(barcode_positions)
    feature_meta = pd.DataFrame(
        [(gene_id, feature_names[gene_id]) for gene_id in features], columns=["gene_id", "gene_name"]
    )
    barcode_meta = pd.DataFrame(
        [(barcode, *barcode_positions[barcode]) for barcode in barcodes], columns=["barcode", "x", "y"]
    )
    for assignment in ("new", "old"):
        write_matrix(counts[assignment], str(matrix_dir / f"{assignment}_umi"), features, barcodes)

    # The row order of features.tsv and barcodes.tsv is the matrix row/column contract.
    feature_meta.to_csv(matrix_dir / "features.tsv", sep="\t", index=False)
    barcode_meta.to_csv(matrix_dir / "barcodes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
