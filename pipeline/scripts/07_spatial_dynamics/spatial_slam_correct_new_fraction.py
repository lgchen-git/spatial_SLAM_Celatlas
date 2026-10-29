#!/usr/bin/env python3
import argparse
import gzip
import io
import json
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.special import gammaln


def write_tsv_gzip(frame, path):
    """Write TSV through parallel pigz when available, with a portable fallback."""
    path = Path(path)
    pigz = shutil.which("pigz")
    if not pigz:
        frame.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 1})
        return
    threads = max(1, min(8, os.cpu_count() or 1))
    with path.open("wb") as output:
        process = subprocess.Popen([pigz, "-1", "-p", str(threads), "-c"], stdin=subprocess.PIPE, stdout=output)
        assert process.stdin is not None
        with io.TextIOWrapper(process.stdin, encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, sep="\t", index=False)
        if process.wait() != 0:
            raise RuntimeError("pigz failed while writing corrected bin-gene counts")


def write_corrected_frame(frame, args):
    if args.output_pickle:
        Path(args.output_pickle).parent.mkdir(parents=True, exist_ok=True)
        frame.to_pickle(args.output_pickle)
    if not args.skip_tsv:
        write_tsv_gzip(frame, Path(args.outdir) / f"{args.sample}.corrected_bin_gene_counts.tsv.gz")


def parse_args():
    parser = argparse.ArgumentParser(description="Correct Spatial SLAM new RNA fractions with a binomial mixture model.")
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--matrix-dir")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bin-size", type=int, default=50)
    parser.add_argument("--bin-size-pixels", type=int, help="Raw pixel width of the requested physical bin")
    parser.add_argument("--official-bin-positions", help="Official Celatlas tissue_positions_list.csv used as the bin grid")
    parser.add_argument("--background-conversion-rate", type=float)
    parser.add_argument("--marked-conversion-rate", type=float)
    parser.add_argument("--initial-pi", type=float)
    parser.add_argument("--initial-p-old", type=float)
    parser.add_argument("--initial-p-new", type=float)
    parser.add_argument("--tc-column", default="tc_events")
    parser.add_argument("--t-column", default="t_covered")
    parser.add_argument("--chunk-size", type=int, default=2_000_000)
    parser.add_argument("--max-em-iter", type=int, default=120)
    parser.add_argument("--tol", type=float, default=1e-7)
    parser.add_argument("--write-corrected-molecules", action="store_true")
    parser.add_argument(
        "--max-estimate-rows",
        type=int,
        default=None,
        help="Optional deterministic prefix cap for model fitting; default scans all molecules",
    )
    parser.add_argument("--output-pickle", help="Optional fast internal DataFrame cache for downstream analysis")
    parser.add_argument("--skip-tsv", action="store_true", help="Skip portable TSV output; intended only for time-critical internal runs")
    return parser.parse_args()


def load_official_grid(path):
    positions = pd.read_csv(path, header=None, names=["barcode", "in_tissue", "row", "col", "pxl_row_in_fullres", "pxl_col_in_fullres"], dtype={"barcode": str})
    positions["row"] = positions["row"].astype(int)
    positions["col"] = positions["col"].astype(int)
    positions["bin_x"] = positions["pxl_col_in_fullres"].astype(float)
    positions["bin_y"] = positions["pxl_row_in_fullres"].astype(float)
    return positions.set_index(["col", "row"])


def assign_bins(frame, args, official_grid):
    pixel_size = int(args.bin_size_pixels or args.bin_size)
    frame["bin_col"] = np.floor(frame["x"].to_numpy(dtype=float) / pixel_size).astype(int)
    frame["bin_row"] = np.floor(frame["y"].to_numpy(dtype=float) / pixel_size).astype(int)
    if official_grid is None:
        frame["bin_x"] = frame["bin_col"].astype(float) * pixel_size + pixel_size / 2.0
        frame["bin_y"] = frame["bin_row"].astype(float) * pixel_size + pixel_size / 2.0
        frame["bin_id"] = frame["bin_x"].astype(int).astype(str) + "_" + frame["bin_y"].astype(int).astype(str)
        return frame
    keys = pd.MultiIndex.from_arrays([frame["bin_col"], frame["bin_row"]])
    frame = frame.loc[keys.isin(official_grid.index)].copy()
    keys = pd.MultiIndex.from_arrays([frame["bin_col"], frame["bin_row"]])
    mapped = official_grid.reindex(keys)
    frame["bin_id"] = mapped["barcode"].to_numpy()
    frame["bin_x"] = mapped["bin_x"].to_numpy()
    frame["bin_y"] = mapped["bin_y"].to_numpy()
    return frame


def read_pair_counts(path, chunk_size, max_rows=None, tc_column="tc_events", t_column="t_covered"):
    pair_counts = Counter()
    total = 0
    informative = 0
    total_tc = 0
    total_t = 0
    usecols = [tc_column, t_column]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=chunk_size):
        if max_rows is not None and total >= max_rows:
            break
        if max_rows is not None and total + len(chunk) > max_rows:
            chunk = chunk.iloc[: max_rows - total].copy()
        chunk[tc_column] = chunk[tc_column].fillna(0).astype(int)
        chunk[t_column] = chunk[t_column].fillna(0).astype(int)
        total += len(chunk)
        sub = chunk[chunk[t_column] > 0]
        informative += len(sub)
        total_tc += int(sub[tc_column].sum())
        total_t += int(sub[t_column].sum())
        counts = sub.groupby([tc_column, t_column]).size()
        for (k, t), value in counts.items():
            pair_counts[(int(k), int(t))] += int(value)
    return pair_counts, {"molecules_total": total, "molecules_informative": informative, "tc_events": total_tc, "t_covered": total_t}


def estimate_assignment_posteriors(path, chunk_size, params, max_rows, tc_column="tc_events", t_column="t_covered"):
    usecols = [tc_column, t_column, "assignment"]
    total = 0
    sums = Counter()
    counts = Counter()
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=chunk_size):
        if max_rows is not None and total >= max_rows:
            break
        if max_rows is not None and total + len(chunk) > max_rows:
            chunk = chunk.iloc[: max_rows - total].copy()
        total += len(chunk)
        posterior = posterior_new_array(
            chunk[tc_column].fillna(0).astype(int).to_numpy(),
            chunk[t_column].fillna(0).astype(int).to_numpy(),
            params["pi"],
            params["p_old"],
            params["p_new"],
        )
        for assignment, values in pd.DataFrame({"assignment": chunk["assignment"].to_numpy(), "posterior": posterior}).groupby("assignment"):
            sums[str(assignment)] += float(values["posterior"].sum())
            counts[str(assignment)] += int(len(values))
    means = {}
    for assignment in ("old", "new"):
        if counts[assignment]:
            means[assignment] = sums[assignment] / counts[assignment]
        else:
            means[assignment] = 0.0 if assignment == "old" else 1.0
    return {"rows": total, "posterior_mean_by_hard_assignment": means}


def log_binom_pmf(k, t, p):
    p = float(np.clip(p, 1e-8, 1 - 1e-8))
    return gammaln(t + 1) - gammaln(k + 1) - gammaln(t - k + 1) + k * np.log(p) + (t - k) * np.log1p(-p)


def posterior_new_array(k, t, pi, p_old, p_new):
    k = np.asarray(k, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    # Molecules without covered thymidines are uninformative about labeling;
    # their posterior must remain at the fitted prior instead of being forced old.
    post = np.full_like(k, float(pi), dtype=np.float64)
    informative = t > 0
    if not np.any(informative):
        return post
    log_old = np.log1p(-pi) + log_binom_pmf(k[informative], t[informative], p_old)
    log_new = np.log(pi) + log_binom_pmf(k[informative], t[informative], p_new)
    max_log = np.maximum(log_old, log_new)
    old_prob = np.exp(log_old - max_log)
    new_prob = np.exp(log_new - max_log)
    post[informative] = new_prob / np.maximum(old_prob + new_prob, 1e-12)
    return post


def fit_mixture(
    pair_counts, user_p_old=None, user_p_new=None, max_iter=120, tol=1e-7,
    initial_pi=None, initial_p_old=None, initial_p_new=None,
):
    rows = np.array([(k, t, c) for (k, t), c in pair_counts.items()], dtype=np.float64)
    if rows.size == 0:
        return {"pi": 0.0, "p_old": 1e-6, "p_new": 1e-3, "iterations": 0, "converged": True}
    k = rows[:, 0]
    t = rows[:, 1]
    counts = rows[:, 2]
    global_rate = float(np.sum(k * counts) / max(np.sum(t * counts), 1.0))
    pi = float(initial_pi) if initial_pi is not None else float(np.clip(np.sum(counts[k > 0]) / np.sum(counts), 0.01, 0.95))
    p_old = (
        float(user_p_old) if user_p_old is not None else
        float(initial_p_old) if initial_p_old is not None else
        max(1e-6, min(global_rate * 0.15, 0.02))
    )
    p_new = (
        float(user_p_new) if user_p_new is not None else
        float(initial_p_new) if initial_p_new is not None else
        min(max(global_rate * 3.0, 0.01), 0.95)
    )
    if p_new <= p_old:
        p_new = min(0.95, p_old + 0.01)

    converged = False
    for iteration in range(1, max_iter + 1):
        old_pi, old_p_old, old_p_new = pi, p_old, p_new
        post = posterior_new_array(k, t, pi, p_old, p_new)
        new_weight = counts * post
        old_weight = counts * (1.0 - post)
        if user_p_old is None:
            p_old = float(np.sum(old_weight * k) / max(np.sum(old_weight * t), 1.0))
        if user_p_new is None:
            p_new = float(np.sum(new_weight * k) / max(np.sum(new_weight * t), 1.0))
        pi = float(np.sum(new_weight) / max(np.sum(counts), 1.0))
        p_old = float(np.clip(p_old, 1e-7, 0.5))
        p_new = float(np.clip(p_new, p_old + 1e-7, 0.999))
        pi = float(np.clip(pi, 1e-6, 0.999999))
        delta = max(abs(pi - old_pi), abs(p_old - old_p_old), abs(p_new - old_p_new))
        if delta < tol:
            converged = True
            break
    return {
        "pi": pi,
        "p_old": p_old,
        "p_new": p_new,
        "global_tc_rate": global_rate,
        "iterations": iteration,
        "converged": converged,
    }


def aggregate_corrected_counts(args, params, official_grid=None):
    outdir = Path(args.outdir)
    aggregate_parts = []
    molecule_handle = None
    if args.write_corrected_molecules:
        molecule_handle = gzip.open(outdir / f"{args.sample}.corrected_molecule_calls.tsv.gz", "wt", encoding="utf-8")
        molecule_handle.write(
            "barcode\tx\ty\tbin_id\tgene_id\tgene_name\tumi\ttc_events\tt_covered\tassignment\tposterior_new\n"
        )
    usecols = ["barcode", "x", "y", "gene_id", "gene_name", "umi", args.tc_column, args.t_column, "assignment"]
    for chunk in pd.read_csv(args.molecule_calls, sep="\t", usecols=usecols, chunksize=args.chunk_size):
        chunk[args.tc_column] = chunk[args.tc_column].fillna(0).astype(int)
        chunk[args.t_column] = chunk[args.t_column].fillna(0).astype(int)
        posterior = posterior_new_array(
            chunk[args.tc_column].to_numpy(),
            chunk[args.t_column].to_numpy(),
            params["pi"],
            params["p_old"],
            params["p_new"],
        )
        chunk["posterior_new"] = posterior
        chunk["posterior_variance"] = posterior * (1.0 - posterior)
        chunk = assign_bins(chunk, args, official_grid)
        if chunk.empty:
            continue
        grouped = chunk.groupby(["bin_id", "bin_x", "bin_y", "gene_id", "gene_name"], sort=False).agg(
            total=("posterior_new", "size"),
            expected_new=("posterior_new", "sum"),
            posterior_variance=("posterior_variance", "sum"),
            tc_events=(args.tc_column, "sum"),
            t_covered=(args.t_column, "sum"),
        )
        aggregate_parts.append(grouped.reset_index())
        if len(aggregate_parts) >= 10:
            merged = pd.concat(aggregate_parts, ignore_index=True).groupby(
                ["bin_id", "bin_x", "bin_y", "gene_id", "gene_name"], sort=False
            ).agg(
                total=("total", "sum"), expected_new=("expected_new", "sum"),
                posterior_variance=("posterior_variance", "sum"),
                tc_events=("tc_events", "sum"), t_covered=("t_covered", "sum"),
            )
            aggregate_parts = [merged.reset_index()]
        if molecule_handle is not None:
            cols = ["barcode", "x", "y", "bin_id", "gene_id", "gene_name", "umi", "tc_events", "t_covered", "assignment", "posterior_new"]
            chunk[cols].to_csv(molecule_handle, sep="\t", index=False, header=False, float_format="%.6g")
    if molecule_handle is not None:
        molecule_handle.close()
    merged = pd.concat(aggregate_parts, ignore_index=True).groupby(
        ["bin_id", "bin_x", "bin_y", "gene_id", "gene_name"], sort=False
    ).agg(
        total=("total", "sum"), expected_new=("expected_new", "sum"),
        posterior_variance=("posterior_variance", "sum"),
        tc_events=("tc_events", "sum"), t_covered=("t_covered", "sum"),
    )
    output = outdir / f"{args.sample}.corrected_bin_gene_counts.tsv.gz"
    columns = [
        "bin_id", "bin_x", "bin_y", "gene_id", "gene_name", "total", "expected_new",
        "expected_old", "posterior_variance", "ntr", "ntr_se", "ntr_ci_lower",
        "ntr_ci_upper", "tc_events", "t_covered",
    ]
    n_rows = 0
    pigz = shutil.which("pigz")
    with output.open("wb") as raw_handle:
        compressor = None
        if pigz:
            compressor = subprocess.Popen([pigz, "-1", "-p", "4", "-c"], stdin=subprocess.PIPE, stdout=raw_handle)
            text_handle = io.TextIOWrapper(compressor.stdin, encoding="utf-8", newline="")
        else:
            text_handle = gzip.open(raw_handle, "wt", encoding="utf-8", compresslevel=1)
        try:
            text_handle.write("\t".join(columns) + "\n")
            for start in range(0, len(merged), 500_000):
                block = merged.iloc[start:start + 500_000].reset_index()
                total = block["total"]
                expected_new = block["expected_new"]
                variance = block["posterior_variance"]
                standard_error = np.sqrt(np.maximum(variance, 0.0)) / total.where(total > 0)
                ntr = expected_new / total.where(total > 0)
                block["expected_old"] = total - expected_new
                block["ntr"] = ntr
                block["ntr_se"] = standard_error
                block["ntr_ci_lower"] = (ntr - 1.96 * standard_error).clip(lower=0.0)
                block["ntr_ci_upper"] = (ntr + 1.96 * standard_error).clip(upper=1.0)
                block.to_csv(text_handle, sep="\t", index=False, header=False, columns=columns, float_format="%.10g")
                n_rows += len(block)
            if compressor is not None:
                text_handle.close()
                compressor.wait()
                if compressor.returncode != 0:
                    raise RuntimeError("pigz failed while writing corrected bin-gene counts")
        finally:
            if compressor is not None and not text_handle.closed:
                text_handle.close()
            if compressor is None:
                text_handle.close()
    if args.output_pickle:
        merged.reset_index().to_pickle(args.output_pickle)
    return n_rows


def aggregate_matrix_counts(args, params, assignment_posteriors):
    matrix_dir = Path(args.matrix_dir)
    feature_meta = pd.read_csv(matrix_dir / "features.tsv", sep="\t", dtype={"gene_id": str})
    features = feature_meta["gene_id"].astype(str).tolist()
    gene_names = dict(zip(feature_meta["gene_id"].astype(str), feature_meta["gene_name"].astype(str)))
    barcode_meta = pd.read_csv(matrix_dir / "barcodes.tsv", sep="\t", dtype={"barcode": str})
    barcodes = barcode_meta["barcode"].astype(str).tolist()
    barcode_meta = barcode_meta.set_index("barcode", drop=False)
    barcode_meta["bin_x"] = (barcode_meta["x"] // args.bin_size).astype(int) * args.bin_size + args.bin_size / 2.0
    barcode_meta["bin_y"] = (barcode_meta["y"] // args.bin_size).astype(int) * args.bin_size + args.bin_size / 2.0
    barcode_meta["bin_id"] = barcode_meta["bin_x"].astype(int).astype(str) + "_" + barcode_meta["bin_y"].astype(int).astype(str)
    bin_ids = sorted(barcode_meta["bin_id"].unique())
    bin_index = {bin_id: idx for idx, bin_id in enumerate(bin_ids)}
    target_bins = barcode_meta["bin_id"].map(bin_index).to_numpy(dtype=int)
    barcode_to_bin = sparse.coo_matrix(
        (
            np.ones(len(barcodes), dtype=np.float32),
            (np.arange(len(barcodes)), target_bins),
        ),
        shape=(len(barcodes), len(bin_ids)),
        dtype=np.float32,
    ).tocsr()
    old = (mmread(matrix_dir / "old_umi.mtx").tocsr() @ barcode_to_bin).T.tocsr()
    new = (mmread(matrix_dir / "new_umi.mtx").tocsr() @ barcode_to_bin).T.tocsr()
    total = (old + new).tocsr()
    p_old_hard = assignment_posteriors["posterior_mean_by_hard_assignment"]["old"]
    p_new_hard = assignment_posteriors["posterior_mean_by_hard_assignment"]["new"]
    corrected = (old.multiply(p_old_hard) + new.multiply(p_new_hard)).tocsr()
    rows = []
    for bin_idx, bin_id in enumerate(bin_ids):
        start, end = total.indptr[bin_idx], total.indptr[bin_idx + 1]
        if start == end:
            continue
        genes = total.indices[start:end]
        totals = total.data[start:end]
        corr_start, corr_end = corrected.indptr[bin_idx], corrected.indptr[bin_idx + 1]
        corr_genes = corrected.indices[corr_start:corr_end]
        corr_values = corrected.data[corr_start:corr_end]
        if not np.array_equal(genes, corr_genes):
            corr_lookup = dict(zip(corr_genes.tolist(), corr_values.tolist()))
            expected_values = [corr_lookup.get(int(gene_idx), 0.0) for gene_idx in genes]
        else:
            expected_values = corr_values
        bin_id = bin_ids[int(bin_idx)]
        bin_x, bin_y = bin_id.split("_")
        for gene_idx, total_value, expected_new in zip(genes, totals, expected_values):
            gene_id = features[int(gene_idx)]
            rows.append(
                {
                    "bin_id": bin_id,
                    "bin_x": float(bin_x),
                    "bin_y": float(bin_y),
                    "gene_id": gene_id,
                    "gene_name": gene_names.get(gene_id, gene_id),
                    "total": float(total_value),
                    "expected_new": float(expected_new),
                    "expected_old": float(total_value) - float(expected_new),
                    "ntr": float(expected_new) / float(total_value) if total_value else 0.0,
                }
            )
    frame = pd.DataFrame(rows)
    write_corrected_frame(frame, args)
    return frame


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    official_grid = load_official_grid(args.official_bin_positions) if args.official_bin_positions else None
    pair_counts, qc = read_pair_counts(
        args.molecule_calls, args.chunk_size, args.max_estimate_rows, args.tc_column, args.t_column
    )
    params = fit_mixture(
        pair_counts,
        user_p_old=args.background_conversion_rate,
        user_p_new=args.marked_conversion_rate,
        max_iter=args.max_em_iter,
        tol=args.tol,
        initial_pi=args.initial_pi,
        initial_p_old=args.initial_p_old,
        initial_p_new=args.initial_p_new,
    )
    assignment_posteriors = estimate_assignment_posteriors(
        args.molecule_calls, args.chunk_size, params, args.max_estimate_rows, args.tc_column, args.t_column
    )
    if args.matrix_dir:
        frame = aggregate_matrix_counts(args, params, assignment_posteriors)
        n_bin_gene_rows = len(frame)
        aggregation_mode = "matrix_bin_gene_with_sampled_posterior_calibration"
    else:
        n_bin_gene_rows = aggregate_corrected_counts(args, params, official_grid)
        aggregation_mode = "full_molecule_posterior"
    summary = {
        "sample": args.sample,
        "bin_size": args.bin_size,
        "bin_size_pixels": int(args.bin_size_pixels or args.bin_size),
        "official_bin_positions": args.official_bin_positions,
        "model": "two_component_binomial_mixture",
        "aggregation_mode": aggregation_mode,
        "tc_column": args.tc_column,
        "t_column": args.t_column,
        "parameters": params,
        "assignment_posteriors": assignment_posteriors,
        "qc": qc,
        "n_bin_gene_rows": int(n_bin_gene_rows),
        "corrected_bin_gene_counts": None if args.skip_tsv else str(outdir / f"{args.sample}.corrected_bin_gene_counts.tsv.gz"),
        "corrected_bin_gene_pickle": args.output_pickle,
        "corrected_molecule_calls": str(outdir / f"{args.sample}.corrected_molecule_calls.tsv.gz")
        if args.write_corrected_molecules
        else None,
    }
    (outdir / "new_fraction_correction_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
