#!/usr/bin/env python3
"""Propagate reviewed counting/model alternatives to bin, domain and candidate rankings."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from spatial_slam_correct_new_fraction import posterior_new_array


def safe_spearman(left, right):
    keep = np.isfinite(left) & np.isfinite(right)
    if keep.sum() < 3 or np.unique(left[keep]).size < 2 or np.unique(right[keep]).size < 2:
        return np.nan
    return float(spearmanr(left[keep], right[keep]).statistic)


def load_methods(path):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    methods = {}
    for name, record in payload["methods"].items():
        params = record["parameters"]
        methods[name] = {
            "tc": record["tc_column"], "t": record["t_column"],
            "pi": params["pi"], "p_old": params["p_old"], "p_new": params["p_new"],
        }
    return payload["default_method"], methods


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls", required=True); parser.add_argument("--model-summary", required=True)
    parser.add_argument("--domains", required=True); parser.add_argument("--positions", required=True)
    parser.add_argument("--domain-gene", required=True); parser.add_argument("--outdir", required=True)
    parser.add_argument("--bin-size-pixels", type=int, default=40); parser.add_argument("--chunk-size", type=int, default=1_000_000)
    parser.add_argument("--candidate-genes", type=int, default=500)
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    default_method, methods = load_methods(args.model_summary)
    domains = pd.read_csv(args.domains, sep="\t", dtype={"bin_id": str}).rename(columns={"cluster": "domain"})
    domains["domain"] = domains.domain.astype(str)
    bin_ids = domains.bin_id.astype(str).tolist(); bin_index = {value: index for index, value in enumerate(bin_ids)}
    bin_domains = domains.domain.to_numpy(str)
    positions = pd.read_csv(args.positions, header=None,
        names=["bin_id", "in_tissue", "row", "col", "pxl_row", "pxl_col"], dtype={"bin_id": str})
    positions = positions[positions.bin_id.isin(bin_index)].copy()
    min_col, max_col = int(positions.col.min()), int(positions.col.max())
    min_row, max_row = int(positions.row.min()), int(positions.row.max())
    grid = np.full((max_row - min_row + 1, max_col - min_col + 1), -1, dtype=np.int32)
    for row in positions.itertuples(index=False):
        grid[int(row.row) - min_row, int(row.col) - min_col] = bin_index[str(row.bin_id)]
    domain_gene = pd.read_csv(args.domain_gene, sep="\t", dtype={"gene_id": str, "domain": str})
    domain_gene["abs_dynamic"] = domain_gene.dynamic_score.abs()
    candidate_ids = set(domain_gene.nlargest(args.candidate_genes, "abs_dynamic").gene_id.astype(str))

    totals = {method: {"n": 0, "new": 0.0} for method in methods}
    bin_total = np.zeros(len(bin_ids), dtype=np.int64)
    bin_new = {method: np.zeros(len(bin_ids), dtype=np.float64) for method in methods}
    candidate_total = defaultdict(int)
    candidate_new = {method: defaultdict(float) for method in methods}
    usecols = sorted({"x", "y", "gene_id", *[item[axis] for item in methods.values() for axis in ("tc", "t")]})
    for chunk in pd.read_csv(args.molecule_calls, sep="\t", usecols=usecols, chunksize=args.chunk_size,
                             dtype={"gene_id": str}):
        cols = (chunk.x.to_numpy(int) // args.bin_size_pixels).astype(int)
        rows = (chunk.y.to_numpy(int) // args.bin_size_pixels).astype(int)
        mapped = np.full(len(chunk), -1, dtype=np.int32)
        inside = (cols >= min_col) & (cols <= max_col) & (rows >= min_row) & (rows <= max_row)
        mapped[inside] = grid[rows[inside] - min_row, cols[inside] - min_col]
        keep = mapped >= 0
        if not np.any(keep):
            continue
        chunk = chunk.loc[keep].reset_index(drop=True); mapped = mapped[keep]
        bin_total += np.bincount(mapped, minlength=len(bin_ids))
        candidate_mask = chunk.gene_id.isin(candidate_ids).to_numpy()
        candidate_frame = None
        candidate_domains = None
        candidate_genes_chunk = None
        if candidate_mask.any():
            candidate_domains = bin_domains[mapped[candidate_mask]]
            candidate_genes_chunk = chunk.loc[candidate_mask, "gene_id"].astype(str).to_numpy()
            candidate_frame = pd.DataFrame({"domain": candidate_domains, "gene_id": candidate_genes_chunk})
            for key, value in candidate_frame.groupby(["domain", "gene_id"], sort=False).size().items():
                candidate_total[key] += value
        cache = {}
        for method, spec in methods.items():
            pair = (spec["tc"], spec["t"])
            if pair not in cache:
                cache[pair] = (chunk[pair[0]].fillna(0).to_numpy(int), chunk[pair[1]].fillna(0).to_numpy(int))
            k, t = cache[pair]
            posterior = posterior_new_array(k, t, spec["pi"], spec["p_old"], spec["p_new"])
            totals[method]["n"] += len(posterior); totals[method]["new"] += float(posterior.sum())
            bin_new[method] += np.bincount(mapped, weights=posterior, minlength=len(bin_ids))
            if candidate_frame is not None:
                candidate_frame["posterior"] = posterior[candidate_mask]
                grouped = candidate_frame.groupby(["domain", "gene_id"], sort=False).posterior.sum()
                for key, value in grouped.items():
                    candidate_new[method][key] += float(value)

    formal_bin = bin_new[default_method] / np.maximum(bin_total, 1)
    formal_domain = {}
    for domain in sorted(set(bin_domains)):
        mask = bin_domains == domain
        formal_domain[domain] = float(bin_new[default_method][mask].sum() / max(bin_total[mask].sum(), 1))
    keys_by_gene = defaultdict(list)
    for key in candidate_total:
        keys_by_gene[key[1]].append(key)
    formal_gene_total = {gene: sum(candidate_total[item] for item in keys) for gene, keys in keys_by_gene.items()}
    formal_gene_new = {gene: sum(candidate_new[default_method][item] for item in keys) for gene, keys in keys_by_gene.items()}
    formal_scores = {}
    for key, count in candidate_total.items():
        domain, gene = key
        formal_scores[key] = (
            candidate_new[default_method][key] / count - formal_gene_new[gene] / formal_gene_total[gene]
        ) * np.log1p(count)

    summary_rows, bin_rows, domain_rows, candidate_rows = [], [], [], []
    for method in methods:
        ntr = bin_new[method] / np.maximum(bin_total, 1)
        for index, bin_id in enumerate(bin_ids):
            bin_rows.append({"method": method, "bin_id": bin_id, "domain": bin_domains[index],
                             "molecules": int(bin_total[index]), "expected_new": bin_new[method][index], "ntr": ntr[index]})
        domain_values = []
        for domain in sorted(set(bin_domains)):
            mask = bin_domains == domain
            value = float(bin_new[method][mask].sum() / max(bin_total[mask].sum(), 1))
            domain_values.append(value)
            domain_rows.append({"method": method, "domain": domain, "molecules": int(bin_total[mask].sum()), "ntr": value})
        scores = {}
        method_gene_new = {gene: sum(candidate_new[method][item] for item in keys) for gene, keys in keys_by_gene.items()}
        for key, count in candidate_total.items():
            domain, gene = key
            value = (
                candidate_new[method][key] / count - method_gene_new[gene] / formal_gene_total[gene]
            ) * np.log1p(count)
            scores[key] = value
            candidate_rows.append({"method": method, "domain": domain, "gene_id": gene,
                                   "molecules": count, "dynamic_score": value})
        summary_rows.append({
            "method": method, "pi": methods[method]["pi"], "p_old": methods[method]["p_old"],
            "p_new": methods[method]["p_new"], "mapped_molecules": totals[method]["n"],
            "mapped_bin_global_ntr": totals[method]["new"] / max(totals[method]["n"], 1),
            "bin_ntr_spearman_vs_default": safe_spearman(ntr, formal_bin),
            "domain_ntr_spearman_vs_default": safe_spearman(np.asarray(domain_values), np.asarray(list(formal_domain.values()))),
            "candidate_score_spearman_vs_default": safe_spearman(np.asarray([scores[k] for k in formal_scores]), np.asarray([formal_scores[k] for k in formal_scores])),
        })
        summary_rows[-1]["counting_scope"] = (
            "reviewed raw-sum proxy; N-base denominator differs" if method == "raw_sum_bq0_unmasked"
            else "reviewed unique genomic-site implementation"
        )
    summary_frame = pd.DataFrame(summary_rows)
    summary_frame.to_csv(outdir / "method_spatial_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(bin_rows).to_csv(outdir / "method_bin_ntr.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(domain_rows).to_csv(outdir / "method_domain_ntr.tsv", sep="\t", index=False)
    pd.DataFrame(candidate_rows).to_csv(outdir / "method_candidate_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    comparison_path = outdir / "molecule_counting_method_comparison.tsv"
    if comparison_path.is_file():
        comparison = pd.read_csv(comparison_path, sep="\t")
        sensitivity = summary_frame.set_index("method")
        comparison = comparison.merge(
            sensitivity[["bin_ntr_spearman_vs_default", "domain_ntr_spearman_vs_default", "candidate_score_spearman_vs_default"]],
            left_on="method", right_index=True, how="left",
        )
        comparison.to_csv(comparison_path, sep="\t", index=False)
    (outdir / "method_spatial_sensitivity_summary.json").write_text(json.dumps({
        "default_method": default_method, "methods": list(methods), "candidate_gene_pool": len(candidate_ids),
        "mapped_molecules": int(bin_total.sum()), "statistics_use_raw_unsmoothed_counts": True,
        "raw_method_scope": "reviewed raw-sum proxy; N-base denominator differs",
    }, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
