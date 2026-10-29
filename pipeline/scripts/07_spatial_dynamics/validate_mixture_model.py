#!/usr/bin/env python3
"""Validate Spatial SLAM mixture fits across counting, BQ and background choices."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from spatial_slam_correct_new_fraction import fit_mixture
from spatial_slam_correct_new_fraction import posterior_new_array


INK = "#25313A"
BLUE = "#3C5488"
CYAN = "#4DBBD5"
RED = "#E64B35"


def parse_methods(frame_columns, thresholds):
    methods = {
        "raw_sum_bq0_unmasked": ("raw_tc_events", "raw_t_covered"),
        "unique_site_bq0_unmasked": ("unique_tc_events_bq0", "unique_t_covered_bq0"),
    }
    for threshold in thresholds:
        tc = f"tc_events_bq{threshold}"
        covered = f"t_covered_bq{threshold}"
        if tc in frame_columns and covered in frame_columns:
            methods[f"unique_site_bq{threshold}_snp_masked"] = (tc, covered)
    return methods


def collect_pair_counts(path, thresholds, chunk_size):
    header = pd.read_csv(path, sep="\t", nrows=0)
    methods = parse_methods(set(header.columns), thresholds)
    usecols = sorted({column for pair in methods.values() for column in pair})
    counts = {method: Counter() for method in methods}
    totals = {method: Counter() for method in methods}
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=chunk_size):
        for method, (tc_col, t_col) in methods.items():
            tc = chunk[tc_col].fillna(0).astype(int)
            covered = chunk[t_col].fillna(0).astype(int)
            totals[method]["molecules"] += len(chunk)
            totals[method]["informative"] += int((covered > 0).sum())
            totals[method]["tc"] += int(tc.sum())
            totals[method]["t"] += int(covered.sum())
            grouped = pd.DataFrame({"k": tc[covered > 0], "t": covered[covered > 0]}).groupby(["k", "t"]).size()
            for (k, t), value in grouped.items():
                counts[method][(int(k), int(t))] += int(value)
    return methods, counts, totals


def parameter_row(method, params, totals, fixed_background=None):
    return {
        "method": method,
        "fixed_background": fixed_background,
        "molecules": totals["molecules"],
        "informative_molecules": totals["informative"],
        "tc_events": totals["tc"],
        "t_covered": totals["t"],
        "global_tc_rate": totals["tc"] / totals["t"] if totals["t"] else 0.0,
        "pi": params["pi"],
        "p_old": params["p_old"],
        "p_new": params["p_new"],
        "iterations": params["iterations"],
        "converged": params["converged"],
    }


def mismatch_background(path, min_base_quality):
    if not path:
        return None
    frame = pd.read_csv(path, sep="\t")
    selected = frame[
        (frame["scope"] == "global")
        & (frame["min_base_quality"] == min_base_quality)
        & (frame["channel"] != "T>C")
    ]
    values = selected.loc[selected["reference_base_evaluations"] > 0, "rate"].to_numpy(dtype=float)
    return float(np.median(values)) if len(values) else None


def initialization_runs(pair_counts):
    baseline = fit_mixture(pair_counts)
    starts = [("default", None, None, None)]
    starts.extend((f"pi_{value}", value, baseline["p_old"], baseline["p_new"]) for value in (0.05, 0.1, 0.2, 0.3, 0.5))
    starts.extend((f"pold_{value}", baseline["pi"], value, baseline["p_new"]) for value in (0.0001, 0.0005, 0.001, 0.002))
    starts.extend((f"pnew_{value}", baseline["pi"], baseline["p_old"], value) for value in (0.03, 0.05, 0.08, 0.12))
    rows = []
    for label, pi0, pold0, pnew0 in starts:
        fit = fit_mixture(pair_counts, initial_pi=pi0, initial_p_old=pold0, initial_p_new=pnew0)
        rows.append({
            "initialization": label,
            "initial_pi": pi0,
            "initial_p_old": pold0,
            "initial_p_new": pnew0,
            **fit,
            "delta_pi_vs_default": fit["pi"] - baseline["pi"],
            "delta_p_old_vs_default": fit["p_old"] - baseline["p_old"],
            "delta_p_new_vs_default": fit["p_new"] - baseline["p_new"],
        })
    return pd.DataFrame(rows)


def poisson_bootstrap(pair_counts, replicates, seed):
    rng = np.random.default_rng(seed)
    pairs = list(pair_counts)
    weights = np.array([pair_counts[pair] for pair in pairs], dtype=np.int64)
    rows = []
    for replicate in range(replicates):
        sampled = rng.poisson(weights)
        boot = {pair: int(value) for pair, value in zip(pairs, sampled) if value > 0}
        fit = fit_mixture(boot)
        rows.append({"replicate": replicate + 1, **fit})
    return pd.DataFrame(rows)


def save_validation_figure(initialization, bootstrap, output):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True, facecolor="white")
    x = np.arange(len(initialization))
    axes[0].plot(x, initialization["pi"], marker="o", color=BLUE, label="pi")
    axes[0].plot(x, initialization["p_new"], marker="o", color=RED, label="p_new")
    axes[0].plot(x, initialization["p_old"], marker="o", color=CYAN, label="p_old")
    axes[0].set_xticks(x, initialization["initialization"], rotation=70, ha="right", fontsize=7)
    axes[0].set_title("A  EM initialization sensitivity", color=INK, fontweight="bold")
    axes[0].set_ylabel("Fitted parameter")
    axes[0].legend(frameon=False, fontsize=8)
    bootstrap[["pi", "p_old", "p_new"]].plot.box(ax=axes[1], color={"boxes": BLUE, "medians": RED})
    axes[1].set_title("B  Poisson frequency-table bootstrap", color=INK, fontweight="bold")
    axes[1].set_ylabel("Bootstrap estimate")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.55)
    for suffix in ("png", "pdf"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def counting_comparison(base_quality, pair_counts, registry, default_method, outdir):
    methods = ["raw_sum_bq0_unmasked", default_method]
    labels = {"raw_sum_bq0_unmasked": "old_raw_summed_sites", default_method: "reviewed_unique_site_bq20_snp_masked"}
    comparison = base_quality[base_quality.method.isin(methods)].copy()
    comparison["counting_method"] = comparison.method.map(labels)
    comparison["comparison_scope"] = comparison.method.map({
        "raw_sum_bq0_unmasked": "reviewed raw-sum proxy; N-base denominator differs",
        default_method: "reviewed unique genomic-site implementation",
    })
    comparison["global_expected_new_fraction"] = comparison["pi"]
    comparison.to_csv(outdir / "molecule_counting_method_comparison.tsv", sep="\t", index=False)
    posterior_rows = []
    for method in methods:
        params = registry[method]["parameters"]
        for (k, t), count in pair_counts[method].items():
            value = float(posterior_new_array(np.array([k]), np.array([t]), params["pi"], params["p_old"], params["p_new"])[0])
            posterior_rows.append({"method": labels[method], "tc_events": k, "t_covered": t,
                                   "posterior_new": value, "molecules": count})
    posterior = pd.DataFrame(posterior_rows)
    posterior.to_csv(outdir / "molecule_posterior_distribution.tsv.gz", sep="\t", index=False, compression="gzip")
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), constrained_layout=True, facecolor="white")
    x = np.arange(len(comparison))
    width = .22
    for offset, column, color in ((-width, "p_old", CYAN), (0, "pi", BLUE), (width, "p_new", RED)):
        axes[0].bar(x + offset, comparison[column], width, label=column, color=color)
    axes[0].set_xticks(x, [labels[value] for value in comparison.method], rotation=15, ha="right", fontsize=8)
    axes[0].set_ylabel("Fitted parameter"); axes[0].legend(frameon=False)
    for method, color in zip(methods, ("#8491B4", RED)):
        sub = posterior[posterior.method == labels[method]]
        axes[1].hist(sub.posterior_new, bins=np.linspace(0, 1, 41), weights=sub.molecules,
                     histtype="step", linewidth=1.8, color=color, label=labels[method])
    axes[1].set_yscale("log"); axes[1].set_xlabel("Molecule posterior P(new)"); axes[1].set_ylabel("Molecules")
    axes[1].legend(frameon=False, fontsize=7)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E5E7EB", linewidth=.5)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"molecule_counting_method_comparison.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)


def substitution_qc_figures(spectrum_path, default_bq, outdir):
    if not spectrum_path or not Path(spectrum_path).is_file():
        return
    metadata = Path(spectrum_path).parent
    spectrum = pd.read_csv(spectrum_path, sep="\t")
    global_spectrum = spectrum[(spectrum.scope == "global") & (spectrum.min_base_quality == default_bq)]
    fig, ax = plt.subplots(figsize=(9.5, 4.8), facecolor="white")
    colors = [RED if value == "T>C" else BLUE for value in global_spectrum.channel]
    ax.bar(global_spectrum.channel, global_spectrum.rate, color=colors)
    ax.set(ylabel="Mismatch rate", title=f"Transcript-oriented substitution spectrum (BQ≥{default_bq})")
    ax.tick_params(axis="x", rotation=45); ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(outdir / f"spatial_substitution_spectrum.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)

    plots = [
        ("tc_by_read_position.tsv", "read_position_1based", "tc_rate", "tc_by_read_position", "Read position", "T>C rate"),
        ("tc_by_base_quality.tsv", "base_quality", "tc_rate", "tc_by_base_quality", "Base quality", "T>C rate"),
        ("lane_tc_background.tsv", "lane", "tc_rate", "lane_tc_background", "Lane", "T>C rate"),
    ]
    for filename, xcol, ycol, stem, xlabel, ylabel in plots:
        path = metadata / filename
        if not path.is_file():
            continue
        frame = pd.read_csv(path, sep="\t")
        if "min_base_quality" in frame:
            frame = frame[frame.min_base_quality == default_bq]
        fig, ax = plt.subplots(figsize=(8.8, 4.8), facecolor="white")
        if xcol == "lane":
            ax.bar(frame[xcol].astype(str), frame[ycol], color=BLUE)
        else:
            ax.plot(frame[xcol], frame[ycol], color=BLUE, linewidth=1.5)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=stem.replace("_", " ").title())
        ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", color="#E5E7EB", linewidth=.5)
        for suffix in ("png", "pdf"):
            fig.savefig(outdir / f"{stem}.{suffix}", dpi=450, bbox_inches="tight")
        plt.close(fig)
    artifact_rows = []
    position = pd.read_csv(metadata / "tc_by_read_position.tsv", sep="\t")
    position = position[position.min_base_quality == default_bq].copy()
    if not position.empty:
        maximum = int(position.read_position_1based.max())
        position["region"] = np.where(
            (position.read_position_1based <= 10) | (position.read_position_1based > maximum - 10), "read_end_10bp", "read_middle"
        )
        for region, sub in position.groupby("region"):
            artifact_rows.append({"diagnostic": "read_position", "group": region,
                                  "evaluations": int(sub.t_evaluations.sum()), "tc_events": int(sub.tc_events.sum()),
                                  "tc_rate": float(sub.tc_events.sum() / max(sub.t_evaluations.sum(), 1))})
    quality = pd.read_csv(metadata / "tc_by_base_quality.tsv", sep="\t")
    quality["group"] = np.where(quality.base_quality < default_bq, f"BQ<{default_bq}", f"BQ>={default_bq}")
    for group, sub in quality.groupby("group"):
        artifact_rows.append({"diagnostic": "base_quality", "group": group,
                              "evaluations": int(sub.t_evaluations.sum()), "tc_events": int(sub.tc_events.sum()),
                              "tc_rate": float(sub.tc_events.sum() / max(sub.t_evaluations.sum(), 1))})
    lanes = pd.read_csv(metadata / "lane_tc_background.tsv", sep="\t")
    lanes = lanes[lanes.min_base_quality == default_bq]
    for row in lanes.itertuples(index=False):
        artifact_rows.append({"diagnostic": "lane", "group": row.lane, "evaluations": int(row.t_evaluations),
                              "tc_events": int(row.tc_events), "tc_rate": float(row.tc_rate)})
    pd.DataFrame(artifact_rows).to_csv(outdir / "substitution_artifact_summary.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule-calls", required=True)
    parser.add_argument("--substitution-spectrum")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--base-quality-thresholds", default="15,20,25,30")
    parser.add_argument("--default-base-quality", type=int, default=20)
    parser.add_argument("--grand-slam-p-old", type=float, default=0.001303)
    parser.add_argument("--bootstrap-replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=2_000_000)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = [int(value) for value in args.base_quality_thresholds.split(",") if value.strip()]
    methods, pair_counts, totals = collect_pair_counts(args.molecule_calls, thresholds, args.chunk_size)

    parameter_rows = []
    registry = {}
    for method, columns in methods.items():
        fit = fit_mixture(pair_counts[method])
        parameter_rows.append(parameter_row(method, fit, totals[method]))
        registry[method] = {"tc_column": columns[0], "t_column": columns[1], "parameters": fit}
    base_quality = pd.DataFrame(parameter_rows)
    base_quality.to_csv(outdir / "base_quality_sensitivity.tsv", sep="\t", index=False)

    default_method = f"unique_site_bq{args.default_base_quality}_snp_masked"
    if default_method not in methods:
        raise ValueError(f"Default method missing from molecule table: {default_method}")
    initialization = initialization_runs(pair_counts[default_method])
    initialization.to_csv(outdir / "em_initialization_sensitivity.tsv", sep="\t", index=False)
    bootstrap = poisson_bootstrap(pair_counts[default_method], args.bootstrap_replicates, args.seed)
    bootstrap.to_csv(outdir / "em_parameter_bootstrap.tsv.gz", sep="\t", index=False, compression="gzip")
    uncertainty_rows = []
    for parameter in ("pi", "p_old", "p_new"):
        values = bootstrap[parameter].to_numpy(dtype=float)
        uncertainty_rows.append({
            "parameter": parameter,
            "estimate": registry[default_method]["parameters"][parameter],
            "bootstrap_mean": float(np.mean(values)),
            "bootstrap_sd": float(np.std(values, ddof=1)),
            "ci_lower_2.5pct": float(np.quantile(values, 0.025)),
            "ci_upper_97.5pct": float(np.quantile(values, 0.975)),
            "replicates": len(values),
            "bootstrap_type": "poisson_weighted_frequency_table",
        })
    uncertainty = pd.DataFrame(uncertainty_rows)
    uncertainty.to_csv(outdir / "em_parameter_uncertainty.tsv", sep="\t", index=False)

    mismatch_p_old = mismatch_background(args.substitution_spectrum, args.default_base_quality)
    fixed_rows = [parameter_row("free_p_old", registry[default_method]["parameters"], totals[default_method], "free")]
    fixed_specs = [("grand_slam_p_old", args.grand_slam_p_old)]
    if mismatch_p_old is not None:
        fixed_specs.append(("mismatch_median_non_tc", mismatch_p_old))
    for label, p_old in fixed_specs:
        fit = fit_mixture(pair_counts[default_method], user_p_old=p_old)
        method = f"fixed_{label}"
        fixed_rows.append(parameter_row(method, fit, totals[default_method], p_old))
        registry[method] = {**registry[default_method], "parameters": fit, "fixed_background": p_old}
    fixed = pd.DataFrame(fixed_rows)
    fixed.to_csv(outdir / "fixed_background_sensitivity.tsv", sep="\t", index=False)
    save_validation_figure(initialization, bootstrap, outdir / "mixture_model_validation")
    counting_comparison(base_quality, pair_counts, registry, default_method, outdir)
    substitution_qc_figures(args.substitution_spectrum, args.default_base_quality, outdir)

    payload = {
        "default_method": default_method,
        "methods": registry,
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_seed": args.seed,
        "mismatch_background_definition": "median rate of the 11 non-T>C transcript-oriented substitution channels",
        "mismatch_background": mismatch_p_old,
    }
    (outdir / "model_validation_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
