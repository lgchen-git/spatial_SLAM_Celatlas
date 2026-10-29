#!/usr/bin/env python3
"""Summarize and plot tissue-only slamdunk and GRAND-SLAM results."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter


NPG_GENERAL = LinearSegmentedColormap.from_list(
    "npg_general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"]
)
INK = "#25313A"
MUTED = "#64727B"
BLUE = "#3C5488"
CYAN = "#4DBBD5"
GREEN = "#00A087"
RED = "#E64B35"
INTERPRETATION_MODE = "standard"


def table_meta(paths):
    result = []
    for path in paths:
        try:
            frame = pd.read_csv(path, sep="\t", comment="#", nrows=5, compression="infer", dtype={"Chromosome": str})
            result.append({"file": str(path), "columns": list(frame.columns), "preview_rows": int(len(frame))})
        except (OSError, UnicodeDecodeError, pd.errors.ParserError):
            result.append({"file": str(path), "columns": [], "preview_rows": 0})
    return result


def find_suffix(columns, suffix):
    return next((column for column in columns if str(column).endswith(suffix)), None)


def read_grand(prefix):
    candidates = [Path(f"{prefix}.tsv.gz"), Path(f"{prefix}.tsv"), Path(f"{prefix}.ext.tsv")]
    path = next((candidate for candidate in candidates if candidate.is_file()), None)
    if path is None:
        return None, pd.DataFrame(), {}
    frame = pd.read_csv(path, sep="\t", compression="infer")
    columns = list(frame.columns)
    mapping = {
        "readcount": find_suffix(columns, " Readcount"),
        "lower": find_suffix(columns, " 0.05 quantile"),
        "mean": find_suffix(columns, " Mean"),
        "map": find_suffix(columns, " MAP"),
        "upper": find_suffix(columns, " 0.95 quantile"),
    }
    return path, frame, mapping


def read_grand_aux(prefix, suffix):
    path = Path(f"{prefix}.{suffix}.tsv")
    if not path.is_file():
        return path, pd.DataFrame()
    frame = pd.read_csv(path, sep="\t")
    frame.columns = [str(column).lstrip(".") for column in frame.columns]
    return path, frame


def read_slamdunk(tcounts):
    frames = []
    for path in tcounts:
        try:
            frame = pd.read_csv(path, sep="\t", comment="#", dtype={"Chromosome": str}, low_memory=False)
            frame["source_file"] = str(path)
            frames.append(frame)
        except (OSError, UnicodeDecodeError, pd.errors.ParserError):
            continue
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def numeric(frame, columns):
    result = frame.copy()
    for column in columns:
        if column and column in result:
            result[column] = pd.to_numeric(result[column], errors="coerce")
    return result


def interval_errors(means, lower, upper):
    """Return non-negative plotting errors while preserving raw estimates elsewhere.

    Some GRAND-SLAM edge-case fits emit a point estimate outside the reported
    quantile endpoints.  Matplotlib rejects the resulting negative error bar.
    For display only, use the envelope of all three values and report which
    rows needed that adjustment; callers retain the unmodified raw values in
    tables, annotations, and JSON provenance.
    """
    means = np.asarray(means, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    plot_lower = np.minimum.reduce([means, lower, upper])
    plot_upper = np.maximum.reduce([means, lower, upper])
    adjusted = (plot_lower != lower) | (plot_upper != upper)
    return means - plot_lower, plot_upper - means, plot_lower, plot_upper, adjusted


def style_axes(axes):
    for axis in np.ravel(axes):
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.55, zorder=0)
        axis.tick_params(labelsize=8, colors=MUTED)
        axis.spines[["top", "right"]].set_visible(False)


def save_figure(fig, output_dir, stem, dpi):
    output_dir.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"{stem}.{suffix}", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def gene_labels(frame):
    symbol = frame["Symbol"] if "Symbol" in frame else pd.Series(index=frame.index, dtype=object)
    gene = frame["Gene"] if "Gene" in frame else pd.Series(frame.index.astype(str), index=frame.index)
    return symbol.replace("", np.nan).fillna(gene).astype(str)


def plot_slamdunk(slam, grand, output_dir, dpi, min_reads, top_genes):
    needed = {"ReadCount", "ConversionRate"}
    if slam.empty or not needed.issubset(slam.columns):
        return None, {}
    slam = numeric(slam, ["ReadCount", "ConversionRate", "CoverageOnTs", "ConversionsOnTs", "TcReadCount"])
    valid = slam[(slam["ReadCount"] > 0) & slam["ConversionRate"].notna()].copy()
    if valid.empty:
        return None, {"min_reads": min_reads, "n_regions_with_reads": 0, "n_regions_ge_min_reads": 0}
    supported = valid[valid["ReadCount"] >= min_reads].copy()
    if supported.empty:
        return None, {"min_reads": min_reads, "n_regions_with_reads": int(len(valid)), "n_regions_ge_min_reads": 0}
    display_ceiling = float(supported["ConversionRate"].quantile(0.98))
    display_ceiling = max(display_ceiling, 1e-4)
    visible = supported[supported["ConversionRate"] <= display_ceiling]

    symbols = {}
    if not grand.empty and {"Gene", "Symbol"}.issubset(grand.columns):
        symbols = grand.dropna(subset=["Gene"]).drop_duplicates("Gene").set_index("Gene")["Symbol"].to_dict()
    coverage_col = "CoverageOnTs" if "CoverageOnTs" in supported else "ReadCount"
    representatives = supported.sort_values(coverage_col, ascending=False)
    if "Name" in representatives:
        representatives = representatives.drop_duplicates("Name")
        representatives["label"] = representatives["Name"].map(symbols).fillna(representatives["Name"])
    else:
        representatives["label"] = representatives.index.astype(str)
    top = representatives.nlargest(top_genes, "ReadCount").sort_values("ConversionRate")

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.1), facecolor="white", constrained_layout=True)
    ax = axes[0]
    bins = np.linspace(0, display_ceiling, 31)
    ax.hist(np.minimum(supported["ConversionRate"], display_ceiling), bins=bins, color=CYAN, edgecolor="white", linewidth=0.35)
    median = float(supported["ConversionRate"].median())
    ax.axvline(median, color=RED, linestyle="--", linewidth=1.4, label=f"median {median:.2%}")
    ax.legend(frameon=False, fontsize=8)
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_title("A  Conversion distribution", fontweight="bold", color=INK)
    ax.set_xlabel(f"T>C conversion rate (display capped at p98 = {display_ceiling:.2%})")
    ax.set_ylabel("UTR regions")

    ax = axes[1]
    hexbin = ax.hexbin(
        np.log10(visible["ReadCount"] + 1), visible["ConversionRate"], gridsize=48,
        mincnt=1, bins="log", cmap=NPG_GENERAL, linewidths=0, rasterized=True,
    )
    colorbar = fig.colorbar(hexbin, ax=ax, fraction=0.045, pad=0.02)
    colorbar.set_label("UTR-region density", fontsize=8)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_title("B  Coverage versus conversion", fontweight="bold", color=INK)
    ax.set_xlabel("log10(read count + 1)")
    ax.set_ylabel("T>C conversion rate")

    ax = axes[2]
    y = np.arange(len(top))
    ax.hlines(y, 0, top["ConversionRate"], color="#B7DDE5", linewidth=2)
    ax.scatter(top["ConversionRate"], y, color=BLUE, s=28, zorder=3)
    ax.set_yticks(y, top["label"], fontsize=8)
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlim(left=0)
    ax.set_title(f"C  Top {len(top)} high-read genes", fontweight="bold", color=INK)
    ax.set_xlabel("T>C conversion rate")
    style_axes(axes)
    fig.suptitle(
        f"slamdunk tissue-only conversion QC · support ≥ {min_reads:,} reads",
        fontsize=13, fontweight="bold", color=INK,
    )
    stem = "02_slamdunk_conversion_qc"
    save_figure(fig, output_dir, stem, dpi)
    return stem, {
        "min_reads": min_reads,
        "n_regions_with_reads": int(len(valid)),
        "n_regions_ge_min_reads": int(len(supported)),
        "median_supported_conversion_rate": median,
        "display_conversion_ceiling_p98": display_ceiling,
    }


def plot_grand_genes(grand, mapping, output_dir, dpi, min_reads, top_genes):
    required = [mapping.get(key) for key in ("readcount", "lower", "mean", "upper")]
    if grand.empty or not all(required):
        return None, {}
    read_col, lower_col, mean_col, upper_col = required
    work = numeric(grand, required).dropna(subset=required).copy()
    supported = work[work[read_col] >= min_reads].copy()
    if supported.empty:
        return None, {"min_reads": min_reads, "n_genes_ge_min_reads": 0}
    supported["interval_width"] = (supported[upper_col] - supported[lower_col]).abs()
    top = supported.nlargest(top_genes, read_col).sort_values(mean_col)

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.1), facecolor="white", constrained_layout=True)
    ax = axes[0]
    values = supported[mean_col].clip(0, 1)
    ax.hist(values, bins=np.linspace(0, 1, 31), color=GREEN, alpha=0.84, edgecolor="white", linewidth=0.35)
    median = float(values.median())
    ax.axvline(median, color=RED, linewidth=1.4, linestyle="--", label=f"median {median:.3f}")
    ax.legend(frameon=False, fontsize=8)
    control = INTERPRETATION_MODE == "unlabeled_control"
    ax.set_title("A  Apparent gene-NTR distribution" if control else "A  Gene NTR distribution", fontweight="bold", color=INK)
    ax.set_xlabel("Apparent posterior NTR" if control else "Posterior mean NTR")
    ax.set_ylabel("Genes")

    ax = axes[1]
    color_limit = float(supported["interval_width"].quantile(0.98))
    scatter = ax.scatter(
        np.log10(supported[read_col] + 1), supported[mean_col],
        c=supported["interval_width"].clip(upper=color_limit), cmap=NPG_GENERAL,
        s=7, alpha=0.72, linewidths=0, rasterized=True,
    )
    colorbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.02)
    colorbar.set_label("90% interval width (capped p98)", fontsize=8)
    ax.set_title("B  Coverage and posterior precision", fontweight="bold", color=INK)
    ax.set_xlabel("log10(read count + 1)")
    ax.set_ylabel("Apparent posterior NTR" if control else "Posterior mean NTR")
    ax.set_ylim(0, 1)

    ax = axes[2]
    labels = gene_labels(top)
    means = top[mean_col].to_numpy(float)
    lower = top[lower_col].to_numpy(float)
    upper = top[upper_col].to_numpy(float)
    err_low, err_high, plot_lower, plot_upper, interval_adjusted = interval_errors(means, lower, upper)
    y = np.arange(len(top))
    ax.errorbar(means, y, xerr=[err_low, err_high], fmt="o", color=BLUE, ecolor=CYAN, capsize=2.5, markersize=4.5)
    ax.set_yticks(y, labels, fontsize=8)
    if len(top):
        pad = max((plot_upper.max() - plot_lower.min()) * 0.12, 0.006)
        ax.set_xlim(max(0, plot_lower.min() - pad), min(1, plot_upper.max() + pad))
    ax.set_title(f"C  Top {len(top)} high-read genes", fontweight="bold", color=INK)
    ax.set_xlabel("Apparent NTR (90% credible interval)" if control else "Posterior mean NTR (90% credible interval)")
    style_axes(axes)
    fig.suptitle(
        (f"GRAND-SLAM unlabeled-control apparent NTR · support ≥ {min_reads:,} reads"
         if control else f"GRAND-SLAM tissue-only gene NTR · support ≥ {min_reads:,} reads"),
        fontsize=13, fontweight="bold", color=INK,
    )
    stem = "03_grand_slam_gene_ntr"
    save_figure(fig, output_dir, stem, dpi)
    return stem, {
        "min_reads": min_reads,
        "n_genes_ge_min_reads": int(len(supported)),
        "median_supported_gene_mean_ntr": median,
        "display_interval_adjustments": int(interval_adjusted.sum()),
    }


def plot_grand_model(prefix, output_dir, dpi):
    ntr_path, ntr = read_grand_aux(prefix, "ntrstat")
    rates_path, rates = read_grand_aux(prefix, "rates")
    mismatch_path, mismatches = read_grand_aux(prefix, "mismatches")
    if ntr.empty or rates.empty or mismatches.empty:
        return None, {}, [str(path) for path in (ntr_path, rates_path, mismatch_path)]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.1), facecolor="white", constrained_layout=True)
    ax = axes[0]
    category_col = "Category" if "Category" in ntr else "Type"
    ntr_columns = {
        "mean": "Mean" if "Mean" in ntr else "ntr",
        "lower": "0.05 quantile" if "0.05 quantile" in ntr else "ntr_lower",
        "upper": "0.95 quantile" if "0.95 quantile" in ntr else "ntr_upper",
    }
    ntr = numeric(ntr, list(ntr_columns.values()))
    ntr = ntr[ntr[category_col].isin(["Exonic", "Intronic"])].copy()
    y = np.arange(len(ntr))
    means = ntr[ntr_columns["mean"]].to_numpy(float)
    lower = ntr[ntr_columns["lower"]].to_numpy(float)
    upper = ntr[ntr_columns["upper"]].to_numpy(float)
    err_low, err_high, plot_lower, plot_upper, interval_adjusted = interval_errors(means, lower, upper)
    ax.errorbar(means, y, xerr=[err_low, err_high], fmt="o", color=BLUE, ecolor=CYAN, capsize=4, markersize=7)
    ax.set_yticks(y, ntr[category_col])
    pad = max((plot_upper.max() - plot_lower.min()) * 0.35, 0.002)
    ax.set_xlim(max(0, plot_lower.min() - pad), min(1, plot_upper.max() + pad))
    for value, position, lo, hi, changed in zip(means, y, lower, upper, interval_adjusted):
        marker = " †" if changed else ""
        ax.annotate(f"{value:.5f} [{lo:.5f}, {hi:.5f}]{marker}", (value, position), xytext=(8, 0), textcoords="offset points", va="center", fontsize=7.5, color=MUTED)
    control = INTERPRETATION_MODE == "unlabeled_control"
    ax.set_title("A  Apparent background-NTR model" if control else "A  Global NTR model", fontweight="bold", color=INK)
    ax.set_xlabel("NTR (90% credible interval)")

    ax = axes[1]
    rate_label_col = "Type" if "Type" in rates else "Rate"
    rate_value_col = next(column for column in rates.columns if column != rate_label_col)
    rates[rate_label_col] = rates[rate_label_col].astype(str).str.replace("single_", "", regex=False)
    rates = numeric(rates, [rate_value_col])
    rates = rates[rates[rate_label_col].isin(["old", "new"])].copy()
    colors = [BLUE if value == "old" else RED for value in rates[rate_label_col]]
    bars = ax.bar(rates[rate_label_col], rates[rate_value_col], color=colors, width=0.58)
    ax.set_yscale("log")
    ax.set_ylabel("Estimated conversion rate")
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    for bar, value in zip(bars, rates[rate_value_col]):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.22, f"{value:.3%}", ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_ylim(max(rates[rate_value_col].min() / 3, 1e-5), min(rates[rate_value_col].max() * 4, 1))
    ax.set_title("B  Background mixture conversion model" if control else "B  Old/new conversion model", fontweight="bold", color=INK)
    ax.set_xlabel("Fitted mixture component" if control else "RNA population")

    mismatches = numeric(mismatches, ["Coverage", "Mismatches"])
    first = mismatches[mismatches["Orientation"].astype(str).str.lower() == "first"].copy()
    first = first[first["Category"].isin(["Exonic", "Intronic"])]
    first["mutation"] = first["Genomic"].astype(str) + ">" + first["Read"].astype(str)
    grouped = first.groupby(["Category", "mutation"], as_index=False)[["Coverage", "Mismatches"]].sum()
    grouped["rate"] = grouped["Mismatches"] / grouped["Coverage"].replace(0, np.nan)
    style_axes(axes)
    fig.suptitle("GRAND-SLAM unlabeled-control model QC" if control else "GRAND-SLAM tissue-only model QC", fontsize=13, fontweight="bold", color=INK)
    if interval_adjusted.any():
        fig.text(0.5, 0.005, "† Input interval did not enclose its point estimate; error bar expanded for display only. Raw values are unchanged.", ha="center", fontsize=7.2, color=MUTED)
    stem = "04_grand_slam_model_qc"
    save_figure(fig, output_dir, stem, dpi)

    model = {
        "global_ntr": [
            {"category": row[category_col], "mean": row[ntr_columns["mean"]], "lower": row[ntr_columns["lower"]], "upper": row[ntr_columns["upper"]]}
            for _, row in ntr.iterrows()
        ],
        "conversion_rates": [
            {"type": row[rate_label_col], "rate": row[rate_value_col]} for _, row in rates.iterrows()
        ],
        "first_read_tc_rates": grouped[grouped["mutation"] == "T>C"][["Category", "rate"]].to_dict("records"),
        "display_interval_adjustments": int(interval_adjusted.sum()),
    }
    return stem, model, [str(path) for path in (ntr_path, rates_path, mismatch_path)]


def plot_tc_mismatch_enrichment(slam, prefix, output_dir, dpi, min_reads):
    """Show whether T>C is enriched over the other substitution classes."""
    required = {"ReadCount", "ConversionRate", "CoverageOnTs", "ConversionsOnTs"}
    if slam.empty or not required.issubset(slam.columns):
        return None, {}
    work = numeric(slam, list(required))
    supported = work[(work["ReadCount"] >= min_reads) & work["ConversionRate"].notna()]
    coverage = float(work["CoverageOnTs"].sum())
    if supported.empty or not coverage:
        return None, {}
    slam_weighted = safe_ratio(work["ConversionsOnTs"].sum(), coverage)
    slam_median = float(supported["ConversionRate"].median())

    mismatch_path, mismatches = read_grand_aux(prefix, "mismatches")
    if mismatches.empty or not {"Category", "Orientation", "Genomic", "Read", "Coverage", "Mismatches"}.issubset(mismatches.columns):
        return None, {}
    mismatches = numeric(mismatches, ["Coverage", "Mismatches"])
    first = mismatches[
        (mismatches["Orientation"].astype(str).str.lower() == "first")
        & mismatches["Category"].isin(["Exonic", "Intronic"])
    ].copy()
    first["mutation"] = first["Genomic"].astype(str) + ">" + first["Read"].astype(str)
    spectrum = first.groupby(["Category", "mutation"], as_index=False)[["Coverage", "Mismatches"]].sum()
    spectrum["rate"] = spectrum["Mismatches"] / spectrum["Coverage"].replace(0, np.nan)
    spectrum = spectrum.dropna(subset=["rate"])
    if spectrum.empty or "T>C" not in set(spectrum["mutation"]):
        return None, {}

    channel_order = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
    spectrum["channel_order"] = spectrum["mutation"].map({channel: index for index, channel in enumerate(channel_order)})
    complete_spectrum = spectrum.sort_values(["Category", "channel_order"])[
        ["Category", "mutation", "Coverage", "Mismatches", "rate"]
    ].rename(columns={"Category": "region", "mutation": "channel", "Coverage": "coverage", "Mismatches": "mismatches", "rate": "mismatch_rate"})
    tc_rates, enrichment = {}, {}
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2), facecolor="white", constrained_layout=True)
    for ax, category, neutral in zip(axes[:2], ["Exonic", "Intronic"], [BLUE, CYAN]):
        data = spectrum[spectrum["Category"] == category].copy()
        tc = float(data.loc[data["mutation"] == "T>C", "rate"].iloc[0])
        background = data[data["mutation"] != "T>C"]["rate"]
        tc_rates[category] = tc
        enrichment[category] = {
            "tc_rate": tc,
            "non_tc_median_rate": float(background.median()),
            "non_tc_max_rate": float(background.max()),
            "fold_over_non_tc_median": float(tc / background.median()),
            "fold_over_non_tc_max": float(tc / background.max()),
            "rank_among_12_substitutions": int(data["rate"].rank(method="min", ascending=False)[data["mutation"] == "T>C"].iloc[0]),
        }
        order = ["T>C"] + data[data["mutation"] != "T>C"].sort_values("rate", ascending=False)["mutation"].tolist()
        data = data.set_index("mutation").reindex(order).reset_index()
        y = np.arange(len(data))
        colors = [RED if mutation == "T>C" else neutral for mutation in data["mutation"]]
        ax.hlines(y, 1e-5, data["rate"], color="#D8E1E5", linewidth=2)
        ax.scatter(data["rate"], y, c=colors, s=[52 if mutation == "T>C" else 28 for mutation in data["mutation"]], zorder=3)
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.set_yticks(y, data["mutation"], fontsize=8)
        ax.invert_yaxis()
        ax.annotate(f"T>C {tc:.3%}", (tc, 0), xytext=(7, 0), textcoords="offset points", va="center", fontsize=8.5, color=RED, fontweight="bold")
        ax.set_xlim(max(data["rate"].min() / 2, 1e-5), tc * 2.2)
        ax.set_title(f"{'A' if category == 'Exonic' else 'B'}  {category} mismatch spectrum", fontweight="bold", color=INK)
        ax.set_xlabel("Mismatch rate (log scale)")

    ax = axes[2]
    categories = ["Exonic", "Intronic"]
    x = np.arange(len(categories))
    median_fold = [enrichment[category]["fold_over_non_tc_median"] for category in categories]
    max_fold = [enrichment[category]["fold_over_non_tc_max"] for category in categories]
    width = 0.34
    bars_a = ax.bar(x - width / 2, median_fold, width, color=RED, label="vs non-T>C median")
    bars_b = ax.bar(x + width / 2, max_fold, width, color="#F39B7F", label="vs largest non-T>C")
    for bars in (bars_a, bars_b):
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.35, f"{bar.get_height():.1f}×", ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(x, categories)
    ax.set_ylim(0, max(median_fold) * 1.22)
    ax.set_ylabel("T>C fold enrichment")
    ax.set_title("C  T>C enrichment over background", fontweight="bold", color=INK)
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.02, 0.04,
        f"Independent slamdunk UTR T>C: weighted {slam_weighted:.3%};\nmedian {slam_median:.3%} (≥{min_reads} reads).",
        transform=ax.transAxes, fontsize=8, color=MUTED, va="bottom",
        bbox={"facecolor": "white", "edgecolor": "#D8E1E5", "boxstyle": "round,pad=0.45"},
    )
    style_axes(axes)
    fig.suptitle("T>C enrichment over the complete substitution spectrum", fontsize=13, fontweight="bold", color=INK)
    stem = "05_tc_mismatch_enrichment_over_background"
    save_figure(fig, output_dir, stem, dpi)
    return stem, {
        "slamdunk_weighted_utr_tc_rate": slam_weighted,
        "slamdunk_supported_utr_median_tc_rate": slam_median,
        "grand_slam_first_read_enrichment": enrichment,
        "complete_12_channel_spectrum": complete_spectrum.to_dict("records"),
        "sources": [str(mismatch_path)],
        "interpretation": "Fold enrichment is descriptive; no replicate-level significance test was performed.",
    }


def safe_ratio(numerator, denominator):
    return float(numerator / denominator) if denominator and np.isfinite(denominator) else None


def main():
    global INTERPRETATION_MODE
    parser = argparse.ArgumentParser()
    parser.add_argument("--slamdunk-dir", required=True)
    parser.add_argument("--grand-prefix", required=True)
    parser.add_argument("--mask-qc", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--figure-dir")
    parser.add_argument("--provenance")
    parser.add_argument("--spectrum-table")
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--slamdunk-min-reads", type=int, default=100)
    parser.add_argument("--grand-min-reads", type=int, default=100)
    parser.add_argument("--top-genes", type=int, default=15)
    parser.add_argument("--interpretation-mode", choices=("standard", "unlabeled_control"), default="standard")
    args = parser.parse_args()
    INTERPRETATION_MODE = args.interpretation_mode

    slamdunk = Path(args.slamdunk_dir)
    prefix = Path(args.grand_prefix)
    tcounts = sorted(slamdunk.rglob("*_tcount.tsv"))
    slam = read_slamdunk(tcounts)
    grand_path, grand, mapping = read_grand(prefix)
    output_dir = Path(args.figure_dir) if args.figure_dir else None
    stems, slam_plot_stats, grand_plot_stats, model_stats, enrichment_stats, aux_paths = [], {}, {}, {}, {}, []
    if output_dir:
        stale_stems = [
            "01_slamdunk_conversion_qc", "02_grand_slam_gene_ntr",
            "03_grand_slam_model_qc", "04_tc_mismatch_enrichment_over_background",
            "04_slamdunk_grand_slam_tc_rate_comparison",
        ]
        for stale_stem in stale_stems:
            for suffix in ("png", "pdf"):
                stale = output_dir / f"{stale_stem}.{suffix}"
                if stale.is_file():
                    stale.unlink()
        stem, slam_plot_stats = plot_slamdunk(slam, grand, output_dir, args.dpi, args.slamdunk_min_reads, args.top_genes)
        if stem:
            stems.append(stem)
        stem, grand_plot_stats = plot_grand_genes(grand, mapping, output_dir, args.dpi, args.grand_min_reads, args.top_genes)
        if stem:
            stems.append(stem)
        stem, model_stats, aux_paths = plot_grand_model(prefix, output_dir, args.dpi)
        if stem:
            stems.append(stem)
        stem, enrichment_stats = plot_tc_mismatch_enrichment(slam, prefix, output_dir, args.dpi, args.slamdunk_min_reads)
        if stem:
            stems.append(stem)

    top_genes_rows = []
    if not grand.empty and mapping.get("readcount") and mapping.get("mean"):
        columns = [column for column in ["Gene", "Symbol", mapping["readcount"], mapping["mean"]] if column in grand]
        top_genes_rows = grand.nlargest(20, mapping["readcount"])[columns].rename(
            columns={mapping["readcount"]: "readcount", mapping["mean"]: "mean_ntr"}
        ).to_dict("records")

    slam_numeric = numeric(slam, ["ReadCount", "ConversionRate", "CoverageOnTs", "ConversionsOnTs", "TcReadCount"])
    read_sum = float(slam_numeric["ReadCount"].sum()) if "ReadCount" in slam_numeric else 0
    coverage_sum = float(slam_numeric["CoverageOnTs"].sum()) if "CoverageOnTs" in slam_numeric else 0
    summary = {
        "tissue_mask": json.loads(Path(args.mask_qc).read_text(encoding="utf-8")),
        "slamdunk": {
            "tcount_files": table_meta(tcounts),
            "snp_vcfs": [str(path) for path in slamdunk.rglob("*.vcf")],
            "n_tcount_rows": int(len(slam)),
            "median_conversion_rate": float(slam_numeric["ConversionRate"].median()) if "ConversionRate" in slam_numeric else None,
            "total_counted_reads": int(read_sum),
            "weighted_conversion_rate": safe_ratio(slam_numeric["ConversionsOnTs"].sum(), coverage_sum) if "ConversionsOnTs" in slam_numeric else None,
            "tc_read_fraction": safe_ratio(slam_numeric["TcReadCount"].sum(), read_sum) if "TcReadCount" in slam_numeric else None,
            **slam_plot_stats,
        },
        "grand_slam": {
            "table": str(grand_path) if grand_path else None,
            "available": grand_path is not None,
            "n_genes": int(len(grand)),
            "median_gene_mean_ntr": float(grand[mapping["mean"]].median()) if mapping.get("mean") else None,
            "top_genes_by_readcount": top_genes_rows,
            "model_qc": model_stats,
            "warning": (
                "Unlabeled negative-control apparent NTR; this is background or false-positive signal, not newly synthesized RNA."
                if INTERPRETATION_MODE == "unlabeled_control" else
                "Exploratory single-sample background estimation; no matched no4sU control or validated error model."
            ),
            **grand_plot_stats,
        },
        "figures": [f"{stem}.png" for stem in stems],
        "tc_mismatch_enrichment": enrichment_stats,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.spectrum_table and enrichment_stats.get("complete_12_channel_spectrum"):
        spectrum_table = Path(args.spectrum_table)
        spectrum_table.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(enrichment_stats["complete_12_channel_spectrum"]).to_csv(spectrum_table, sep="\t", index=False)

    provenance = {
        "scope": "tissue-only bulk SLAM",
        "interpretation_mode": INTERPRETATION_MODE,
        "raw_statistics_modified": False,
        "figures": [
            {
                "file": "02_slamdunk_conversion_qc.png",
                "sources": [str(path) for path in tcounts] + ([str(grand_path)] if grand_path else []),
                "filters": {"minimum_read_count": args.slamdunk_min_reads},
                "display_only": {"conversion_rate_ceiling": "98th percentile; statistics use all supported rows"},
                "aggregation": "one representative UTR per gene for panel C, selected by maximum CoverageOnTs",
            },
            {
                "file": "03_grand_slam_gene_ntr.png",
                "sources": [str(grand_path)] if grand_path else [],
                "filters": {"minimum_read_count": args.grand_min_reads},
                "display_only": {"interval_width_color_ceiling": "98th percentile"},
                "aggregation": "gene-level GRAND-SLAM posterior output",
            },
            {
                "file": "04_grand_slam_model_qc.png",
                "sources": aux_paths,
                "filters": {},
                "display_only": {"conversion_rate_axis": "log scale"},
                "aggregation": "GRAND-SLAM global NTR posterior and modeled old/new conversion rates",
            },
            {
                "file": "05_tc_mismatch_enrichment_over_background.png",
                "sources": [str(path) for path in tcounts] + enrichment_stats.get("sources", []),
                "filters": {"slamdunk_minimum_read_count_for_median": args.slamdunk_min_reads},
                "display_only": {"mismatch_spectrum_axis": "log scale; T>C highlighted in red"},
                "aggregation": "GRAND-SLAM first-read mismatch counts and coverage summed by annotation and substitution; T>C fold enrichment calculated against the median and maximum of the other 11 substitutions. slamdunk contributes an independent UTR T>C estimate only because it does not output a 12-channel mismatch spectrum",
            },
        ],
    }
    if args.provenance:
        path = Path(args.provenance)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
