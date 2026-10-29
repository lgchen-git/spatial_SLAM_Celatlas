#!/usr/bin/env python3
"""Build reviewed recovery-bias, marker interpretation and cross-check outputs."""
from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pearsonr, spearmanr


DOMAIN_COLORS = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                 "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"]
INK = "#25313A"

MARKER_PANELS = {
    "proximal tubule-like": ["Slc34a1", "Lrp2", "Cubn", "Slc5a2", "Aqp1"],
    "thick ascending limb-like": ["Umod", "Slc12a1", "Kcnj1", "Clcnka"],
    "distal convoluted tubule-like": ["Slc12a3", "Pvalb", "Calb1"],
    "collecting duct principal-like": ["Aqp2", "Fxyd4", "Scnn1g", "Avpr2"],
    "podocyte-like": ["Nphs1", "Nphs2", "Podxl", "Synpo"],
    "endothelial-like": ["Pecam1", "Kdr", "Emcn", "Klf2", "Vwf"],
    "immune-like": ["Ptprc", "Lst1", "Tyrobp", "C1qa", "Cd74"],
    "stromal-like": ["Col1a1", "Col3a1", "Dcn", "Pdgfra", "Col1a2"],
    "erythroid-like": ["Hbb-bs", "Hbb-bt", "Hba-a1", "Alas2"],
    "injury-associated epithelial-like": ["Havcr1", "Lcn2", "Krt8", "Krt18", "Spp1"],
}


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def correlations(x, y):
    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 3:
        return {"n": int(finite.sum()), "spearman_rho": np.nan, "pearson_r": np.nan}
    return {
        "n": int(finite.sum()),
        "spearman_rho": float(spearmanr(x[finite], y[finite]).statistic),
        "pearson_r": float(pearsonr(x[finite], y[finite]).statistic),
    }


def recovery_qc(domains, bin_summary, positions_path, tables, figures):
    frame = domains[["bin_id", "bin_x", "bin_y", "domain", "umi_count"]].copy()
    frame = frame.merge(bin_summary[["bin_id", "total", "corrected_ntr"]], on="bin_id", how="left")
    frame[["total", "corrected_ntr"]] = frame[["total", "corrected_ntr"]].fillna(0)
    frame = frame.rename(columns={"umi_count": "celatlas_total_umi", "total": "slam_molecules"})
    frame["relative_slam_recovery_index"] = (
        (frame["slam_molecules"] / max(frame["slam_molecules"].sum(), 1))
        / (frame["celatlas_total_umi"] / max(frame["celatlas_total_umi"].sum(), 1)).replace(0, np.nan)
    )
    frame["log2_relative_slam_recovery_index"] = np.log2(frame["relative_slam_recovery_index"].clip(lower=1e-6))
    positions = pd.read_csv(positions_path, header=None,
        names=["bin_id", "in_tissue", "row", "col", "pxl_row", "pxl_col"], dtype={"bin_id": str})
    coords = {(int(r.row), int(r.col)) for r in positions.itertuples(index=False)}
    edge = {}
    for r in positions.itertuples(index=False):
        neighbors = sum((int(r.row) + dr, int(r.col) + dc) in coords for dr, dc in ((-1,0),(1,0),(0,-1),(0,1)))
        edge[str(r.bin_id)] = "edge" if neighbors < 4 else "center"
    frame["tissue_location"] = frame["bin_id"].map(edge).fillna("unknown")
    frame["celatlas_depth_quartile"] = pd.qcut(
        frame["celatlas_total_umi"].rank(method="first"), 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]
    )
    frame.to_csv(tables / "bin_slam_recovery_qc.tsv.gz", sep="\t", index=False, compression="gzip")
    by_domain = frame.groupby("domain", as_index=False).agg(
        n_bins=("bin_id", "size"), celatlas_total_umi=("celatlas_total_umi", "sum"),
        slam_molecules=("slam_molecules", "sum"), median_recovery_index=("relative_slam_recovery_index", "median"),
        median_ntr=("corrected_ntr", "median"),
    )
    by_domain.to_csv(tables / "domain_slam_recovery_qc.tsv", sep="\t", index=False)
    frame.groupby(["tissue_location", "celatlas_depth_quartile"], observed=True, as_index=False).agg(
        n_bins=("bin_id", "size"), median_celatlas_umi=("celatlas_total_umi", "median"),
        median_slam_molecules=("slam_molecules", "median"),
        median_relative_recovery=("relative_slam_recovery_index", "median"), median_ntr=("corrected_ntr", "median"),
    ).to_csv(tables / "recovery_edge_depth_summary.tsv", sep="\t", index=False)
    tests = []
    for label, x in (("relative_recovery_vs_ntr", frame["relative_slam_recovery_index"].to_numpy()),
                     ("celatlas_depth_vs_ntr", frame["celatlas_total_umi"].to_numpy()),
                     ("slam_depth_vs_ntr", frame["slam_molecules"].to_numpy())):
        tests.append({"comparison": label, **correlations(x, frame["corrected_ntr"].to_numpy())})
    pd.DataFrame(tests).to_csv(tables / "recovery_ntr_correlations.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True, facecolor="white")
    specs = [("celatlas_total_umi", "Celatlas total UMI", "viridis"),
             ("slam_molecules", "SLAM molecule density", "magma"),
             ("log2_relative_slam_recovery_index", "Relative SLAM recovery index (log2)", "coolwarm")]
    for ax, (column, title, cmap) in zip(axes, specs):
        values = frame[column].to_numpy(float)
        if column != "log2_relative_slam_recovery_index":
            values = np.log1p(values)
        vmax = np.nanquantile(np.abs(values), .99) if column.startswith("log2") else np.nanquantile(values, .99)
        kwargs = {"vmin": -vmax, "vmax": vmax} if column.startswith("log2") else {"vmin": 0, "vmax": vmax}
        artist = ax.scatter(frame["bin_x"], frame["bin_y"], c=values, s=4, cmap=cmap, linewidths=0, **kwargs)
        ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(title, fontweight="bold", color=INK); ax.axis("off")
        fig.colorbar(artist, ax=ax, fraction=.03, pad=.01)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"07_barcode_recovery_bias.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)
    return frame, by_domain, tests


def read_gene_info(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        frame = pd.read_csv(handle, sep="\t", dtype=str)
    frame = frame[frame["#tax_id"] == "10090"]
    mapping = dict(zip(frame["Symbol"].str.upper(), frame["GeneID"]))
    for row in frame.itertuples(index=False):
        synonyms = getattr(row, "Synonyms", "-")
        if synonyms and synonyms != "-":
            for synonym in synonyms.split("|"):
                mapping.setdefault(synonym.upper(), str(row.GeneID))
    return mapping


def read_gmt(path):
    pathways = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            pathways[fields[0]] = set(fields[2:])
    return pathways


def enrich(groups, universe, pathways):
    universe = set(universe)
    rows = []
    for group, selected in groups.items():
        selected = set(selected) & universe
        for pathway, genes in pathways.items():
            genes = genes & universe
            overlap = selected & genes
            if len(overlap) < 2 or len(genes) < 5 or not selected:
                continue
            pvalue = hypergeom.sf(len(overlap) - 1, len(universe), len(genes), len(selected))
            rows.append({"group": group, "pathway": pathway, "overlap": len(overlap),
                         "selected_genes": len(selected), "pathway_genes": len(genes),
                         "gene_ratio": len(overlap) / len(selected), "pvalue": pvalue,
                         "overlap_entrez": ",".join(sorted(overlap, key=lambda x: int(x)))})
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["group", "pathway", "overlap", "selected_genes", "pathway_genes", "gene_ratio", "pvalue", "adjusted_pvalue", "overlap_entrez"])
    result["adjusted_pvalue"] = np.nan
    for _, index in result.groupby("group").groups.items():
        result.loc[index, "adjusted_pvalue"] = bh_adjust(result.loc[index, "pvalue"])
    return result.sort_values(["group", "adjusted_pvalue", "pvalue"])


def markers_annotations_pathways(markers, domain_gene, marker_universe, gene_info_path, gmt_path, tables, figures):
    eligible = markers[(markers.adjusted_pvalue < .05) & (markers.logfoldchange > .5) &
                       (markers.pct_in_domain >= .05) & (markers.pct_difference >= .05)].copy()
    annotations = []
    for domain, subset in markers.groupby("domain"):
        scores = {}
        found = {}
        for label, panel in MARKER_PANELS.items():
            panel_upper = {gene.upper(): gene for gene in panel}
            evidence = subset[
                subset.gene_name.str.upper().isin(panel_upper)
                & (subset.adjusted_pvalue < .05) & (subset.logfoldchange > .25)
            ].sort_values("logfoldchange", ascending=False).drop_duplicates(subset=["gene_name"])
            hits = [panel_upper[value] for value in evidence.gene_name.str.upper()]
            found[label] = hits
            scores[label] = float(evidence.logfoldchange.clip(lower=0).sum())
        best = max(scores, key=scores.get)
        hits = found[best]
        if not hits:
            label, confidence, notes = f"D{domain} unresolved expression domain", "low", "No pre-specified compartment panel had positive marker support."
        else:
            label = best
            confidence = "high" if len(hits) >= 3 else "medium" if len(hits) >= 2 else "low"
            notes = "Putative label from a pre-specified kidney compartment marker panel; not a validated cell type."
        annotations.append({"domain": str(domain), "putative_annotation": label,
                            "supporting_markers": ", ".join(hits), "supporting_pathways": "",
                            "confidence": confidence, "notes": notes})

    symbol_to_entrez = read_gene_info(gene_info_path)
    pathways = read_gmt(gmt_path)
    markers["entrez"] = markers.gene_name.str.upper().map(symbol_to_entrez)
    eligible["entrez"] = eligible.gene_name.str.upper().map(symbol_to_entrez)
    domain_gene["entrez"] = domain_gene.gene_name.str.upper().map(symbol_to_entrez)
    marker_universe["entrez"] = marker_universe.gene_name.str.upper().map(symbol_to_entrez)
    universe = set(marker_universe.entrez.dropna())
    marker_groups = {f"D{domain}": set(sub.entrez.dropna()) for domain, sub in eligible.groupby("domain")}
    candidates = domain_gene[(domain_gene.total >= 500) & (domain_gene.n_bins >= 10)].copy()
    candidates["abs_dynamic_score"] = candidates.dynamic_score.abs()
    candidate_groups = {f"D{domain}": set(sub.nlargest(200, "abs_dynamic_score").entrez.dropna())
                        for domain, sub in candidates.groupby("domain")}
    marker_enrichment = enrich(marker_groups, universe, pathways)
    ntr_enrichment = enrich(candidate_groups, universe, pathways)
    marker_enrichment.to_csv(tables / "domain_marker_pathways.tsv", sep="\t", index=False)
    ntr_enrichment.to_csv(tables / "domain_ntr_shift_pathways.tsv", sep="\t", index=False)
    top_pathway = marker_enrichment[marker_enrichment.adjusted_pvalue < .1].groupby("group").head(2)
    annotations_frame = pd.DataFrame(annotations)
    annotations_frame["supporting_pathways"] = annotations_frame.domain.map(
        top_pathway.groupby(top_pathway.group.str.removeprefix("D"))["pathway"].apply(lambda x: "; ".join(x)).to_dict()
    ).fillna("")
    annotations_frame.to_csv(tables / "domain_annotation_evidence.tsv", sep="\t", index=False)
    plot_enrichment(marker_enrichment, figures / "10_domain_marker_pathways", "Domain marker Reactome enrichment")
    plot_enrichment(ntr_enrichment, figures / "11_domain_ntr_shift_pathways", "Candidate NTR-shift Reactome enrichment")
    return annotations_frame, marker_enrichment, ntr_enrichment, candidates


def plot_enrichment(frame, output, title):
    selected = frame[frame.adjusted_pvalue < .1].groupby("group").head(3) if not frame.empty else frame
    if selected.empty:
        selected = frame.head(20)
    if selected.empty:
        return
    selected = selected.tail(30).copy()
    selected["label"] = selected["group"] + " | " + selected["pathway"].str.replace("_", " ").str.slice(0, 58)
    fig, ax = plt.subplots(figsize=(10, max(5, .28 * len(selected) + 1.5)), facecolor="white")
    values = -np.log10(selected.adjusted_pvalue.clip(lower=1e-300))
    artist = ax.scatter(values, np.arange(len(selected)), s=25 + 35 * selected.overlap,
                        c=selected.gene_ratio, cmap="viridis", linewidths=0)
    ax.set_yticks(np.arange(len(selected)), selected.label, fontsize=7)
    ax.set_xlabel("-log10 adjusted P"); ax.set_title(title, fontweight="bold", color=INK)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="x", color="#E5E7EB", lw=.5)
    fig.colorbar(artist, ax=ax, label="Gene ratio", fraction=.025, pad=.02)
    for suffix in ("png", "pdf"):
        fig.savefig(output.with_suffix(f".{suffix}"), dpi=450, bbox_inches="tight")
    plt.close(fig)


def classify_story_types(candidates):
    """Assign A/B/C labels only when the quantitative pattern supports them."""
    candidates = candidates.copy()
    candidates["mean_total_per_bin"] = candidates.total / candidates.n_bins.clip(lower=1)
    candidates["mean_new_per_bin"] = candidates.expected_new / candidates.n_bins.clip(lower=1)
    gene_patterns = candidates.groupby("gene_id").agg(
        domain_total_mean=("mean_total_per_bin", "mean"), domain_total_sd=("mean_total_per_bin", "std"),
        domain_total_median=("mean_total_per_bin", "median"),
        domain_new_median=("mean_new_per_bin", "median"),
        domain_ntr_min=("ntr", "min"), domain_ntr_max=("ntr", "max"),
    )
    gene_patterns["domain_total_cv"] = gene_patterns.domain_total_sd.fillna(0) / gene_patterns.domain_total_mean.clip(lower=1)
    gene_patterns["domain_ntr_range"] = gene_patterns.domain_ntr_max - gene_patterns.domain_ntr_min
    candidates = candidates.merge(
        gene_patterns[["domain_total_cv", "domain_total_median", "domain_new_median", "domain_ntr_range"]],
        left_on="gene_id", right_index=True,
    )
    candidates["total_enrichment_ratio"] = (
        candidates.mean_total_per_bin / candidates.domain_total_median.clip(lower=1e-12)
    )
    candidates["new_enrichment_ratio"] = (
        candidates.mean_new_per_bin / candidates.domain_new_median.clip(lower=1e-12)
    )
    type_c = (candidates.domain_total_cv <= .35) & (candidates.domain_ntr_range >= .10)
    type_a = (
        (candidates.total_enrichment_ratio >= 1.25)
        & (candidates.new_enrichment_ratio >= 1.25)
        & (candidates.delta_ntr >= .05)
    )
    type_b = (
        (candidates.total_enrichment_ratio >= 1.25)
        & (candidates.new_enrichment_ratio < candidates.total_enrichment_ratio)
        & (candidates.delta_ntr <= -.05)
    )
    candidates["story_type"] = np.select(
        [type_c, type_a, type_b],
        [
            "C: similar Total RNA, different NTR",
            "A: Total/New/NTR enriched",
            "B: Total enriched with lower NTR",
        ],
        default="unclassified candidate NTR shift",
    )
    return candidates


def select_representatives(candidates, tables):
    candidates = classify_story_types(candidates)
    candidates["relative_uncertainty"] = (candidates.ntr_ci_upper - candidates.ntr_ci_lower).fillna(1)
    candidates["eligible"] = (candidates.total >= 500) & (candidates.n_bins >= 10) & (candidates.relative_uncertainty <= .25)
    candidates["selection_score"] = candidates.dynamic_score.abs() * np.log1p(candidates.n_bins) / (1 + candidates.relative_uncertainty)
    candidates["story_priority"] = np.select(
        [candidates.story_type.str.startswith("C:"), candidates.story_type.str.match(r"^[AB]:")],
        [2, 1], default=0,
    )
    candidates = candidates.sort_values(["story_priority", "selection_score"], ascending=False)
    candidates["candidate_rank"] = np.arange(1, len(candidates) + 1)
    chosen, used = [], set()
    for row in candidates[candidates.eligible].itertuples(index=False):
        if row.gene_id in used or any(str(item.domain) == str(row.domain) for item in chosen):
            continue
        chosen.append(row); used.add(row.gene_id)
        if len(chosen) == 4:
            break
    selected = pd.DataFrame(chosen)
    if not selected.empty:
        selected["molecule_support"] = selected["total"]
        selected["spatial_coverage_bins"] = selected["n_bins"]
        selected["reason_selected"] = selected.story_type.map(
            lambda story: f"{story}; high uncertainty-adjusted |delta NTR| with >=500 molecules, >=10 bins, and one gene per domain"
        )
    selected.to_csv(tables / "representative_gene_selection.tsv", sep="\t", index=False)
    return selected


def representative_figure(raw, smoothed, selected, domains, figures):
    if selected.empty:
        return
    genes = selected.gene_id.astype(str).tolist()
    raw = raw[raw.gene_id.astype(str).isin(genes)].copy()
    smoothed = smoothed[smoothed.gene_id.astype(str).isin(genes)].copy()
    missing_raw = sorted(set(genes) - set(raw.gene_id.astype(str)))
    missing_smoothed = sorted(set(genes) - set(smoothed.gene_id.astype(str)))
    if missing_raw or missing_smoothed:
        raise ValueError(
            "Representative genes are missing from per-bin display inputs: "
            f"raw={missing_raw}, smoothed={missing_smoothed}"
        )
    for mode, frame, prefix in (("raw", raw, ""), ("smoothed", smoothed, "smoothed_")):
        if frame.empty:
            continue
        fig, axes = plt.subplots(len(genes), 3, figsize=(11.8, 3.4 * len(genes)), constrained_layout=True, facecolor="white")
        axes = np.atleast_2d(axes)
        for row_index, gene in enumerate(genes):
            sub = frame[frame.gene_id.astype(str) == gene]
            name = str(sub.gene_name.iloc[0]) if len(sub) else gene
            columns = ((f"{prefix}total", "Total RNA", "viridis"),
                       (f"{prefix}expected_new", "Expected New RNA", "magma"),
                       (f"{prefix}ntr", "NTR", "plasma"))
            for ax, (column, label, cmap) in zip(axes[row_index], columns):
                ax.scatter(domains.bin_x, domains.bin_y, c="#E5E7EB", s=3, linewidths=0, rasterized=True)
                if sub.empty or column not in sub:
                    ax.axis("off"); continue
                values = sub[column].to_numpy(float)
                vmax = 1 if label == "NTR" else max(float(np.nanquantile(values, .99)), 1e-6)
                artist = ax.scatter(sub.bin_x, sub.bin_y, c=values, s=4, cmap=cmap, vmin=0, vmax=vmax, linewidths=0)
                ax.invert_yaxis(); ax.set_aspect("equal"); ax.axis("off")
                ax.set_title(f"{name} | {label} ({mode})", fontsize=9, fontweight="bold", color=INK)
                fig.colorbar(artist, ax=ax, fraction=.025, pad=.01)
        for suffix in ("png", "pdf"):
            fig.savefig(figures / f"12_representative_total_new_ntr_{mode}.{suffix}", dpi=450, bbox_inches="tight")
        plt.close(fig)


def grand_comparison(global_gene, grand_path, tables, figures):
    spatial = global_gene[["gene_id", "gene_name", "total", "ntr"]].rename(
        columns={"total": "global_total", "ntr": "global_ntr"}
    )
    grand = pd.read_csv(grand_path, sep="\t", compression="infer")
    read_col = next(col for col in grand if col.endswith(" Readcount"))
    ntr_col = next(col for col in grand if col.endswith(" Mean"))
    grand = grand[["Gene", "Symbol", read_col, ntr_col]].rename(columns={"Gene": "gene_id", "Symbol": "grand_symbol", read_col: "grand_reads", ntr_col: "grand_ntr"})
    merged = spatial.merge(grand, on="gene_id", how="inner")
    merged.to_csv(tables / "spatial_vs_grandslam_gene_ntr.tsv", sep="\t", index=False)
    stats = []
    for threshold in (100, 500, 1000):
        sub = merged[(merged.global_total >= threshold) & (merged.grand_reads >= threshold)]
        corr = correlations(sub.grand_ntr.to_numpy(float), sub.global_ntr.to_numpy(float))
        stats.append({"minimum_support_each_method": threshold, **corr,
                      "mae": float(np.mean(np.abs(sub.global_ntr - sub.grand_ntr))) if len(sub) else np.nan,
                      "median_bias_spatial_minus_grand": float(np.median(sub.global_ntr - sub.grand_ntr)) if len(sub) else np.nan})
    pd.DataFrame(stats).to_csv(tables / "spatial_vs_grandslam_gene_ntr_statistics.tsv", sep="\t", index=False)
    plot = merged[(merged.global_total >= 100) & (merged.grand_reads >= 100)].copy()
    if not plot.empty:
        fig, ax = plt.subplots(figsize=(6.4, 5.8), facecolor="white")
        sizes = 5 + 30 * np.log1p(np.minimum(plot.global_total, plot.grand_reads)) / np.log1p(np.minimum(plot.global_total, plot.grand_reads).max())
        ax.scatter(plot.grand_ntr, plot.global_ntr, s=sizes, c="#3C5488", alpha=.45, linewidths=0)
        ax.plot([0,1], [0,1], ls="--", lw=1, color="#9CA3AF")
        label_offsets = [(4, 5), (4, -10), (-42, 5), (-42, -10), (4, 15), (-42, 15)]
        for offset, row in zip(
            label_offsets,
            plot.assign(delta=(plot.global_ntr-plot.grand_ntr).abs()).nlargest(len(label_offsets), "delta").itertuples(index=False),
        ):
            ax.annotate(
                row.gene_name,
                (row.grand_ntr, row.global_ntr),
                xytext=offset,
                textcoords="offset points",
                fontsize=6,
                ha="left" if offset[0] > 0 else "right",
                va="bottom" if offset[1] > 0 else "top",
            )
        ax.set(xlim=(0,1), ylim=(0,1), xlabel="GRAND-SLAM gene NTR", ylabel="Spatial model gene NTR")
        ax.set_title("Independent model-based cross-check", fontweight="bold", color=INK)
        ax.spines[["top","right"]].set_visible(False); ax.grid(color="#E5E7EB", lw=.5)
        for suffix in ("png", "pdf"):
            fig.savefig(figures / f"13_spatial_vs_grandslam_gene_ntr.{suffix}", dpi=450, bbox_inches="tight")
        plt.close(fig)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domains", required=True); parser.add_argument("--bin-summary", required=True)
    parser.add_argument("--domain-gene", required=True); parser.add_argument("--global-gene", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--marker-universe", required=True)
    parser.add_argument("--raw-bin-gene", required=True); parser.add_argument("--smoothed-bin-gene", required=True)
    parser.add_argument("--positions", required=True); parser.add_argument("--gene-info", required=True)
    parser.add_argument("--reactome-gmt", required=True); parser.add_argument("--grand-gene-ntr", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir); tables = output / "01_tables"; figures = output / "02_figures"
    tables.mkdir(parents=True, exist_ok=True); figures.mkdir(parents=True, exist_ok=True)
    domains = pd.read_csv(args.domains, sep="\t", dtype={"bin_id": str}).rename(columns={"cluster": "domain"})
    domains["domain"] = domains.domain.astype(str)
    bin_summary = pd.read_csv(args.bin_summary, sep="\t", compression="infer", dtype={"bin_id": str})
    _, domain_recovery, recovery_tests = recovery_qc(domains, bin_summary, args.positions, tables, figures)
    markers = pd.read_csv(args.markers, sep="\t", dtype={"domain": str, "gene_id": str})
    marker_universe = pd.read_csv(args.marker_universe, sep="\t", dtype={"gene_id": str})
    domain_gene = pd.read_csv(args.domain_gene, sep="\t", dtype={"domain": str, "gene_id": str})
    global_gene = pd.read_csv(args.global_gene, sep="\t", dtype={"gene_id": str})
    annotations, marker_paths, ntr_paths, candidates = markers_annotations_pathways(
        markers, domain_gene, marker_universe, args.gene_info, args.reactome_gmt, tables, figures)
    candidates.to_csv(tables / "candidate_ntr_shift_genes.tsv", sep="\t", index=False)
    selected = select_representatives(candidates, tables)
    raw = pd.read_csv(args.raw_bin_gene, sep="\t", compression="infer", dtype={"gene_id": str, "bin_id": str})
    smoothed = pd.read_csv(args.smoothed_bin_gene, sep="\t", compression="infer", dtype={"gene_id": str, "bin_id": str})
    representative_figure(raw, smoothed, selected, domains, figures)
    grand_stats = grand_comparison(global_gene, args.grand_gene_ntr, tables, figures)
    summary = {"recovery_correlations": recovery_tests, "domain_annotations": len(annotations),
               "marker_pathway_rows": len(marker_paths), "ntr_pathway_rows": len(ntr_paths),
               "candidate_rows": len(candidates), "representative_genes": selected.gene_name.tolist() if not selected.empty else [],
               "grand_comparison": grand_stats,
               "annotation_scope": "putative expression-domain interpretation, not validated cell types"}
    (output / "05_metadata/reviewed_biology_summary.json").parent.mkdir(parents=True, exist_ok=True)
    (output / "05_metadata/reviewed_biology_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
