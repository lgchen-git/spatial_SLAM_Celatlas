#!/usr/bin/env python3
"""Create a delivery-ready total-only Spatial SLAM binned figure package."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm


NPG_CATEGORICAL = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#6F99AD", "#8A7AA9",
    "#C9A227", "#6E7F80", "#D4537E", "#378ADD",
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255",
    "#AA4499", "#661100", "#6699CC", "#AA4466",
    "#4477AA", "#228833", "#EE6677", "#BBBBBB",
]
NPG_GENERAL = LinearSegmentedColormap.from_list("npg_general", ["#E5E7EB", "#4DBBD5", "#3C5488", "#E64B35"])
NPG_ENRICHMENT = LinearSegmentedColormap.from_list("npg_enrichment", ["#E5E7EB", "#91D1C2", "#00A087", "#DC0000"])
NPG_SIGNED = LinearSegmentedColormap.from_list("npg_signed", ["#3C5488", "#F7F7F7", "#E64B35"])
INK = "#25313A"
MUTED = "#64727B"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument("--corrected-bin-gene", required=True)
    parser.add_argument("--bin-summary", help="Pre-aggregated bin summary from spatial_slam_kinetics.py")
    parser.add_argument("--domain-gene-kinetics", required=True)
    parser.add_argument("--smoothed-bin-gene", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bin-size-raw", type=int, default=50)
    parser.add_argument("--coordinate-um", type=float, default=0.5)
    parser.add_argument("--bin-size-um", type=float, default=None)
    parser.add_argument("--spatial-point-size", type=float, default=14.0)
    parser.add_argument("--umap-point-size", type=float, default=2.2)
    parser.add_argument("--interpretation-mode", choices=("standard", "unlabeled_control"), default="standard")
    parser.add_argument("--dpi", type=int, default=450)
    return parser.parse_args()


def save(fig, figures, stem, dpi):
    figures.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"{stem}.{suffix}", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def cluster_order(values):
    return sorted(pd.Series(values).astype(str).unique(), key=lambda value: int(value) if value.isdigit() else value)


def named_colors(levels):
    if len(levels) > len(NPG_CATEGORICAL):
        raise ValueError("More than 32 spatial domains require an explicit palette")
    return dict(zip(levels, NPG_CATEGORICAL))


def clip(values, lower=2, upper=98):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return tuple(np.percentile(finite, [lower, upper])) if finite.size else (0.0, 1.0)


def style_physical(ax, title=None):
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", color=INK, pad=7)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xlabel("Spatial X (raw chip units)", fontsize=8, color=MUTED)
    ax.set_ylabel("Spatial Y (raw chip units; downward)", fontsize=8, color=MUTED)
    ax.tick_params(labelsize=7, colors=MUTED, length=2)
    for spine in ax.spines.values():
        spine.set_visible(False)


def spatial_categories(ax, frame, colors, title, point_size=14.0):
    levels = list(colors)
    code = {level: idx for idx, level in enumerate(levels)}
    values = frame["domain"].astype(str).map(code)
    ax.scatter(frame["bin_x"], frame["bin_y"], c=values, cmap=ListedColormap([colors[level] for level in levels]),
               vmin=-0.5, vmax=len(levels) - 0.5, marker="s", s=point_size, linewidths=0, rasterized=True)
    style_physical(ax, title)


def spatial_continuous(ax, frame, column, title, cmap, label, transform=None, point_size=14.0):
    values = frame[column].to_numpy(dtype=float)
    shown = transform(values) if transform else values
    order = np.argsort(shown)
    vmin, vmax = clip(shown)
    scatter = ax.scatter(frame.iloc[order]["bin_x"], frame.iloc[order]["bin_y"], c=shown[order], cmap=cmap,
                         vmin=vmin, vmax=vmax, marker="s", s=point_size, linewidths=0, rasterized=True)
    style_physical(ax, title)
    cbar = ax.figure.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)


def load_inputs(args):
    domains = pd.read_csv(args.domains, sep="\t")
    domains = domains.rename(columns={"cluster": "domain"})
    domains["bin_id"] = domains["bin_id"].astype(str)
    domains["domain"] = domains["domain"].astype(str)
    if args.bin_summary:
        bins = pd.read_csv(args.bin_summary, sep="\t", compression="gzip")
    else:
        corrected = pd.read_csv(
            args.corrected_bin_gene,
            sep="\t",
            compression="gzip",
            usecols=["bin_id", "bin_x", "bin_y", "total", "expected_new", "expected_old"],
        )
        corrected["bin_id"] = corrected["bin_id"].astype(str)
        bins = corrected.groupby(["bin_id", "bin_x", "bin_y"], as_index=False).agg(
            total=("total", "sum"), expected_new=("expected_new", "sum"), expected_old=("expected_old", "sum")
        )
        bins["corrected_ntr"] = bins["expected_new"] / bins["total"].clip(lower=1e-12)
    bins["bin_id"] = bins["bin_id"].astype(str)
    bins = bins.merge(domains[["bin_id", "domain", "UMAP1", "UMAP2"]], on="bin_id", how="inner")
    kinetics = pd.read_csv(args.domain_gene_kinetics, sep="\t")
    kinetics["domain"] = kinetics["domain"].astype(str)
    smoothed = pd.read_csv(args.smoothed_bin_gene, sep="\t", compression="gzip")
    smoothed["domain"] = smoothed["domain"].astype(str)
    return domains, bins, kinetics, smoothed


def heatmap_data(kinetics, domains, top_per_domain=4, max_genes=28):
    chosen = []
    for domain in domains:
        sub = kinetics[kinetics["domain"] == domain]
        chosen.extend(sub.assign(abs_delta=sub["delta_ntr"].abs()).nlargest(top_per_domain, "abs_delta")["gene_id"])
    genes = list(dict.fromkeys(chosen))[:max_genes]
    labels = kinetics.drop_duplicates("gene_id").set_index("gene_id")["gene_name"].to_dict()
    pivot = kinetics[kinetics["gene_id"].isin(genes)].pivot_table(index="gene_id", columns="domain", values="delta_ntr", aggfunc="mean")
    pivot = pivot.reindex(index=genes, columns=domains).fillna(0.0)
    return pivot, labels


def draw_heatmap(ax, pivot, labels, title, colorbar_label="ΔNTR vs gene mean"):
    values = pivot.to_numpy(dtype=float)
    vmax = max(float(np.nanmax(np.abs(values))), 0.05)
    im = ax.imshow(values, aspect="auto", cmap=NPG_SIGNED, norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_title(title, fontsize=11, fontweight="bold", color=INK, pad=7)
    ax.set_xticks(range(len(pivot.columns)), [f"D{x}" for x in pivot.columns], fontsize=7)
    ax.set_yticks(range(len(pivot.index)), [labels.get(x, x) for x in pivot.index], fontsize=6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(colorbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)


def draw_domain_box(ax, bins, colors, domains, interpretation_mode="standard"):
    data = [bins.loc[bins["domain"] == domain, "corrected_ntr"].to_numpy() for domain in domains]
    box = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.68, medianprops={"color": INK, "linewidth": 1.1})
    for patch, domain in zip(box["boxes"], domains):
        patch.set_facecolor(colors[domain]); patch.set_alpha(0.82); patch.set_edgecolor("white")
    ax.set_xticks(range(1, len(domains) + 1), [f"D{x}" for x in domains], fontsize=8)
    control = interpretation_mode == "unlabeled_control"
    ax.set_ylabel("Apparent NTR" if control else "Corrected NTR", fontsize=8)
    ax.set_title("Domain-resolved background signal" if control else "Domain-resolved RNA dynamics", fontsize=11, fontweight="bold", color=INK, pad=7)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.6)
    ax.tick_params(axis="y", labelsize=7, colors=MUTED)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def main():
    args = parse_args()
    bin_label = f"bin{args.bin_size_raw}"
    control = args.interpretation_mode == "unlabeled_control"
    root = Path(args.output_dir)
    tables, figures, metadata = root / "01_tables", root / "02_figures", root / "05_metadata"
    for path in (tables, figures, metadata):
        path.mkdir(parents=True, exist_ok=True)
    domains_frame, bins, kinetics, smoothed = load_inputs(args)
    domains = cluster_order(bins["domain"])
    colors = named_colors(domains)
    fast_gzip = {"method": "gzip", "compresslevel": 1}
    bins.to_csv(tables / f"{bin_label}_total_old_new_ntr.tsv.gz", sep="\t", index=False, compression=fast_gzip)
    kinetics.to_csv(tables / f"{bin_label}_domain_gene_kinetics.tsv.gz", sep="\t", index=False, compression=fast_gzip)
    pivot, labels = heatmap_data(kinetics, domains)

    fig, ax = plt.subplots(figsize=(7.4, 6.8), facecolor="white")
    spatial_categories(ax, bins, colors, f"Total RNA expression domains — raw {bin_label}", point_size=args.spatial_point_size)
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=colors[x], label=f"D{x}", markersize=6) for x in domains]
    ax.legend(handles=handles, title="Domain", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)
    save(fig, figures, f"01_total_rna_spatial_domains_{bin_label}", args.dpi)

    fig, ax = plt.subplots(figsize=(7.4, 6.2), facecolor="white")
    for domain in domains:
        sub = bins[bins["domain"] == domain]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=args.umap_point_size, alpha=0.86, linewidths=0, color=colors[domain], label=f"D{domain}", rasterized=True)
        ax.text(sub["UMAP1"].median(), sub["UMAP2"].median(), f"D{domain}", ha="center", va="center", fontsize=8, fontweight="bold",
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 1.4})
    ax.set_title(f"Total RNA expression embedding — raw {bin_label}", fontsize=12, fontweight="bold", color=INK)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_aspect("equal", adjustable="datalim")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False, fontsize=8)
    ax.grid(False); ax.spines[["top", "right"]].set_visible(False)
    save(fig, figures, f"02_total_rna_umap_{bin_label}", args.dpi)

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.0), facecolor="white", constrained_layout=True)
    decomposition_titles = (
        ["Total RNA", "Expected background-old component", "Apparent new / T>C-positive component"]
        if control else ["Total RNA", "Expected old RNA", "Expected new RNA"]
    )
    for index, (ax, column, title) in enumerate(zip(axes, ["total", "expected_old", "expected_new"], decomposition_titles)):
        spatial_continuous(ax, bins, column, title, NPG_GENERAL, "log1p expected molecules", transform=np.log1p, point_size=args.spatial_point_size)
        if index:
            ax.set_ylabel("")
    fig.suptitle(
        "Unlabeled-control decomposition quantifies background-positive signal"
        if control else "SLAM-aware decomposition preserves abundance and adds RNA age",
        fontsize=13, fontweight="bold", color=INK,
    )
    save(fig, figures, f"03_total_old_new_spatial_decomposition_{bin_label}", args.dpi)

    fig, ax = plt.subplots(figsize=(7.4, 6.8), facecolor="white")
    spatial_continuous(
        ax, bins, "corrected_ntr",
        f"Apparent NTR in unlabeled control — raw {bin_label}" if control else f"Corrected new-to-total RNA ratio — raw {bin_label}",
        NPG_ENRICHMENT, "Apparent NTR" if control else "Corrected NTR", point_size=args.spatial_point_size,
    )
    save(fig, figures, f"04_corrected_ntr_spatial_{bin_label}", args.dpi)

    fig, axes = plt.subplots(1, 2, figsize=(14.2, max(6.0, 0.19 * len(pivot) + 2.3)), facecolor="white", gridspec_kw={"width_ratios": [1, 1.55]})
    draw_domain_box(axes[0], bins, colors, domains, args.interpretation_mode)
    draw_heatmap(
        axes[1], pivot, labels,
        "Domain-specific background deviation" if control else "Domain-specific RNA turnover",
        "Δ apparent NTR vs gene mean" if control else "ΔNTR vs gene mean",
    )
    save(fig, figures, f"05_domain_ntr_and_dynamic_gene_heatmap_{bin_label}", args.dpi)

    manifest = []
    titles = {
        f"01_total_rna_spatial_domains_{bin_label}": "Total RNA spatial domains",
        f"02_total_rna_umap_{bin_label}": "Total RNA UMAP",
        f"03_total_old_new_spatial_decomposition_{bin_label}": "Total/background/apparent-new spatial decomposition" if control else "Total/old/new spatial decomposition",
        f"04_corrected_ntr_spatial_{bin_label}": "Apparent NTR negative-control spatial map" if control else "Corrected NTR spatial map",
        f"05_domain_ntr_and_dynamic_gene_heatmap_{bin_label}": "Domain background NTR and background-high genes" if control else "Domain NTR and dynamic genes",
    }
    for stem, title in titles.items():
        manifest.append({"stem": stem, "title": title, "png": f"02_figures/{stem}.png", "pdf": f"02_figures/{stem}.pdf", "classification": "expression-derived"})
    pd.DataFrame(manifest).to_csv(metadata / "figure_manifest.tsv", sep="\t", index=False)
    provenance = {
        "sample": args.sample,
        "interpretation_mode": args.interpretation_mode,
        "bin_size_raw": args.bin_size_raw,
        "coordinate_um_per_raw_unit": args.coordinate_um,
        "bin_size_um": args.bin_size_um if args.bin_size_um is not None else args.bin_size_raw * args.coordinate_um,
        "bin_size_pixels": int(round(args.bin_size_raw / args.coordinate_um)) if args.coordinate_um else None,
        "spatial_point_size": args.spatial_point_size,
        "umap_point_size": args.umap_point_size,
        "canonical_orientation": "X increases right; Y increases downward",
        "coordinate_frame": "Celatlas tissue_positions_list.csv full-resolution pixel coordinates",
        "histology_background": False,
        "image_classification": "expression-derived; no histology or microscopy background",
        "source_artifacts": {
            "domains": str(Path(args.domains).resolve()),
            "corrected_bin_gene": str(Path(args.corrected_bin_gene).resolve()),
            "bin_summary": str(Path(args.bin_summary).resolve()) if args.bin_summary else None,
            "domain_gene_kinetics": str(Path(args.domain_gene_kinetics).resolve()),
            "smoothed_bin_gene": str(Path(args.smoothed_bin_gene).resolve()),
        },
        "interpretation": (
            "Unlabeled negative-control figures quantify residual T>C-positive and apparent-NTR background; they do not represent newly synthesized RNA or turnover."
            if control else
            "Figures demonstrate additional RNA-age and turnover information; no matched external method comparison was available."
        ),
        "n_bins": int(len(bins)),
        "total_molecules": int(round(bins["total"].sum())),
        "assignment_totals": {
            "old": float(bins["expected_old"].sum()),
            "new": float(bins["expected_new"].sum()),
        },
        "clustering": {
            "total": {
                "algorithm": "Scanpy normalize_total/log1p/HVG/PCA/neighbors/Leiden/UMAP",
                "coordinate_weight": 0.0,
            }
        },
        "cluster_method_note": "Scanpy normalize_total/log1p/HVG/PCA/neighbors/Leiden；坐标权重为 0",
        "kinetics": {
            "mode": "unlabeled_control_apparent_ntr" if control else "single_sample_posterior_calibrated_ntr",
            "half_life_reported": False,
            "reason": "no matched time course or pulse/chase design",
        },
        "domains": domains,
        "cluster_colors": colors,
    }
    (metadata / "coordinate_and_figure_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
