#!/usr/bin/env python3
"""Build the single customer-facing, offline Spatial SLAM report."""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import os
from pathlib import Path


FIGURE_ORDER = {
    "01_total_rna_spatial_domains_bin50.png": 1,
    "02_total_rna_umap_bin50.png": 2,
    "03_total_old_new_spatial_decomposition_bin50.png": 3,
    "04_corrected_ntr_spatial_bin50.png": 4,
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png": 5,
    "06_representative_dynamic_genes_old_new_ntr_bin50.png": 6,
    "02_slamdunk_conversion_qc.png": 7,
    "03_grand_slam_gene_ntr.png": 8,
    "04_grand_slam_model_qc.png": 9,
    "05_tc_mismatch_enrichment_over_background.png": 10,
}

TITLES = {
    "01_total_rna_spatial_domains_bin50.png": "Total RNA bin50 空间表达域",
    "02_total_rna_umap_bin50.png": "Total RNA bin50 UMAP",
    "03_total_old_new_spatial_decomposition_bin50.png": "Total / old / new RNA 空间分解",
    "04_corrected_ntr_spatial_bin50.png": "校正后新生 RNA 比例（NTR）",
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png": "空间域 RNA 动力学与动态基因",
    "06_representative_dynamic_genes_old_new_ntr_bin50.png": "代表性动态基因：Old / New / NTR 同基因对照",
    "02_slamdunk_conversion_qc.png": "slamdunk：T>C 转换与覆盖 QC",
    "03_grand_slam_gene_ntr.png": "GRAND-SLAM：基因 NTR 后验估计",
    "04_grand_slam_model_qc.png": "GRAND-SLAM：全局 NTR 与转换率模型 QC",
    "05_tc_mismatch_enrichment_over_background.png": "T>C 相对其他突变类型的富集",
}

CAPTIONS = {
    "01_total_rna_spatial_domains_bin50.png": "基于官方 Celatlas bin50 filtered total-RNA 矩阵独立进行 Scanpy 聚类；空间位置直接采用 Celatlas tissue_positions_list.csv 的 full-resolution pixel 坐标，X 向右、Y 向下。",
    "02_total_rna_umap_bin50.png": "UMAP 仅表示 total RNA 表达相似性，不代表物理组织坐标。",
    "03_total_old_new_spatial_decomposition_bin50.png": "在同一 total RNA 空间结构上分解 expected old 与 expected new RNA；old/new 不单独聚类。各面板显示色阶按 p2–p98 裁剪，统计值未修改。",
    "04_corrected_ntr_spatial_bin50.png": "模型校正后的 new-to-total RNA ratio 用于展示组织内相对 T>C-positive/new 信号，不解释为绝对更新速率。显示色阶按 p2–p98 裁剪，统计值未修改。",
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png": "比较 total RNA 空间域的校正 NTR，并展示相对全局偏移的候选动态基因程序；结果为单样本探索性证据。",
    "06_representative_dynamic_genes_old_new_ntr_bin50.png": "每行是同一基因；old/new 使用基因内共享 log1p 色阶，NTR 固定为 0–1，以区分丰度与更新比例。",
    "02_slamdunk_conversion_qc.png": "仅使用表达推断组织区域 reads；展示有足够 read 支持的 UTR 转换率分布、覆盖关系和高 read-count 基因。p98 截断只用于显示，不改变统计值。",
    "03_grand_slam_gene_ntr.png": "展示 GRAND-SLAM 基因 NTR 后验分布、覆盖与不确定性的关系，以及高 read-count 基因的 90% 可信区间。",
    "04_grand_slam_model_qc.png": "展示 GRAND-SLAM 全局外显子/内含子 NTR 及 old/new 转换率模型；完整错配谱在下一图单独呈现。",
    "05_tc_mismatch_enrichment_over_background.png": "在完整 12 类碱基替换谱中突出 T>C，并报告其相对其余 11 类替换中位数和最大值的倍数；slamdunk UTR T>C 率作为独立验证。",
}

CONTROL_TITLES = {
    "03_total_old_new_spatial_decomposition_bin50.png": "Total / 背景 old / 表观 new 空间分解",
    "04_corrected_ntr_spatial_bin50.png": "无标记阴性对照的表观 NTR",
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png": "空间域背景信号与背景偏高基因",
    "06_representative_dynamic_genes_old_new_ntr_bin50.png": "代表性背景偏高基因：组分与表观 NTR",
    "03_grand_slam_gene_ntr.png": "GRAND-SLAM：基因表观 NTR 后验",
    "04_grand_slam_model_qc.png": "GRAND-SLAM：阴性对照背景模型 QC",
}

CONTROL_CAPTIONS = {
    "03_total_old_new_spatial_decomposition_bin50.png": "普通空转未进行代谢标记；expected new 仅表示流程在无标记数据中识别出的 T>C 背景阳性组分，不代表新生 RNA。",
    "04_corrected_ntr_spatial_bin50.png": "无标记阴性对照中的表观 NTR 用于量化流程残余背景与假阳性上限，不能解释为 RNA 更新活性。",
    "05_domain_ntr_and_dynamic_gene_heatmap_bin50.png": "比较不同 total RNA 空间域的表观 NTR 背景；偏高基因是背景诊断对象，不是动态基因。",
    "06_representative_dynamic_genes_old_new_ntr_bin50.png": "每行是同一背景偏高基因；组分与表观 NTR 仅用于检查背景信号的空间异质性。",
    "03_grand_slam_gene_ntr.png": "GRAND-SLAM 在无标记阴性对照上输出的后验作为零假设诊断，不代表真实新生 RNA 比例。",
    "04_grand_slam_model_qc.png": "展示 GRAND-SLAM 在无标记数据上拟合出的背景模型；结果不用于半衰期或转录速率推断。",
}

ACTIVE_BIN_LABEL = "bin50"
ACTIVE_BIN_SIZE_UM = 50.0


def configure_bin(bin_label, bin_size_um):
    """Switch the six spatial assets and visible labels to the requested bin scale."""
    global ACTIVE_BIN_LABEL, ACTIVE_BIN_SIZE_UM
    ACTIVE_BIN_LABEL = str(bin_label)
    ACTIVE_BIN_SIZE_UM = float(bin_size_um)
    if ACTIVE_BIN_LABEL == "bin50":
        return
    for mapping in (FIGURE_ORDER, TITLES, CAPTIONS, CONTROL_TITLES, CONTROL_CAPTIONS):
        replacements = {}
        for key, value in list(mapping.items()):
            if "bin50" not in key:
                continue
            replacements[key.replace("bin50", ACTIVE_BIN_LABEL)] = (
                value.replace("bin50", ACTIVE_BIN_LABEL) if isinstance(value, str) else value
            )
            del mapping[key]
        mapping.update(replacements)

REPORT_SECTIONS = [
    ("section-01", "1 项目概览与核心结论"),
    ("section-02", "2 守恒数据流与关键损失"),
    ("section-03", "3 数据质量与 Celatlas 基线"),
    ("section-04", "4 Total RNA 空间结构"),
    ("section-05", "5 RNA 年龄分解与校正 NTR"),
    ("section-06", "6 空间动力学与代表基因"),
    ("section-07", "7 Tissue-only bulk SLAM"),
    ("section-08", "8 方法、限制与运行信息"),
]

REFERENCES = [
    "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology. 2018;19:15. DOI: 10.1186/s13059-017-1382-0.",
    "McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv. 2018. DOI: 10.48550/arXiv.1802.03426.",
    "Ståhl PL, Salmén F, Vickovic S, et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science. 2016;353:78-82. DOI: 10.1126/science.aaf2403.",
    "Herzog VA, Reichholf B, Neumann T, et al. Thiol-linked alkylation of RNA to assess expression dynamics. Nature Methods. 2017;14:1198-1204. DOI: 10.1038/nmeth.4435.",
    "Neumann T, Herzog VA, Muhar M, et al. Quantification of experimentally induced nucleotide conversions in high-throughput sequencing datasets. BMC Bioinformatics. 2019;20:258. DOI: 10.1186/s12859-019-2849-7.",
    "Jürges C, Dölken L, Erhard F. Dissecting newly transcribed and old RNA using GRAND-SLAM. Bioinformatics. 2018;34:i218-i226. DOI: 10.1093/bioinformatics/bty256.",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bin-label", default="bin50")
    parser.add_argument("--bin-size-um", type=float, default=50.0)
    parser.add_argument("--figure-dir", required=True)
    parser.add_argument("--spatial-metrics", required=True)
    parser.add_argument("--barcode-qc")
    parser.add_argument("--call-qc")
    parser.add_argument("--bulk-summary")
    parser.add_argument("--bulk-figure-dir")
    parser.add_argument("--bulk-spectrum-table")
    parser.add_argument("--celatlas-metrics")
    parser.add_argument("--chip-assignment-qc", help="Project-level full-read chip identity QC JSON")
    parser.add_argument("--sample-metadata", help="Sample-sheet TSV used for project/species/chemistry metadata")
    parser.add_argument("--domain-kinetics", help="Gzipped domain-gene kinetics TSV")
    parser.add_argument("--resource-summary", help="Consolidated pipeline resource TSV")
    parser.add_argument("--software-versions", help="Bulk tool software-version TSV")
    parser.add_argument("--interpretation-mode", choices=("standard", "unlabeled_control"), default="standard")
    parser.add_argument("--analysis-markdown", help="Optional analysis-only Markdown included in the delivery manifest")
    parser.add_argument("--read-flow-section", help="Validated Sankey report fragment without Plotly.js")
    parser.add_argument("--read-flow-plotly", help="Local Plotly.js asset to inline once")
    parser.add_argument("--read-flow-manifest", help="Portable Sankey module manifest")
    parser.add_argument("--footer-logo", help="PNG logo embedded at the end of the offline report")
    parser.add_argument("--output", help="Final report path; defaults to <sample>_report.html")
    parser.add_argument("--manifest")
    return parser.parse_args()


def read_json(path):
    if not path or not Path(path).is_file():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    if not path or not Path(path).is_file():
        return []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_sample_metadata(path, sample):
    return next((row for row in read_tsv(path) if row.get("sample_id") == sample), {})


def fmt_count(value):
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError):
        return "NA"


def fmt_float(value, digits=3):
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "NA"


def fmt_pct_value(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number * 100:.2f}%" if abs(number) <= 1 else f"{number:.2f}%"


def fmt_duration(seconds):
    try:
        seconds = int(round(float(seconds)))
    except (TypeError, ValueError):
        return "NA"
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours}h {minutes}m {seconds}s" if hours else f"{minutes}m {seconds}s"


def apply_scanpy_summary(spatial_metrics, figure_dir):
    summary = read_json(Path(figure_dir) / "scanpy_spatial_umap_summary.json")
    if not summary:
        return spatial_metrics
    updated = dict(spatial_metrics)
    clustering = {}
    for assignment, values in summary.get("assignments", {}).items():
        clustering[assignment] = {
            "algorithm": values.get("algorithm", "scanpy"),
            "scanpy_version": values.get("scanpy_version", "NA"),
            "n_bins": values.get("n_bins", "NA"),
            "n_clusters": values.get("n_clusters", "NA"),
            "n_hvg": values.get("n_hvg", "NA"),
            "n_pcs": values.get("n_pcs", "NA"),
            "neighbors": values.get("neighbors", "NA"),
            "leiden_resolution": values.get("leiden_resolution", "NA"),
        }
    if clustering:
        updated["clustering"] = clustering
        updated["n_bins"] = max(
            [int(v.get("n_bins", 0)) for v in clustering.values() if str(v.get("n_bins", "0")).isdigit()]
            or [updated.get("n_bins", 0)]
        )
    updated["cluster_method_note"] = "Scanpy normalize_total/log1p/HVG/PCA/neighbors/Leiden/UMAP"
    return updated


def apply_kinetics_summary(spatial_metrics, figure_dir):
    summary = read_json(Path(figure_dir) / "kinetics_summary.json")
    if not summary:
        return spatial_metrics
    updated = dict(spatial_metrics)
    updated["kinetics"] = {
        "mode": summary.get("mode", "single_sample_corrected_ntr"),
        "label_hours": summary.get("label_hours", "NA"),
        "half_life_reported": summary.get("half_life_reported", False),
    }
    return updated


def select_top_dynamic_genes(path, per_domain=5):
    grouped = {}
    for row in read_tsv(path):
        try:
            score = float(row.get("dynamic_score", "nan"))
        except ValueError:
            continue
        if score != score:
            continue
        grouped.setdefault(str(row.get("domain", "NA")), []).append((score, row))
    selected = []
    for domain in sorted(grouped, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value)):
        selected.extend(row for _, row in sorted(grouped[domain], key=lambda item: item[0], reverse=True)[:per_domain])
    return selected


def select_resource_rows(rows):
    wanted = {
        "02_celatlas_full": "Celatlas 完整基线",
        f"scanpy_total_{ACTIVE_BIN_LABEL}": f"Total RNA Scanpy {ACTIVE_BIN_LABEL}",
        f"correction_{ACTIVE_BIN_LABEL}": "NTR 校正",
        "slamdunk_all": "slamdunk",
        "grand_slam": "GRAND-SLAM",
        "09": "最终报告",
    }
    selected = []
    for row in rows:
        step = row.get("step", "")
        if step in wanted:
            selected.append({**row, "display_step": wanted[step]})
    return selected


def collect_figures(figure_dirs):
    figures = {}
    for directory in figure_dirs:
        if not directory or not Path(directory).is_dir():
            continue
        for path in Path(directory).glob("*.png"):
            if path.name in FIGURE_ORDER:
                figures[path.name] = path
    return figures


def encode_png(path):
    return base64.b64encode(path.read_bytes()).decode("ascii")


def figure_card(figures, name, layout="full"):
    title = TITLES[name]
    path = figures.get(name)
    if not path:
        return f'<article class="figure-card {layout}"><h3>{html.escape(title)}</h3><div class="empty-state">该可选图件未生成。</div></article>'
    return f"""
    <article class="figure-card {layout}">
      <h3>{html.escape(title)}</h3>
      <img class="expandable" onclick="openFigure(this)" src="data:image/png;base64,{encode_png(path)}" alt="{html.escape(title)}" title="点击放大">
      <p>{html.escape(CAPTIONS[name])}</p>
    </article>"""


def metric_cards(sample, sample_meta, spatial_metrics, interpretation_mode="standard"):
    assignment = spatial_metrics.get("assignment_totals", {})
    total = float(spatial_metrics.get("total_molecules", 0) or 0)
    new = float(assignment.get("new", 0) or 0)
    domains = spatial_metrics.get("domains", [])
    control = interpretation_mode == "unlabeled_control"
    cards = [
        ("样本", sample),
        ("物种", sample_meta.get("species", "Mus musculus").replace("_", " ")),
        ("分析尺度", f"bin{spatial_metrics.get('bin_size_raw', 50)} · {spatial_metrics.get('bin_size_um', 25)} µm"),
        ("有效空间 bins", fmt_count(spatial_metrics.get("n_bins"))),
        ("Total RNA 空间域", fmt_count(len(domains))),
        ("总分子数", fmt_count(total)),
        ("表观 new 背景分子" if control else "Expected new RNA", fmt_count(new)),
        ("表观 NTR" if control else "校正后 new fraction", f"{new / total:.2%}" if total else "NA"),
    ]
    return "".join(
        f'<div class="metric-card"><div class="metric-label">{html.escape(str(label))}</div><div class="metric-value">{html.escape(str(value))}</div></div>'
        for label, value in cards
    )


def chip_assignment_for_sample(payload, sample):
    for row in payload.get("sample_decisions", []):
        if row.get("sample_id") == sample:
            return row
    return {}


def qc_rows(barcode_qc, call_qc, celatlas, chip_assignment):
    official_barcode = summary_mapping(celatlas, "barcode_summary")
    star = summary_mapping(celatlas, "star_summary")
    feature = summary_mapping(celatlas, "featureCounts_summary")
    count = summary_mapping(celatlas, "count_summary")
    total = barcode_qc.get("total_reads")
    kept = barcode_qc.get("kept_reads")
    rows = [
        (
            "输入身份",
            "当前芯片",
            chip_assignment.get("declared_chip_id", "NA"),
            "当前样本声明芯片",
        ),
        (
            "输入身份",
            "当前芯片 exact+PolyT 命中率",
            fmt_pct_value(chip_assignment.get("declared_exact_poly_t_rate")),
            "结构性输入 QC",
        ),
        (
            "输入身份",
            "当前芯片 1-mismatch-corrected+PolyT 命中率",
            fmt_pct_value(chip_assignment.get("declared_corrected_poly_t_rate")),
            "结构性输入 QC",
        ),
        ("SLAM 主链", "Raw reads", fmt_count(total), "输入 FASTQ reads"),
        ("SLAM 主链", "Barcode/PolyT 过滤后 reads", fmt_count(kept), fmt_pct_value((kept / total) if total else None)),
        ("SLAM 主链", "R1 / R2 平均质量", f"{fmt_float(barcode_qc.get('mean_r1_quality'))} / {fmt_float(barcode_qc.get('mean_r2_quality'))}", "Phred"),
        ("SLAM 主链", "进入分子判定的 BAM reads", fmt_count(call_qc.get("total_bam_reads")), "独立分子判定主链"),
        ("SLAM 主链", "Old / new molecules", f"{fmt_count(call_qc.get('molecules_old'))} / {fmt_count(call_qc.get('molecules_new'))}", "离散分子判定"),
        ("Celatlas 1.8", "Valid reads", fmt_count(official_barcode.get("Valid Reads")), fmt_pct_value(official_barcode.get("Valid Reads Fraction"))),
        ("Celatlas 1.8", "唯一比对 reads", fmt_count(star.get("Uniquely Mapped Reads")), fmt_pct_value(star.get("Uniquely Mapped Reads Fraction"))),
        ("Celatlas 1.8", "Exonic reads", fmt_count(feature.get("Reads Assigned To Exonic Regions")), fmt_pct_value(feature.get("Reads Assigned To Exonic Regions Fraction"))),
        ("Celatlas 1.8", "检测基因", fmt_count(count.get("Total Genes")), "官方 count"),
        ("Celatlas 1.8", "测序饱和度", fmt_pct_value(count.get("Saturation")), "官方 denominator"),
    ]
    return rows


def simple_table(headers, rows, classes=""):
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table class="data-table {classes}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def summary_mapping(payload, name):
    summary = payload.get(name, {})
    if isinstance(summary, dict) and "metric_list" in summary:
        return {item.get("name"): item.get("value") for item in summary["metric_list"]}
    return summary


def paginated_table(table_id, headers, rows, row_domains=None):
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body_rows = []
    for index, row in enumerate(rows):
        domain = "" if row_domains is None else str(row_domains[index])
        cells = "".join(f"<td>{html.escape(str(value))}</td>" for value in row)
        body_rows.append(f'<tr data-domain="{html.escape(domain)}">{cells}</tr>')
    return f"""
    <div class="table-wrap"><table id="{table_id}" class="data-table paginated-table pagination-table" data-page-size="10">
      <thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody>
    </table></div>
    <div class="pager" data-pager-for="{table_id}">
      <button type="button" data-page-prev>上一页</button><span data-page-status>1 / 1</span><button type="button" data-page-next>下一页</button>
    </div>"""


def dynamic_gene_table(rows, interpretation_mode="standard"):
    control = interpretation_mode == "unlabeled_control"
    if not rows:
        label = "背景偏高基因" if control else "动态基因"
        return f'<div class="empty-state">未提供 domain-gene kinetics 表；Top {label}表未生成。</div>'
    domains = sorted({str(row.get("domain", "NA")) for row in rows})
    options = "".join(f'<option value="{html.escape(domain)}">Domain {html.escape(domain)}</option>' for domain in domains)
    values = [
        (
            row.get("domain", "NA"), row.get("gene_name") or row.get("gene_id", "NA"), fmt_count(row.get("total")),
            fmt_float(row.get("ntr")), fmt_float(row.get("global_ntr")), fmt_float(row.get("delta_ntr")), fmt_float(row.get("dynamic_score")),
        )
        for row in rows
    ]
    selector = f'<div class="table-tools"><label>空间域 <select data-filter-for="dynamic-gene-table"><option value="all">全部</option>{options}</select></label><span>每个 domain 按 {"background score" if control else "dynamic score"} 展示 Top 5</span></div>'
    return selector + paginated_table(
        "dynamic-gene-table", ["Domain", "基因", "Total", "表观 NTR" if control else "NTR", "全局表观 NTR" if control else "Global NTR", "Δ表观 NTR" if control else "ΔNTR", "Background score" if control else "Dynamic score"], values,
        [row.get("domain", "NA") for row in rows],
    )


def grand_gene_table(bulk_summary):
    rows = bulk_summary.get("grand_slam", {}).get("top_genes_by_readcount", [])
    if not rows:
        return '<div class="empty-state">GRAND-SLAM Top 基因表未生成。</div>'
    values = [
        (row.get("Symbol") or row.get("Gene", "NA"), row.get("Gene", "NA"), fmt_count(row.get("readcount")), fmt_float(row.get("mean_ntr"), 4))
        for row in rows
    ]
    return paginated_table("grand-gene-table", ["基因", "Ensembl ID", "Read count", "Mean NTR"], values)


def mismatch_spectrum_table(bulk_summary):
    rows = bulk_summary.get("tc_mismatch_enrichment", {}).get("complete_12_channel_spectrum", [])
    if not rows:
        return '<div class="empty-state">完整 12 通道突变谱表未生成。</div>'
    values = [
        (
            row.get("region", "NA"), row.get("channel", "NA"), fmt_count(row.get("coverage")),
            fmt_count(row.get("mismatches")), fmt_pct_value(row.get("mismatch_rate")),
            "T>C" if row.get("channel") == "T>C" else "背景通道",
        )
        for row in rows
    ]
    return paginated_table(
        "mismatch-spectrum-table",
        ["区域", "替换通道", "Coverage", "Mismatch count", "Mismatch rate", "类型"],
        values,
    )


def build_report(sample, sample_meta, spatial_metrics, barcode_qc, call_qc, celatlas, chip_qc, bulk_summary, figures, dynamic_rows, resources, software, interpretation_mode="standard", read_flow_section="", plotly_js="", footer_logo=""):
    control = interpretation_mode == "unlabeled_control"
    assignment = spatial_metrics.get("assignment_totals", {})
    total = float(spatial_metrics.get("total_molecules", 0) or 0)
    new = float(assignment.get("new", 0) or 0)
    mask = bulk_summary.get("tissue_mask", {})
    slam = bulk_summary.get("slamdunk", {})
    grand = bulk_summary.get("grand_slam", {})
    cluster = spatial_metrics.get("clustering", {}).get("total", {})
    nav = "".join(f'<a href="#{section_id}">{title}</a>' for section_id, title in REPORT_SECTIONS)
    chip_assignment = chip_assignment_for_sample(chip_qc, sample)
    qc_table = simple_table(
        ["分析链", "指标", "结果", "说明"],
        qc_rows(barcode_qc, call_qc, celatlas, chip_assignment),
    )
    resource_rows = [
        (row.get("display_step", row.get("step", "NA")), fmt_duration(row.get("wall_seconds")), f"{fmt_float(row.get('max_rss_gb'))} GB")
        for row in select_resource_rows(resources)
    ]
    resource_table = simple_table(["步骤", "Wall time", "峰值内存"], resource_rows) if resource_rows else '<div class="empty-state">本次运行未采集统一资源表。</div>'
    version_table = simple_table(["工具", "版本"], [(row.get("tool", "NA"), row.get("version", "NA")) for row in software]) if software else '<div class="empty-state">软件版本表未提供。</div>'
    dynamic_table = dynamic_gene_table(dynamic_rows, interpretation_mode)
    bulk_table = grand_gene_table(bulk_summary)
    spectrum_table = mismatch_spectrum_table(bulk_summary)
    label_hours = spatial_metrics.get("kinetics", {}).get("label_hours")
    no_half_life = label_hours in (None, "", "NA")
    sample_project = sample_meta.get("project_id", "NA")
    chemistry = sample_meta.get("chemistry", summary_mapping(celatlas, "sample_summary").get("Chemistry", "NA"))
    apparent_fraction = new / total if total else 0.0
    report_title = "Spatial SLAM 分析报告" if ACTIVE_BIN_LABEL == "bin50" else f"Spatial SLAM {ACTIVE_BIN_LABEL} 分析报告"
    eyebrow = "Unlabeled spatial-RNA negative control" if control else "Spatial transcriptome · RNA dynamics"
    hero_text = (
        f"以 total RNA {ACTIVE_BIN_LABEL} 空间表达结构为主线，量化普通空转在无代谢标记条件下产生的 T>C 背景、表观 NTR 和模型假阳性。"
        if control else
        f"以 total RNA {ACTIVE_BIN_LABEL} 空间表达结构为主线，整合 old/new RNA 分解、校正 NTR、动态基因以及 tissue-only bulk SLAM 结果。"
    )
    if control:
        findings = (
            f'<div class="finding"><b>Total RNA 定义空间结构</b>{fmt_count(spatial_metrics.get("n_bins"))} 个 bin 被划分为 {fmt_count(len(spatial_metrics.get("domains",[])))} 个表达域，坐标不参与聚类。</div>'
            f'<div class="finding"><b>表观 new 是背景阳性信号</b>在 {fmt_count(total)} 个分子中模型分配 {fmt_count(new)} 个表观 new 组分（{apparent_fraction:.2%}）；它们不代表新生 RNA。</div>'
            '<div class="finding"><b>空间差异用于检查假阳性</b>空间域和基因层面的表观 NTR 偏移用于识别背景异质性，不作转录动力学解释。</div>'
            '<div class="finding"><b>Bulk 模型是零假设诊断</b>slamdunk 与 GRAND-SLAM 输出用于估计无标记背景，不用于半衰期、转录速率或标记效率推断。</div>'
        )
        section04_title = "无标记背景分解与表观 NTR"
        section04_lead = "在同一空间坐标上比较 total、背景 old 组分、表观 new 组分和表观 NTR。"
        section05_title = "空间背景异质性与代表基因"
        section05_lead = "从空间域到基因层面检查无标记背景的异质性；所有偏高结果均是诊断信号。"
        gene_heading = "各空间域 Top 背景偏高基因"
        bulk_lead = "仅使用表达推断组织区域；slamdunk 与 GRAND-SLAM 用于量化无标记条件下的错配和表观 NTR 背景。"
        bulk_notice = "本样本为无标记阴性对照；T>C 富集、GRAND-SLAM 表观 NTR 和区间均不代表新生 RNA。"
        model_limit = "样本未进行代谢标记，所有 hard-new、expected new、表观 NTR 和 GRAND-SLAM NTR 都是背景或假阳性诊断；不报告半衰期、转录速率或标记效率，也不作为非匹配标记样本的定量校准。"
    else:
        findings = (
            f'<div class="finding"><b>Total RNA 定义稳定空间结构</b>{fmt_count(spatial_metrics.get("n_bins"))} 个 bin 被划分为 {fmt_count(len(spatial_metrics.get("domains",[])))} 个表达域，坐标不参与聚类。</div>'
            f'<div class="finding"><b>Old/new 提供探索性信息轴</b>在 {fmt_count(total)} 个分子上估计 {fmt_count(new)} 个 expected new RNA（{apparent_fraction:.2%}），用于描述相对 T&gt;C-positive 信号。</div>'
            '<div class="finding"><b>候选动态信号具有空间异质性</b>代表基因在相同基因内呈现不同的 old、new 与 NTR 空间模式，但不等同于已校准动力学。</div>'
            '<div class="finding"><b>Bulk 结果提供正交证据</b>Tissue-only slamdunk 转换计数与 GRAND-SLAM 后验模型从不同角度描述 T&gt;C 信号。</div>'
        )
        section04_title = "RNA 年龄分解与校正 NTR"
        section04_lead = "在同一空间坐标上比较 total、expected old、expected new 和校正 NTR。"
        section05_title = "空间动力学与代表基因"
        section05_lead = "从空间域汇总到同基因 old/new/NTR 对照，展示转录更新的空间异质性。"
        gene_heading = "各空间域 Top 动态基因"
        bulk_lead = "仅使用表达推断组织区域；slamdunk 给出转换计数，GRAND-SLAM 给出基因后验与全局模型估计。"
        bulk_notice = "T&gt;C 富集倍数为单样本描述性结果，不是生物学重复层面的显著性检验。"
        model_limit = f"{'未提供标记时长，因此不报告半衰期。' if no_half_life else '标记时长已提供。'} 未提供匹配的 no-4sU 阴性对照或经验证的误差模型，因此 Mean NTR 仅作为探索性相对信号解释。完整替换谱来自 slamdunk-filtered tissue BAM。由于没有匹配外部方法和生物学重复，不作定量优越性或重复层面的显著性结论。"
    if ACTIVE_BIN_LABEL != "bin50":
        bulk_notice += f" 本节复用同一样本既有 tissue-only bulk 结果；该 read-level 分析不依赖 {ACTIVE_BIN_LABEL} 空间聚合尺度，未重复运行 SLAM-DUNK 或 GRAND-SLAM。"
    references = "".join(f"<li>{html.escape(reference)}</li>" for reference in REFERENCES)
    footer_logo_html = (
        f'<div class="footer-logo"><img src="data:image/png;base64,{footer_logo}" alt="Epi logo"></div>'
        if footer_logo else ""
    )
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(report_title)}</title>{f'<script src="{plotly_js}"></script>' if plotly_js else ''}
  <style>
    :root {{--bg:#f3f6f7;--panel:#fff;--ink:#1f3038;--muted:#61727b;--line:#dbe3e6;--nav:#20343e;--accent:#c65a3a;--teal:#177f78;--soft:#eef6f4;--warn:#fff6e7;}}
    *{{box-sizing:border-box}} html{{scroll-behavior:smooth}} body{{margin:0;font-family:"Segoe UI","PingFang SC","Microsoft YaHei",Arial,sans-serif;background:var(--bg);color:var(--ink)}}
    .layout{{display:grid;grid-template-columns:270px minmax(0,1fr);min-height:100vh}} .sidebar{{position:sticky;top:0;height:100vh;overflow:auto;background:var(--nav);padding:27px 18px;color:#fff}}
    .brand{{font-weight:760;font-size:1.08rem;margin:0 8px 8px}} .brand-sub{{font-size:.75rem;color:#aebdc4;margin:0 8px 22px}} .sidebar a{{display:block;color:#dce6ea;text-decoration:none;padding:9px 10px;margin:2px 0;border-radius:7px;font-size:.82rem;line-height:1.35;border-left:3px solid transparent}}
    .sidebar a:hover,.sidebar a.active{{background:#304a55;color:#fff;border-left-color:#e68362}} .main{{min-width:0}} .wrap{{max-width:1480px;margin:auto;padding:34px 30px 58px}}
    .report-section{{scroll-margin-top:18px;margin:0 0 28px}} .section-shell{{background:var(--panel);border:1px solid var(--line);border-radius:13px;padding:24px;box-shadow:0 12px 32px rgba(36,58,67,.055)}}
    .hero{{background:linear-gradient(125deg,#20343e 0%,#2b5660 62%,#28776f 100%);color:#fff;border-radius:15px;padding:30px 32px;margin-bottom:20px}} .eyebrow{{font-size:.75rem;letter-spacing:.11em;text-transform:uppercase;color:#b9d6d4;margin-bottom:10px}}
    h1{{font-size:2.15rem;margin:0 0 10px;line-height:1.18}} .hero p{{max-width:980px;color:#e1ebed;line-height:1.7;margin:0}} .meta-line{{display:flex;flex-wrap:wrap;gap:8px;margin-top:17px}} .meta-pill{{border:1px solid rgba(255,255,255,.24);background:rgba(255,255,255,.08);border-radius:99px;padding:6px 10px;font-size:.78rem}}
    .section-heading{{display:flex;gap:13px;align-items:flex-start;margin-bottom:18px}} .section-number{{display:grid;place-items:center;min-width:35px;height:35px;border-radius:9px;background:var(--nav);color:#fff;font-weight:750}} h2{{font-size:1.35rem;margin:4px 0 3px}} .section-lead{{color:var(--muted);line-height:1.55;margin:0;font-size:.91rem}}
    .metric-grid{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin:18px 0}} .metric-card{{border:1px solid var(--line);border-radius:10px;padding:14px 15px;background:#fbfdfd}} .metric-label{{font-size:.78rem;color:var(--muted);margin-bottom:6px}} .metric-value{{font-size:1.18rem;font-weight:760;color:var(--teal);overflow-wrap:anywhere}}
    .finding-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}} .finding{{border-left:4px solid var(--teal);background:var(--soft);border-radius:8px;padding:14px 16px;line-height:1.6}} .finding b{{display:block;margin-bottom:3px}}
    .notice{{border:1px solid #ecd7ae;background:var(--warn);border-radius:9px;padding:13px 15px;color:#665537;line-height:1.58;margin:14px 0}} .notice.info{{background:#eef6f8;border-color:#cfe0e5;color:#3e5964}}
    .figure-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:17px}} .figure-card{{border:1px solid var(--line);border-radius:11px;background:#fff;padding:15px;margin:0 0 17px}} .figure-card.full{{grid-column:1/-1}} .figure-card.paired{{grid-column:auto}}
    .figure-card h3{{font-size:1.03rem;margin:5px 0 12px}} .figure-card img{{display:block;width:100%;height:auto;border:1px solid #e7ecee;background:#fff;cursor:zoom-in}} .figure-card p{{font-size:.88rem;color:var(--muted);line-height:1.6;margin:11px 2px 1px}}
    .subheading{{font-size:1.02rem;margin:22px 0 10px}} .table-wrap{{overflow:auto;border:1px solid var(--line);border-radius:9px}} .data-table{{border-collapse:collapse;width:100%;font-size:.82rem;background:#fff}} .data-table th{{position:sticky;top:0;background:#edf3f4;color:#334b55;text-align:left;padding:9px 10px;border-bottom:1px solid var(--line);white-space:nowrap}} .data-table td{{padding:8px 10px;border-bottom:1px solid #edf1f2;white-space:nowrap}} .data-table tr:nth-child(even) td{{background:#fafcfc}}
    .table-tools,.pager{{display:flex;align-items:center;justify-content:space-between;gap:12px;color:var(--muted);font-size:.8rem;margin:8px 0}} select,button{{font:inherit;border:1px solid #cbd7db;background:#fff;color:var(--ink);border-radius:6px;padding:6px 9px}} button:disabled{{opacity:.42}} .pager{{justify-content:flex-end}}
    details{{border:1px solid var(--line);border-radius:9px;background:#fbfdfd;margin:9px 0}} summary{{cursor:pointer;font-weight:680;padding:12px 14px}} .details-body{{padding:0 14px 14px;color:var(--muted);font-size:.88rem;line-height:1.65}} .method-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}}
    .empty-state{{border:1px dashed #c9d5d9;border-radius:8px;padding:18px;text-align:center;color:var(--muted);background:#fafcfc}} .reference-list{{margin:0;padding-left:22px}} .reference-list li{{margin:0 0 9px}} .footer{{text-align:center;color:var(--muted);font-size:.78rem;border-top:1px solid var(--line);padding-top:18px;margin-top:30px}} .footer-logo{{display:flex;justify-content:center;padding:24px 0 8px}} .footer-logo img{{display:block;width:min(520px,72vw);height:auto}}
    .modal{{display:none;position:fixed;z-index:1000;inset:0;background:rgba(13,24,29,.94);padding:27px;align-items:center;justify-content:center}} .modal.open{{display:flex}} .modal img{{max-width:95vw;max-height:91vh;object-fit:contain;background:#fff}} .modal-close{{position:absolute;right:22px;top:18px;border-color:#6f858e;background:#263b45;color:#fff;font-size:1rem}}
    @media(max-width:980px){{.metric-grid{{grid-template-columns:repeat(2,1fr)}}.figure-grid,.finding-grid,.method-grid{{grid-template-columns:1fr}}.figure-card.full{{grid-column:auto}}}}
    @media(max-width:720px){{.layout{{display:block}}.sidebar{{position:relative;height:auto}}.wrap{{padding:20px 13px 40px}}.metric-grid{{grid-template-columns:1fr 1fr}}.section-shell{{padding:17px}}h1{{font-size:1.65rem}}}}
    @media print{{.sidebar,.modal,.pager,.table-tools{{display:none!important}}.layout{{display:block}}.wrap{{max-width:none;padding:0}}.report-section,.figure-card{{break-inside:avoid}}.section-shell{{box-shadow:none}}.paginated-table tr[hidden]{{display:table-row!important}}details .details-body{{display:block!important}}}}
  </style>
</head>
<body><div class="layout">
<nav class="sidebar"><div class="brand">Celatlas Spatial SLAM</div><div class="brand-sub">Offline Analysis Report</div>{nav}</nav>
<main class="main"><div class="wrap">
<section class="report-section" id="section-01">
  <div class="hero"><div class="eyebrow">{html.escape(eyebrow)}</div><h1>{html.escape(report_title)}</h1><p>{html.escape(hero_text)}</p><div class="meta-line"><span class="meta-pill">项目 {html.escape(sample_project)}</span><span class="meta-pill">物种 {html.escape(sample_meta.get('species','Mus_musculus').replace('_',' '))}</span><span class="meta-pill">Chemistry {html.escape(str(chemistry))}</span><span class="meta-pill">X 向右 · Y 向下</span></div></div>
  <div class="section-shell"><div class="section-heading"><div class="section-number">1</div><div><h2>项目概览与核心结论</h2><p class="section-lead">先给出本次样本最重要的可展示结论，再在后续章节提供技术证据。</p></div></div>
  <div class="metric-grid">{metric_cards(sample, sample_meta, spatial_metrics, interpretation_mode)}</div>
  <div class="finding-grid">{findings}</div></div>
</section>

<section class="report-section" id="section-02"><div class="section-shell"><div class="section-heading"><div class="section-number">2</div><div><h2>守恒数据流与关键损失</h2><p class="section-lead">按 Celatlas 主链展示 FASTQ read pairs 从 barcode 识别、剪切、比对到 feature assignment 的守恒流向。</p></div></div>{read_flow_section}<div class="notice info">各分支沿用 Celatlas 官方指标及其 denominator，不与分子或 bin 数量混算。</div></div></section>

<section class="report-section" id="section-03"><div class="section-shell"><div class="section-heading"><div class="section-number">3</div><div><h2>数据质量与 Celatlas 基线</h2><p class="section-lead">按原始数据、barcode、比对和分子支持顺序汇总关键 QC。</p></div></div>{qc_table}<div class="notice info">Celatlas 官方链与 Spatial SLAM 分子主链采用不同过滤规则和 denominator；表中结果分别报告，不直接混算比例。</div></div></section>

<section class="report-section" id="section-04"><div class="section-shell"><div class="section-heading"><div class="section-number">4</div><div><h2>Total RNA 空间结构</h2><p class="section-lead">total RNA 是后续动力学比较的统一空间结构主线；本报告尺度为 {ACTIVE_BIN_LABEL}（{ACTIVE_BIN_SIZE_UM:g} µm）。</p></div></div><div class="figure-grid">{figure_card(figures,f'01_total_rna_spatial_domains_{ACTIVE_BIN_LABEL}.png','paired')}{figure_card(figures,f'02_total_rna_umap_{ACTIVE_BIN_LABEL}.png','paired')}</div></div></section>

<section class="report-section" id="section-05"><div class="section-shell"><div class="section-heading"><div class="section-number">5</div><div><h2>{section04_title}</h2><p class="section-lead">{section04_lead}</p></div></div><div class="figure-grid">{figure_card(figures,f'03_total_old_new_spatial_decomposition_{ACTIVE_BIN_LABEL}.png','paired')}{figure_card(figures,f'04_corrected_ntr_spatial_{ACTIVE_BIN_LABEL}.png','paired')}</div></div></section>

<section class="report-section" id="section-06"><div class="section-shell"><div class="section-heading"><div class="section-number">6</div><div><h2>{section05_title}</h2><p class="section-lead">{section05_lead}</p></div></div><div class="figure-grid">{figure_card(figures,f'05_domain_ntr_and_dynamic_gene_heatmap_{ACTIVE_BIN_LABEL}.png','full')}{figure_card(figures,f'06_representative_dynamic_genes_old_new_ntr_{ACTIVE_BIN_LABEL}.png','full')}</div><h3 class="subheading">{gene_heading}</h3>{dynamic_table}</div></section>

<section class="report-section" id="section-07"><div class="section-shell"><div class="section-heading"><div class="section-number">7</div><div><h2>Tissue-only bulk SLAM</h2><p class="section-lead">{bulk_lead}</p></div></div><div class="metric-grid"><div class="metric-card"><div class="metric-label">Tissue barcodes</div><div class="metric-value">{fmt_count(mask.get('n_tissue_barcodes',mask.get('n_tissue_spots')))}</div></div><div class="metric-card"><div class="metric-label">Tissue molecules</div><div class="metric-value">{fmt_count(mask.get('tissue_total_molecules'))}</div></div><div class="metric-card"><div class="metric-label">slamdunk counted reads</div><div class="metric-value">{fmt_count(slam.get('total_counted_reads'))}</div></div><div class="metric-card"><div class="metric-label">GRAND-SLAM genes</div><div class="metric-value">{fmt_count(grand.get('n_genes'))}</div></div></div><h3 class="subheading">slamdunk 转换计数 QC</h3>{figure_card(figures,'02_slamdunk_conversion_qc.png','full')}<h3 class="subheading">GRAND-SLAM 基因后验与模型 QC</h3>{figure_card(figures,'03_grand_slam_gene_ntr.png','full')}{figure_card(figures,'04_grand_slam_model_qc.png','full')}<h3 class="subheading">T&gt;C 相对其他突变类型是否突出</h3>{figure_card(figures,'05_tc_mismatch_enrichment_over_background.png','full')}<h3 class="subheading">完整 12 通道突变谱</h3>{spectrum_table}<h3 class="subheading">GRAND-SLAM 高 read-count 基因</h3>{bulk_table}<div class="notice">{bulk_notice}</div></div></section>

<section class="report-section" id="section-08"><div class="section-shell"><div class="section-heading"><div class="section-number">8</div><div><h2>方法、限制与运行信息</h2><p class="section-lead">技术细节默认折叠，便于复核而不干扰主结果阅读。</p></div></div><div class="method-grid"><details><summary>空间分析方法</summary><div class="details-body">尺度：bin{html.escape(str(spatial_metrics.get('bin_size_raw',50)))}（{html.escape(str(spatial_metrics.get('bin_size_um',25)))} µm）。聚类：{html.escape(str(spatial_metrics.get('cluster_method_note','Scanpy')))}。Total RNA 空间域数：{fmt_count(len(spatial_metrics.get('domains',[])))}。物理坐标只用于展示。</div></details><details><summary>图像与 tissue mask 来源</summary><div class="details-body">空间背景分类：{html.escape(str(spatial_metrics.get('image_classification','expression-derived')))}。Bulk mask：{html.escape(str(mask.get('method','NA')))}；保留 {fmt_pct_value(mask.get('tissue_barcode_fraction'))} barcodes。Tissue mask 来自 RNA 密度，不是组织学染色；所有空间图统一 X 向右、Y 向下。</div></details><details><summary>模型解释限制</summary><div class="details-body">{model_limit} 报告不将 Celatlas 与 SLAM 主链 denominator 混算，也不将 old/new 组分单独聚类。</div></details><details><summary>运行资源</summary><div class="details-body">{resource_table}</div></details><details><summary>软件版本</summary><div class="details-body">{version_table}</div></details><details><summary>参考文献</summary><div class="details-body"><ol class="reference-list">{references}</ol></div></details></div></div></section>

<div class="footer">Celatlas Spatial SLAM · Bioinformatics Analysis · Self-contained offline report</div>
{footer_logo_html}
</div></main></div>
<div class="modal" id="figure-modal" onclick="closeFigure(event)"><button class="modal-close" type="button" onclick="closeFigure(event)">关闭 ×</button><img id="modal-image" alt="放大图"></div>
<script>
  function openFigure(img){{const modal=document.getElementById('figure-modal');document.getElementById('modal-image').src=img.src;modal.classList.add('open');}}
  function closeFigure(event){{if(event&&event.target&&event.target.id==='modal-image')return;document.getElementById('figure-modal').classList.remove('open');}}
  document.addEventListener('keydown',event=>{{if(event.key==='Escape')closeFigure();}});
  const tableState={{}};
  function visibleRows(table){{const filter=document.querySelector(`[data-filter-for="${{table.id}}"]`);const value=filter?filter.value:'all';return [...table.tBodies[0].rows].filter(row=>value==='all'||row.dataset.domain===value);}}
  function renderPage(table){{const state=tableState[table.id]||{{page:0}};const rows=[...table.tBodies[0].rows];const shown=visibleRows(table);const size=Number(table.dataset.pageSize||10);const pages=Math.max(1,Math.ceil(shown.length/size));state.page=Math.min(state.page,pages-1);rows.forEach(row=>row.hidden=true);shown.slice(state.page*size,(state.page+1)*size).forEach(row=>row.hidden=false);tableState[table.id]=state;const pager=document.querySelector(`[data-pager-for="${{table.id}}"]`);if(pager){{pager.querySelector('[data-page-status]').textContent=`${{state.page+1}} / ${{pages}} · ${{shown.length}} 行`;pager.querySelector('[data-page-prev]').disabled=state.page===0;pager.querySelector('[data-page-next]').disabled=state.page>=pages-1;}}}}
  document.querySelectorAll('.paginated-table').forEach(table=>{{tableState[table.id]={{page:0}};const pager=document.querySelector(`[data-pager-for="${{table.id}}"]`);pager.querySelector('[data-page-prev]').onclick=()=>{{tableState[table.id].page--;renderPage(table)}};pager.querySelector('[data-page-next]').onclick=()=>{{tableState[table.id].page++;renderPage(table)}};const filter=document.querySelector(`[data-filter-for="${{table.id}}"]`);if(filter)filter.onchange=()=>{{tableState[table.id].page=0;renderPage(table)}};renderPage(table);}});
  const navLinks=[...document.querySelectorAll('.sidebar a')];const observer=new IntersectionObserver(entries=>entries.forEach(entry=>{{if(entry.isIntersecting){{navLinks.forEach(link=>link.classList.toggle('active',link.getAttribute('href')==='#'+entry.target.id));}}}}),{{rootMargin:'-20% 0px -65% 0px'}});document.querySelectorAll('.report-section').forEach(section=>observer.observe(section));
</script></body></html>"""


def write_manifest(path, report, figure_dirs, extra_files=None, read_flow_manifest=None):
    manifest = Path(path)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest_base = manifest.parent.resolve()

    def delivery_row(category, description, target, status="up_to_date"):
        resolved = Path(target).resolve()
        portable = os.path.relpath(resolved, manifest_base)
        return {
            "category": category,
            "description": description,
            "source": portable,
            "target": portable,
            "bytes": resolved.stat().st_size,
            "status": status,
        }

    rows = []
    late_rows = []
    for extra in extra_files or []:
        extra_path = Path(extra)
        if extra_path.is_file() and extra_path.stat().st_size:
            if extra_path.suffix.lower() == ".md":
                category, description = "09_report", "Analysis-only Chinese Markdown summary"
            else:
                category, description = "01_tables", "Complete 12-channel tissue mismatch spectrum"
            target_rows = late_rows if category == "09_report" else rows
            target_rows.append(delivery_row(category, description, extra_path))
    figure_files = {}
    for directory in figure_dirs:
        if not directory or not Path(directory).is_dir():
            continue
        for figure in Path(directory).iterdir():
            png_name = figure.stem + ".png"
            if figure.suffix.lower() not in {".png", ".pdf"} or png_name not in FIGURE_ORDER:
                continue
            figure_files[(png_name, figure.suffix.lower())] = figure
    for png_name in sorted(FIGURE_ORDER, key=FIGURE_ORDER.get):
        for suffix in (".png", ".pdf"):
            figure = figure_files.get((png_name, suffix))
            if figure is None:
                continue
            rows.append(delivery_row("02_figures", TITLES[png_name], figure))
    if read_flow_manifest:
        module_manifest = Path(read_flow_manifest)
        for row in read_tsv(module_manifest):
            target = (module_manifest.parent / row["target"]).resolve()
            rows.append(delivery_row("02_read_flow", row["description"], target, row["status"]))
    rows.extend(late_rows)
    rows.append(delivery_row("09_report", "Final offline Spatial SLAM report", report))
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["category", "description", "source", "target", "bytes", "status"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    args = parse_args()
    configure_bin(args.bin_label, args.bin_size_um)
    if args.interpretation_mode == "unlabeled_control":
        TITLES.update(CONTROL_TITLES)
        CAPTIONS.update(CONTROL_CAPTIONS)
        REPORT_SECTIONS[4] = ("section-05", "5 无标记背景分解与表观 NTR")
        REPORT_SECTIONS[5] = ("section-06", "6 空间背景异质性与代表基因")
    figure_dir = Path(args.figure_dir)
    spatial_metrics = apply_kinetics_summary(apply_scanpy_summary(read_json(args.spatial_metrics), figure_dir), figure_dir)
    barcode_qc = read_json(args.barcode_qc)
    call_qc = read_json(args.call_qc)
    bulk_summary = read_json(args.bulk_summary)
    celatlas = read_json(args.celatlas_metrics)
    chip_qc = read_json(args.chip_assignment_qc)
    sample_meta = read_sample_metadata(args.sample_metadata, args.sample)
    dynamic_rows = select_top_dynamic_genes(args.domain_kinetics)
    resources = read_tsv(args.resource_summary)
    software = read_tsv(args.software_versions)
    figure_dirs = [figure_dir] + ([Path(args.bulk_figure_dir)] if args.bulk_figure_dir else [])
    figures = collect_figures(figure_dirs)
    read_flow_section = Path(args.read_flow_section).read_text(encoding="utf-8") if args.read_flow_section else '<div class="empty-state">守恒数据流模块未提供。</div>'
    plotly_js = (
        "data:text/javascript;base64," + base64.b64encode(Path(args.read_flow_plotly).read_bytes()).decode("ascii")
        if args.read_flow_plotly else ""
    )
    footer_logo_path = Path(args.footer_logo) if args.footer_logo else None
    if footer_logo_path is not None and (not footer_logo_path.is_file() or footer_logo_path.stat().st_size == 0):
        raise FileNotFoundError(f"footer logo is missing or empty: {footer_logo_path}")
    footer_logo = encode_png(footer_logo_path) if footer_logo_path else ""
    output = Path(args.output) if args.output else Path(f"{args.sample}_report.html")
    if not output.is_absolute():
        output = figure_dir / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        build_report(
            args.sample,
            sample_meta,
            spatial_metrics,
            barcode_qc,
            call_qc,
            celatlas,
            chip_qc,
            bulk_summary,
            figures,
            dynamic_rows,
            resources,
            software,
            interpretation_mode=args.interpretation_mode,
            read_flow_section=read_flow_section,
            plotly_js=plotly_js,
            footer_logo=footer_logo,
        ),
        encoding="utf-8",
    )
    if args.manifest:
        extras = [value for value in (args.bulk_spectrum_table, args.analysis_markdown) if value]
        write_manifest(args.manifest, output, figure_dirs, extras, args.read_flow_manifest)
    print(output)


if __name__ == "__main__":
    main()
