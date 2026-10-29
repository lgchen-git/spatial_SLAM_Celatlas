#!/usr/bin/env python3
"""Build the project analysis catalog as Markdown and a self-contained offline HTML report.

Both artifacts are generated from the same structured catalog data so the HTML and
Markdown cannot drift. The HTML follows the Bio Report v4 contract: inline CSS/JS,
sticky ordered sidebar, static bounded tables, notices, collapsed methods and print CSS.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path


PROJECT = "Celatlas Spatial SLAM"
SAMPLES = ("example_sample_a", "example_sample_b")
STAGE_RANGE = "01–09"

BASIC = [
    ("B1", "输入质控与芯片身份交叉验证", "01",
     "完整 R1 双芯片白名单交叉匹配（exact + 1 mismatch + PolyT），声明芯片 exact PolyT 命中率需 ≥2 倍胜出",
     "run_chip_assignment_qc.py、chip_assignment_scan.rs → chip_assignment_qc.json"),
    ("B2", "官方 Celatlas 1.8 基线", "02",
     "复用并校验官方 total RNA bin50 filtered 矩阵与自动分割结果",
     "validate_reused_celatlas.py"),
    ("B3", "空间 barcode / UMI 提取", "03",
     "按 C4L15C4L15C4U10T18 模式提取空间 barcode 与 UMI",
     "extract_spatial_barcodes.py、run_parallel_barcode.py"),
    ("B4", "修剪与 STAR 比对", "04",
     "序列修剪、STAR 比对，保留 MD/NH 标签供转换感知计数",
     "pipeline/Snakefile（trim_align rule）"),
    ("B5", "基因归属与 old/new 分子判定", "05",
     "外显子/内含子注释构建、T>C 分子级判定（new_threshold=1）",
     "build_exon_intron_annotation.py、call_new_old_molecules.py"),
    ("B6", "空间表达矩阵构建", "06",
     "构建 total / old / new 原始空间表达矩阵",
     "build_old_new_matrices.py"),
    ("B7", "Total RNA 空间域与 UMAP", "07",
     "官方 bin50（50 µm / 100 pixels）独立 Leiden 聚类（HVG 3,000、25 PCs、15 neighbors、resolution 1.0）",
     "spatial_slam_scanpy_analysis.py → 空间域与 UMAP 输出"),
    ("B8", "组织区域识别与 TissUUmaps", "08",
     "表达派生 tissue mask（自动）或人工 GeoJSON ROI 覆盖；TissUUmaps 交互导出",
     "build_tissue_mask.py、build_expression_tissue_mask.py、export_tissuumaps.py"),
    ("B9", "整合报告与交付校验", "09",
     "离线 HTML 报告（7 章）、交付 manifest 与必需章节校验",
     "build_spatial_slam_combined_report.py、validate_delivery_report.py"),
]

ADVANCED = [
    ("A1", "Expected old/new 空间丰度分解", "07",
     "基于两组分 binomial mixture 校正模型，在同一官方 bin 网格上分解 expected old / expected new",
     "spatial_slam_correct_new_fraction.py → old/new 空间分解输出"),
    ("A2", "校正 bin 级 NTR 图谱", "07",
     "EM 校正后的 bin 级 NTR 空间分布",
     "spatial_slam_correct_new_fraction.py → 校正 NTR 空间输出"),
    ("A3", "空间域 NTR 分布与动态基因 ΔNTR 热图", "07",
     "domain 级 NTR 汇总、动态基因筛选与 ΔNTR 热图（min_gene_total_umis=100）",
     "spatial_slam_kinetics.py → 空间域动力学输出"),
    ("A4", "代表性动态基因 old/new/NTR 空间图谱", "07",
     "同基因 old RNA、new RNA、NTR 三面板空间图",
     "plot_representative_gene_old_new_ntr.py → 代表基因空间输出"),
    ("A5", "Tissue-only bulk slamdunk 转换 QC", "08",
     "仅组织 barcode 的 FASTQ，conversion-aware tcount、碱基质量与 SNP 过滤",
     "filter_tissue_fastq.py、run_tissue_bulk_tools.sh → slamdunk QC 输出"),
    ("A6", "GRAND-SLAM 基因级 NTR 后验与模型 QC", "08",
     "基因级 NTR 后验、read-count 支持与模型诊断",
     "run_grand_slam.sh、summarize_bulk_slam.py → GRAND-SLAM 后验与模型 QC"),
    ("A7", "T>C 富集与 12 通道突变谱", "08/09",
     "T>C 相对其他突变类型是否突出、完整突变谱背景诊断",
     "summarize_bulk_slam.py → T>C 富集与完整谱表"),
    ("A8", "无标记背景/模型收敛诊断（探索性）", "07/08",
     "无 matched no-4sU 对照时的背景估计、EM 收敛诊断；不用于绝对半衰期",
     "new_fraction_correction_summary.json、grand_slam 模型 QC"),
]

FIGURES = [
    ("01", "total RNA 空间域", "07_spatial_dynamics", "01_total_rna_spatial_domains_bin50", "基础"),
    ("02", "total RNA UMAP", "07_spatial_dynamics", "02_total_rna_umap_bin50", "基础"),
    ("03", "total / old / new 空间分解", "07_spatial_dynamics", "03_total_old_new_spatial_decomposition_bin50", "高级"),
    ("04", "校正 NTR 空间图谱", "07_spatial_dynamics", "04_corrected_ntr_spatial_bin50", "高级"),
    ("05", "空间域 NTR 与动态基因 ΔNTR 热图", "07_spatial_dynamics", "05_domain_ntr_and_dynamic_gene_heatmap_bin50", "高级"),
    ("06", "代表性动态基因 old / new / NTR", "07_spatial_dynamics", "06_representative_dynamic_genes_old_new_ntr_bin50", "高级"),
    ("02", "slamdunk 转换 QC", "08_tissue_bulk", "02_slamdunk_conversion_qc", "高级"),
    ("03", "GRAND-SLAM 基因 NTR 后验", "08_tissue_bulk", "03_grand_slam_gene_ntr", "高级"),
    ("04", "GRAND-SLAM 模型 QC", "08_tissue_bulk", "04_grand_slam_model_qc", "高级"),
    ("05", "T>C 突变相对背景富集", "08_tissue_bulk", "05_tc_mismatch_enrichment_over_background", "高级"),
]

PARAMETERS = [
    ("空间尺度", "官方 Celatlas bin50；50 µm / 100 raw pixels"),
    ("聚类设置", "Seurat HVG 3,000、25 PCs、15 neighbors、Leiden resolution 1.0、random state 0、2 iterations"),
    ("分子判定", "barcode C4L15C4L15C4U10T18；new_threshold=1；min_bin_umis=20"),
    ("slamdunk", "min base quality 27；min SNP coverage 10；min SNP fraction 0.8"),
    ("GRAND-SLAM", "error_mode=exploratory_sample_background；require MD/NH tags"),
]

PROVENANCE = [
    ("样本", "example_sample_a（IAA）、example_sample_b（PEG），独立运行、不跨样本合并"),
    ("结果根目录", "outputs/"),
    ("配置", "config/analysis.yaml、config/samples.tsv、config/libraries.tsv"),
    ("复现入口", "snakemake --snakefile pipeline/Snakefile --configfile config/<run>.yaml"),
    ("报告生成器", "pipeline/scripts/09_report/build_analysis_catalog_report.py"),
]

CSV = ","


def esc(value: str) -> str:
    return html.escape(str(value), quote=True)


def md_table(rows, headers):
    lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return "\n".join(lines)


def render_markdown() -> str:
    basic_rows = [(code, name, stage, note, scripts) for code, name, stage, note, scripts in BASIC]
    advanced_rows = [(code, name, stage, note, scripts) for code, name, stage, note, scripts in ADVANCED]
    return f"""# 当前项目分析条目清单

> 项目：{PROJECT}（样本 {CSV.join(SAMPLES)}）
> 流程：`01_input_qc` → `09_report`，参数统一由 `config/analysis.yaml` 管理。
> 归类原则：**基础分析** = 从原始数据到标准空间结构的通用空间转录组分析；
> **高级分析** = 依赖 SLAM 代谢标记的定量、动力学与生物学解释。

## 基础分析

{md_table(basic_rows, ["#", "分析条目", "阶段", "内容说明", "关键脚本 / 输出"])}

## 高级分析

{md_table(advanced_rows, ["#", "分析条目", "阶段", "内容说明", "关键脚本 / 输出"])}

## 归类与解释边界

- 基础分析（B1–B9）覆盖“原始数据 → 标准空间结构”的全链路，任何空间转录组项目都通用，可直接作为空间结果交付。
- 高级分析（A1–A8）全部依赖 SLAM T>C 标记，属于定量/动力学层；其中 A8 为探索性诊断，不输出绝对转录速率或半衰期。
- old/new RNA 只用于丰度分解、NTR 校正与 turnover，不单独聚类；空间域编号为各样本独立 Leiden 结果，跨样本需用 marker/表达程序重新对齐。
- tissue mask 是 RNA 表达密度推断区域，不是组织学染色边界；空间图统一 X 向右、Y 向下，UMAP 明确为非物理坐标。
"""


def catalog_table(rows, headers) -> str:
    cells = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        tds = "".join(f"<td>{esc(cell)}</td>" for cell in row)
        body.append(f"<tr>{tds}</tr>")
    return (
        f'<div class="table-wrap"><table><thead><tr>{cells}</tr></thead>'
        f"<tbody>{''.join(body)}</tbody></table></div>"
    )


def key_value_table(rows) -> str:
    body = []
    for key, value in rows:
        body.append(f'<tr><th>{esc(key)}</th><td>{esc(value)}</td></tr>')
    return f'<div class="table-wrap"><table class="kv"><tbody>{"".join(body)}</tbody></table></div>'


def render_html() -> str:
    basic_table = catalog_table(
        [(code, name, stage, note, scripts) for code, name, stage, note, scripts in BASIC],
        ["#", "分析条目", "阶段", "内容说明", "关键脚本 / 输出"],
    )
    advanced_table = catalog_table(
        [(code, name, stage, note, scripts) for code, name, stage, note, scripts in ADVANCED],
        ["#", "分析条目", "阶段", "内容说明", "关键脚本 / 输出"],
    )
    figure_rows = [(stage, name, cls) for _num, stage, name, _stem, cls in FIGURES]
    figure_table = catalog_table(figure_rows, ["阶段目录", "文件名", "归类"])
    parameter_table = key_value_table(PARAMETERS)
    provenance_table = key_value_table(PROVENANCE)

    metric_cards = (
        f'<div class="metric-card"><div class="metric-value">2</div><div class="metric-label">独立样本</div>'
        f'<div class="metric-note">{esc(CSV.join(SAMPLES))}</div></div>'
        f'<div class="metric-card"><div class="metric-value">9</div><div class="metric-label">流程阶段</div>'
        f'<div class="metric-note">{esc(STAGE_RANGE)} 全链路</div></div>'
        f'<div class="metric-card"><div class="metric-value">9</div><div class="metric-label">基础分析条目</div>'
        f'<div class="metric-note">B1–B9 标准空间结构</div></div>'
        f'<div class="metric-card"><div class="metric-value">8</div><div class="metric-label">高级分析条目</div>'
        f'<div class="metric-note">A1–A8 SLAM 动力学</div></div>'
    )

    finding_cards = (
        '<div class="finding-card"><h3>芯片身份交叉验证通过</h3>'
        "<p>项目级 exact PolyT 对角/非对角比为 6.01 倍；A1 与 B1 分别 5.81、6.11 倍，支持当前声明映射。"
        "证据见基础分析 B1。</p></div>"
        '<div class="finding-card"><h3>空间结构以官方 Celatlas bin50 为统一骨架</h3>'
        "<p>Total RNA 空间域（50 µm / 100 pixels）作为唯一聚类主线，old/new 与 SLAM 分子映射到同一官方网格，"
        "避免低深度自建矩阵造成的结构破碎。证据见 B7 与 A1–A4。</p></div>"
        '<div class="finding-card"><h3>SLAM 动力学为单样本探索性结果</h3>'
        "<p>两个样本均无 matched no-4sU 对照，校正 EM 未报告完全收敛；NTR 可展示但不出绝对转录速率或半衰期。"
        "证据见 A2–A8 与解释边界。</p></div>"
    )

    css = """
:root{--bg:#f3f6f7;--panel:#fff;--ink:#1f3038;--muted:#61727b;--line:#dbe3e6;--nav:#20343e;--accent:#c65a3a;--teal:#177f78;--soft:#eef6f4;--warn:#fff6e7;}
*{box-sizing:border-box}body{margin:0;font-family:"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;color:var(--ink);background:var(--bg)}
.layout{display:grid;grid-template-columns:250px 1fr;min-height:100vh}
.sidebar{position:sticky;top:0;height:100vh;background:var(--nav);padding:28px 20px;color:#fff;overflow-y:auto}
.sidebar h3{font-size:.82rem;letter-spacing:.08em;text-transform:uppercase;color:#a9c0c7;margin:0 0 18px}
.sidebar a{display:block;color:#dce6ea;text-decoration:none;padding:9px 12px;border-radius:8px;font-size:.92rem}
.sidebar a:hover{background:rgba(255,255,255,.08)}
.sidebar a.active{background:rgba(255,255,255,.14);color:#fff}
.wrap{padding:34px 40px 48px;max-width:1380px}
.hero{background:linear-gradient(125deg,#20343e,#28776f);color:#fff;border-radius:13px;padding:26px 28px;margin-bottom:22px}
.hero h1{margin:0 0 10px;font-size:1.7rem}.hero p{margin:4px 0;color:#dcebea;line-height:1.6}
.report-section{scroll-margin-top:18px;margin:0 0 28px}
.section-shell{background:var(--panel);border:1px solid var(--line);border-radius:13px;padding:24px;box-shadow:0 12px 32px rgba(36,58,67,.055)}
.section-heading{display:flex;gap:13px;align-items:flex-start;margin-bottom:18px}
.section-number{display:grid;place-items:center;min-width:35px;height:35px;border-radius:9px;background:var(--nav);color:#fff;font-weight:750}
h2{font-size:1.35rem;margin:4px 0 3px}.section-lead{color:var(--muted);line-height:1.55;margin:0;font-size:.91rem}
.metric-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:18px 0}
.metric-card{background:var(--soft);border:1px solid #cfe3de;border-radius:11px;padding:16px}
.metric-value{font-size:1.65rem;font-weight:800;color:var(--teal)}
.metric-label{font-weight:650;margin:2px 0 4px}.metric-note{color:var(--muted);font-size:.8rem;line-height:1.45}
.finding-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:18px 0}
.finding-card{background:#fff;border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:10px;padding:14px 16px}
.finding-card h3{margin:0 0 7px;font-size:1rem}
.finding-card p{margin:0;font-size:.88rem;line-height:1.65;color:#3c4c54}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:.86rem}
th,td{padding:10px 12px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:#edf3f4;white-space:nowrap}
td{white-space:normal}tr:last-child td{border-bottom:0}
table.kv th{width:170px;background:#f6f9f9}
.notice{border-left:4px solid var(--accent);background:var(--warn);padding:14px 16px;border-radius:0 10px 10px 0;margin:18px 0;font-size:.9rem;line-height:1.65}
.notice.info{border-left-color:var(--teal);background:var(--soft)}
details{background:#f8fafb;border:1px solid var(--line);border-radius:10px;padding:13px 16px;margin:12px 0}
summary{cursor:pointer;font-weight:650;color:#20343e}
details .details-body{margin-top:12px}
h3.subheading{font-size:1rem;margin:20px 0 10px}
.footer{color:var(--muted);border-top:1px solid var(--line);padding-top:18px;margin-top:28px;font-size:.84rem;display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap}
@media(max-width:860px){.layout{display:block}.sidebar{position:relative;height:auto}.wrap{padding:20px 14px 40px}.metric-grid{grid-template-columns:1fr 1fr}}
@media print{.sidebar{display:none!important}.layout{display:block}.wrap{max-width:none;padding:0}.report-section,.section-shell{break-inside:avoid;box-shadow:none}details .details-body{display:block!important}details{background:#fff;border:1px solid #ccc}}
"""

    js = """
const links=[...document.querySelectorAll('.sidebar a')];
const observer=new IntersectionObserver(entries=>entries.forEach(entry=>{if(entry.isIntersecting){links.forEach(link=>link.classList.toggle('active',link.getAttribute('href')==='#'+entry.target.id));}}),{rootMargin:'-20% 0px -65% 0px'});
document.querySelectorAll('.report-section').forEach(section=>observer.observe(section));
"""

    nav_links = "".join(
        f'<a href="#{sid}">{title}</a>'
        for sid, title in (
            ("section-01", "1 项目概览与核心结论"),
            ("section-02", "2 基础分析条目"),
            ("section-03", "3 高级分析条目"),
            ("section-04", "4 解释边界与复现信息"),
        )
    )

    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Celatlas Spatial SLAM 分析条目清单</title>
<style>{css}</style>
</head>
<body>
<div class="layout">
<nav class="sidebar"><h3>Analysis Catalog</h3>{nav_links}</nav>
<main class="wrap">
<section class="report-section" id="section-01">
  <div class="section-shell">
    <div class="hero">
      <h1>Celatlas Spatial SLAM 分析条目清单</h1>
      <p>项目：{esc(PROJECT)} · 样本：{esc(CSV.join(SAMPLES))} · 流程：{esc(STAGE_RANGE)}</p>
      <p>归类原则：基础分析 = 原始数据到标准空间结构的通用链路；高级分析 = 依赖 SLAM 代谢标记的定量与动力学。</p>
    </div>
    <div class="metric-grid">{metric_cards}</div>
    <div class="finding-grid">{finding_cards}</div>
  </div>
</section>
<section class="report-section" id="section-02">
  <div class="section-shell">
    <div class="section-heading"><div class="section-number">2</div><div><h2>基础分析条目</h2>
    <p class="section-lead">从原始数据到标准空间结构的通用空间转录组分析，共 9 条（B1–B9）。</p></div></div>
    {basic_table}
    <div class="notice info">基础分析不依赖代谢标记，任何空间转录组项目均可独立交付。</div>
  </div>
</section>
<section class="report-section" id="section-03">
  <div class="section-shell">
    <div class="section-heading"><div class="section-number">3</div><div><h2>高级分析条目</h2>
    <p class="section-lead">依赖 SLAM T&gt;C 标记的定量、动力学与模型诊断，共 8 条（A1–A8）。</p></div></div>
    {advanced_table}
    <div class="notice">高级分析均为单样本探索性结果；无 matched no-4sU 对照时，不输出绝对转录速率或半衰期。</div>
  </div>
</section>
<section class="report-section" id="section-04">
  <div class="section-shell">
    <div class="section-heading"><div class="section-number">4</div><div><h2>解释边界与复现信息</h2>
    <p class="section-lead">技术参数默认折叠，便于复核而不干扰主结果阅读。</p></div></div>
    <div class="notice">old/new RNA 只用于丰度分解、NTR 校正与 turnover，不单独聚类；空间域编号为各样本独立 Leiden 结果，跨样本需用 marker/表达程序重新对齐。tissue mask 是 RNA 表达密度推断区域，不是组织学染色边界；空间图统一 X 向右、Y 向下，UMAP 明确为非物理坐标。</div>
    <details><summary>分析参数</summary><div class="details-body">{parameter_table}</div></details>
    <details><summary>图件索引（07 与 08 阶段，PNG/PDF 各一份）</summary><div class="details-body">{figure_table}</div></details>
    <details><summary>复现与运行信息</summary><div class="details-body">{provenance_table}</div></details>
  </div>
</section>
<div class="footer"><span>Celatlas Spatial SLAM · 分析条目清单</span><span>Self-contained offline report · {esc(STAGE_RANGE)}</span></div>
</main>
</div>
<script>{js}</script>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("docs"))
    parser.add_argument("--md", action="store_true", help="write analysis_catalog.md")
    parser.add_argument("--html", action="store_true", help="write analysis_catalog.html")
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.md:
        (out_dir / "analysis_catalog.md").write_text(render_markdown(), encoding="utf-8")
        print(f"wrote {out_dir / 'analysis_catalog.md'}")
    if args.html:
        (out_dir / "analysis_catalog.html").write_text(render_html(), encoding="utf-8")
        print(f"wrote {out_dir / 'analysis_catalog.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
