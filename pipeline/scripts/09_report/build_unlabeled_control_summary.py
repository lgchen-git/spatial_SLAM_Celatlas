#!/usr/bin/env python3
"""Build analysis-only Markdown and a project-level unlabeled-control report."""
from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path


def read_json(path: Path) -> dict:
    if not path.is_file() or not path.stat().st_size:
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pct(value) -> str:
    try:
        return f"{float(value) * 100:.3f}%"
    except (TypeError, ValueError):
        return "NA"


def num(value, digits=4) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "NA"


def count(value) -> str:
    try:
        return f"{int(float(value)):,}"
    except (TypeError, ValueError):
        return "NA"


def collect_sample(results_root: Path, sample: str) -> dict:
    root = results_root / sample
    barcode = read_json(root / f"03_barcode/05_metadata/{sample}.barcode_qc.json")
    calls = read_json(root / f"05_molecule_call/05_metadata/{sample}.qc.json")
    correction = read_json(root / "07_spatial_dynamics/04_tool_output/correction/new_fraction_correction_summary.json")
    spatial = read_json(root / "07_spatial_dynamics/05_metadata/coordinate_and_figure_provenance.json")
    bulk = read_json(root / "09_report/01_tables/bulk_slam_summary.json")
    celatlas = read_json(root / f"02_celatlas/04_tool_output/celatlas_full/results" / next_project(root) / sample / ".metrics.json")
    report_validation = read_json(root / "09_report/05_metadata/report_validation.json")
    total_reads = barcode.get("total_reads")
    kept_reads = barcode.get("kept_reads")
    molecules = calls.get("molecules_total")
    hard_new = calls.get("molecules_new")
    params = correction.get("parameters", {})
    assignment = spatial.get("assignment_totals", {})
    total_expected = float(spatial.get("total_molecules", 0) or 0)
    apparent_new = float(assignment.get("new", 0) or 0)
    return {
        "sample": sample,
        "raw_reads": total_reads,
        "kept_reads": kept_reads,
        "kept_fraction": (float(kept_reads) / float(total_reads)) if total_reads and kept_reads is not None else None,
        "unique_mapping_fraction": celatlas.get("star_summary", {}).get("Uniquely Mapped Reads Fraction"),
        "molecules_total": molecules,
        "hard_new": hard_new,
        "hard_new_fraction": (float(hard_new) / float(molecules)) if molecules and hard_new is not None else None,
        "global_tc_rate": params.get("global_tc_rate"),
        "mixture_pi": params.get("pi"),
        "p_old": params.get("p_old"),
        "p_new": params.get("p_new"),
        "converged": params.get("converged"),
        "iterations": params.get("iterations"),
        "apparent_ntr": apparent_new / total_expected if total_expected else None,
        "n_bins": spatial.get("n_bins"),
        "n_domains": len(spatial.get("domains", [])),
        "slamdunk_weighted_tc": bulk.get("slamdunk", {}).get("weighted_conversion_rate"),
        "grand_median_ntr": bulk.get("grand_slam", {}).get("median_gene_mean_ntr"),
        "report_ok": report_validation.get("ok"),
    }


def next_project(root: Path) -> str:
    projects = sorted((root / "02_celatlas/04_tool_output/celatlas_full/results").glob("*"))
    projects = [path.name for path in projects if path.is_dir()]
    if len(projects) != 1:
        raise RuntimeError(f"Expected one reused Celatlas project under {root}, found {projects}")
    return projects[0]


def sample_markdown(row: dict) -> str:
    return f"""# {row['sample']} 普通空转 Spatial SLAM 阴性对照分析

## 分析定位

本样本未进行代谢标记。流程输出的 hard-new、expected new、表观 NTR 和 GRAND-SLAM NTR 均用于量化无标记背景与假阳性，不能解释为新生 RNA、转录速率或半衰期。

## 核心结果

| 指标 | 结果 |
|---|---:|
| Raw reads | {count(row['raw_reads'])} |
| Barcode/PolyT 保留比例 | {pct(row['kept_fraction'])} |
| Celatlas 唯一比对比例 | {num(row['unique_mapping_fraction'], 2)}% |
| 分子总数 | {count(row['molecules_total'])} |
| 至少一次 T>C 的背景阳性分子比例 | {pct(row['hard_new_fraction'])} |
| 全局 T>C / covered-T | {pct(row['global_tc_rate'])} |
| 两组分模型表观 new 先验 | {pct(row['mixture_pi'])} |
| 空间表观 NTR | {pct(row['apparent_ntr'])} |
| Total RNA bin50 / 空间域 | {count(row['n_bins'])} / {count(row['n_domains'])} |
| slamdunk 加权 UTR T>C | {pct(row['slamdunk_weighted_tc'])} |
| GRAND-SLAM 基因表观 NTR 中位数 | {pct(row['grand_median_ntr'])} |

## 判读

- Total RNA 空间域来自官方 Celatlas bin50 矩阵，可作为表达结构结果独立解释。
- T>C 与表观 NTR 反映无标记条件下的测序、比对、残余变异和模型分配背景。
- 模型收敛状态为 `{row['converged']}`，迭代次数为 `{row['iterations']}`；无论是否收敛，都不改变阴性对照的解释边界。
- 空间背景图为表达与分割派生图，不是 H&E 或其他组织学图像。
"""


def write_sample(results_root: Path, sample: str) -> Path:
    row = collect_sample(results_root, sample)
    output = results_root / sample / "09_report" / f"{sample}_analysis_results.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(sample_markdown(row), encoding="utf-8")
    return output


def summary_table(rows: list[dict]) -> str:
    body = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(value)}</td>"
            for value in (
                row["sample"], count(row["raw_reads"]), pct(row["kept_fraction"]),
                pct(row["hard_new_fraction"]), pct(row["global_tc_rate"]), pct(row["apparent_ntr"]),
                pct(row["slamdunk_weighted_tc"]), pct(row["grand_median_ntr"]),
                "通过" if row["report_ok"] else "未通过",
            )
        )
        + "</tr>"
        for row in rows
    )
    return f"""<div class="table-wrap"><table><thead><tr><th>样本</th><th>Raw reads</th><th>Barcode 保留</th><th>hard-new 背景</th><th>全局 T&gt;C</th><th>空间表观 NTR</th><th>slamdunk T&gt;C</th><th>GRAND 表观 NTR</th><th>样本报告</th></tr></thead><tbody>{body}</tbody></table></div>"""


def project_markdown(rows: list[dict]) -> str:
    lines = [
        "# 普通空转 Spatial SLAM 阴性对照汇总", "", "## 分析定位", "",
        "两份普通空间转录组样本均未进行代谢标记。本汇总比较流程在无标记数据上的背景表现，不进行生物学差异检验，也不作为非匹配标记样本的定量校准。", "",
        "## 核心指标", "",
        "| 样本 | Raw reads | Barcode 保留 | hard-new 背景 | 全局 T>C | 空间表观 NTR | slamdunk T>C | GRAND 表观 NTR |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample']} | {count(row['raw_reads'])} | {pct(row['kept_fraction'])} | {pct(row['hard_new_fraction'])} | "
            f"{pct(row['global_tc_rate'])} | {pct(row['apparent_ntr'])} | {pct(row['slamdunk_weighted_tc'])} | {pct(row['grand_median_ntr'])} |"
        )
    lines.extend([
        "", "## 解释", "",
        "- 两样本差异反映测序深度、组织组成、比对及技术背景的共同影响，不能视为处理效应。",
        "- 可把 hard-new、全局 T>C、表观 NTR 和 bulk 模型输出作为后续同平台实验的背景参考范围，但不是 matched no-4sU 校准参数。",
        "- Total RNA 空间结构可以独立使用；所有 SLAM 派生量均保持阴性对照解释。",
        "- 所有空间背景均为表达与分割派生图，不是组织学图像。", "",
    ])
    return "\n".join(lines)


def project_html(rows: list[dict]) -> str:
    table = summary_table(rows)
    return f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>普通空转 Spatial SLAM 阴性对照汇总</title><style>
:root{{--ink:#21323a;--muted:#62747c;--nav:#20343e;--line:#dbe3e6;--panel:#fff;--bg:#f3f6f7;--accent:#b85c38}}*{{box-sizing:border-box}}body{{margin:0;font-family:"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;color:var(--ink);background:var(--bg)}}.layout{{display:grid;grid-template-columns:250px 1fr;min-height:100vh}}.sidebar{{position:sticky;top:0;height:100vh;background:var(--nav);padding:28px 20px;color:white}}.sidebar a{{display:block;color:#dce6ea;text-decoration:none;padding:9px 0}}main{{padding:32px;max-width:1380px}}.hero,.card{{background:var(--panel);border:1px solid var(--line);border-radius:13px;padding:24px;margin-bottom:20px}}.hero{{background:linear-gradient(125deg,#20343e,#28776f);color:white}}h1{{margin:0 0 10px}}p,li{{line-height:1.7}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:.86rem}}th,td{{padding:10px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}th{{background:#edf3f4}}.notice{{border-left:4px solid var(--accent);background:#fff5e9;padding:14px 16px}}.footer{{color:var(--muted);border-top:1px solid var(--line);padding-top:18px;margin-top:28px}}@media(max-width:760px){{.layout{{display:block}}.sidebar{{height:auto;position:relative}}main{{padding:14px}}}}</style></head><body><div class="layout"><nav class="sidebar"><h3>Spatial SLAM Control</h3><a href="#overview">1 项目定位</a><a href="#metrics">2 背景指标</a><a href="#limits">3 解释边界</a></nav><main><section class="hero" id="overview"><h1>普通空转 Spatial SLAM 阴性对照汇总</h1><p>两份无代谢标记样本用于评估流程背景与假阳性，不用于真实新生 RNA 或动力学推断。</p></section><section class="card" id="metrics"><h2>背景指标</h2>{table}</section><section class="card" id="limits"><h2>解释边界</h2><div class="notice">hard-new、expected new、表观 NTR 和 GRAND-SLAM NTR 均为无标记背景诊断。Total RNA 空间结构可独立解释；SLAM 派生量不能用于半衰期、转录速率或标记效率推断。</div><ul><li>样本间仅做描述性背景比较，不作生物学差异检验。</li><li>本批不是任何标记样本的 matched no-4sU 对照。</li><li>空间背景为表达与分割派生图，不是 H&amp;E。</li></ul></section><div class="footer">Celatlas Spatial SLAM · Self-contained offline report</div></main></div></body></html>"""


def write_project(results_root: Path, samples: list[str]) -> Path:
    rows = [collect_sample(results_root, sample) for sample in samples]
    stage = results_root / "_project/09_report"
    tables, metadata = stage / "01_tables", stage / "05_metadata"
    for path in (tables, stage / "02_figures", stage / "03_logs", stage / "04_tool_output", metadata):
        path.mkdir(parents=True, exist_ok=True)
    table_path = tables / "unlabeled_control_summary.tsv"
    fields = list(rows[0]) if rows else ["sample"]
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    markdown = stage / "unlabeled_control_analysis_results.md"
    report = stage / "unlabeled_control_summary_report.html"
    markdown.write_text(project_markdown(rows), encoding="utf-8")
    report.write_text(project_html(rows), encoding="utf-8")
    manifest = metadata / "delivery_manifest.tsv"
    manifest_rows = [
        {"category": "01_tables", "description": "Unlabeled-control summary metrics", "source": str(table_path), "target": str(table_path.resolve()), "status": "up_to_date"},
        {"category": "09_report", "description": "Analysis-only Chinese Markdown summary", "source": str(markdown), "target": str(markdown.resolve()), "status": "up_to_date"},
        {"category": "09_report", "description": "Final offline unlabeled-control report", "source": str(report), "target": str(report.resolve()), "status": "up_to_date"},
    ]
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest_rows)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--mode", choices=("sample", "project"), required=True)
    parser.add_argument("--sample")
    args = parser.parse_args()
    rows = read_samples(args.samples)
    sample_ids = [row["sample_id"] for row in rows if row.get("labeling_status") == "unlabeled_control"]
    if args.mode == "sample":
        if not args.sample or args.sample not in sample_ids:
            raise SystemExit("--sample must name an unlabeled_control sample from the sample sheet")
        print(write_sample(args.results_root, args.sample))
    else:
        if not sample_ids:
            raise SystemExit("No unlabeled_control samples found")
        print(write_project(args.results_root, sample_ids))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
